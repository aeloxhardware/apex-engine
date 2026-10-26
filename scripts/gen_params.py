"""
Genere, depuis les dataclasses Python (source de verite unique) :
  - config/params.schema.json   : schema JSON du contrat UI TypeScript -> moteur ;
  - cpp/include/apex/params.hpp : structures C++ equivalentes, valeurs par
                                  defaut comprises, pour le noyau natif.

Usage : .venv/Scripts/python.exe scripts/gen_params.py [--check]
--check : n'ecrit rien, sort en erreur si les fichiers ne sont pas a jour
(utilise par tests/test_config_schema.py).
"""
from __future__ import annotations

import dataclasses
import enum
import json
import os
import sys
import typing

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from config import schema  # noqa: E402
from config.params import ROOTS  # noqa: E402

SCHEMA_PATH = os.path.join(ROOT, "config", "params.schema.json")
HPP_PATH = os.path.join(ROOT, "cpp", "include", "apex", "params.hpp")


def build_schema() -> str:
    doc = {"$schema": schema.SCHEMA_DIALECT, "title": "Apex engine parameters",
           "description": "Genere par scripts/gen_params.py depuis les dataclasses Python. NE PAS EDITER.",
           "x-apex-contract": schema.APEX_CONTRACT_VERSION, "$defs": {}}
    roots = {}
    for name, cls in sorted(ROOTS.items()):
        s = schema.json_schema(cls)
        doc["$defs"].update(s.pop("$defs"))
        s.pop("$schema", None)
        doc["$defs"][name] = s
        roots[name] = {"$ref": f"#/$defs/{name}"}
    doc["type"] = "object"
    doc["properties"] = roots
    doc["additionalProperties"] = False
    return json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=False) + "\n"


# ---------------------------------------------------------------------------
# C++
# ---------------------------------------------------------------------------

def _cpp_literal(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        r = repr(float(v))
        return r if ("." in r or "e" in r or "inf" in r) else r + ".0"
    raise TypeError(v)


def _cpp_type(t, enums: dict, structs: list) -> str:
    opt, inner = schema._is_optional(t)
    if opt:
        return f"std::optional<{_cpp_type(inner, enums, structs)}>"
    if dataclasses.is_dataclass(t):
        if t not in structs:
            structs.append(t)
        return t.__name__
    if isinstance(t, type) and issubclass(t, enum.Enum):
        enums[t.__name__] = t
        return t.__name__
    if t is np.ndarray:
        return "std::vector<double>"
    origin, args = typing.get_origin(t), typing.get_args(t)
    if origin is list or t is list:
        return f"std::vector<{_cpp_type(args[0], enums, structs)}>" if args else "std::vector<double>"
    if origin is tuple or t is tuple:
        if args and args[-1] is not Ellipsis and all(a is float for a in args):
            return f"std::array<double, {len(args)}>"
        return "std::vector<double>"
    if origin is dict or t is dict:
        return "std::map<int, int>"
    return {float: "double", int: "int", bool: "bool", str: "std::string"}[t]


def cpp_value(value, t, enums, structs) -> str:
    """Expression C++ complete de `value` pour le type Python `t`."""
    opt, inner = schema._is_optional(t)
    if value is None:
        # champ non optionnel dont le defaut Python est None (renseigne a
        # l'execution, ex. MachineProfile.platform_origin) : valeur vide
        return "std::nullopt" if opt else ""
    t = inner
    ct = _cpp_type(t, enums, structs)
    if dataclasses.is_dataclass(t):
        hints = typing.get_type_hints(t)
        parts = []
        for f in dataclasses.fields(t):
            if schema.runtime_only(hints[f.name]):
                continue
            parts.append(f".{f.name} = {cpp_value(getattr(value, f.name), hints[f.name], enums, structs)}")
        return f"{ct}{{" + ", ".join(parts) + "}"
    if isinstance(value, enum.Enum):
        return f"{type(value).__name__}::{value.name}"
    origin, args = typing.get_origin(t), typing.get_args(t)
    if t is np.ndarray:
        return "std::vector<double>{" + ", ".join(_cpp_literal(float(x)) for x in np.ravel(value)) + "}"
    if origin is list or t is list:
        elem = args[0] if args else None
        if elem is None:
            elem = type(value[0]) if value and dataclasses.is_dataclass(value[0]) else float
        return f"{ct}{{" + ", ".join(cpp_value(v, elem, enums, structs) for v in value) + "}"
    if origin is tuple or t is tuple:
        return f"{ct}{{" + ", ".join(_cpp_literal(float(x)) for x in value) + "}"
    if origin is dict or t is dict:
        return f"{ct}{{" + ", ".join(f"{{{int(k)}, {int(v)}}}" for k, v in value.items()) + "}"
    if t is str:
        return json.dumps(value)
    if t is float:
        return _cpp_literal(float(value))
    if t is int:
        return _cpp_literal(int(value))
    return _cpp_literal(value)


def build_hpp() -> str:
    enums: dict = {}
    structs: list = []
    for cls in ROOTS.values():
        structs.append(cls)
    # parcours pour decouvrir toutes les structures et enums
    i = 0
    while i < len(structs):
        cls = structs[i]
        hints = typing.get_type_hints(cls)
        for f in dataclasses.fields(cls):
            if not schema.runtime_only(hints[f.name]):
                _cpp_type(hints[f.name], enums, structs)
        i += 1
    # ordre : dependances d'abord
    ordered, placed = [], set()

    def place(c):
        if c in placed:
            return
        hints = typing.get_type_hints(c)
        for f in dataclasses.fields(c):
            _, inner = schema._is_optional(hints[f.name])
            for cand in [inner, *typing.get_args(inner)]:
                if dataclasses.is_dataclass(cand) and cand is not c:
                    place(cand)
        placed.add(c)
        ordered.append(c)

    for c in structs:
        place(c)

    out = ["// GENERE par scripts/gen_params.py depuis les dataclasses Python. NE PAS EDITER.",
           "// Contrat de parametres du moteur (meme contenu que config/params.schema.json).",
           "// Les champs Python `typing.Any` (ex. MachineProfile.nozzle_model) sont",
           "// volontairement absents : donnees d'execution, pas parametres.",
           "#pragma once", "", "#include <array>", "#include <map>", "#include <optional>",
           "#include <string>", "#include <vector>", "", "namespace apex::params {", ""]
    for name in sorted(enums):
        e = enums[name]
        out.append(f"enum class {name} {{ " + ", ".join(m.name for m in e) + " };")
        out.append(f"inline const char* to_string({name} v) {{")
        out.append("    switch (v) {")
        for m in e:
            out.append(f"        case {name}::{m.name}: return {json.dumps(m.value)};")
        out.append("    }")
        out.append('    return "";')
        out.append("}")
        out.append("")
    for cls in ordered:
        hints = typing.get_type_hints(cls)
        try:
            inst = cls()
        except TypeError:
            inst = None
        out.append(f"struct {cls.__name__} {{")
        for f in dataclasses.fields(cls):
            t = hints[f.name]
            if schema.runtime_only(t):
                out.append(f"    // {f.name} : donnee d'execution, non transcrite")
                continue
            ct = _cpp_type(t, enums, structs)
            dv = cpp_value(getattr(inst, f.name), t, enums, structs) if inst is not None else ""
            out.append(f"    {ct} {f.name}{{{dv}}};")
        out.append("};")
        out.append("")
    out.append("}  // namespace apex::params")
    return "\n".join(out) + "\n"


def main():
    check = "--check" in sys.argv
    stale = []
    for path, content in ((SCHEMA_PATH, build_schema()), (HPP_PATH, build_hpp())):
        cur = open(path, encoding="utf-8").read() if os.path.exists(path) else None
        if cur != content:
            stale.append(path)
            if not check:
                with open(path, "w", encoding="utf-8", newline="\n") as f:
                    f.write(content)
    if check and stale:
        print("pas a jour : " + ", ".join(stale))
        sys.exit(1)
    print("a jour" if not stale else "ecrit : " + ", ".join(stale))


if __name__ == "__main__":
    main()
