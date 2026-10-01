"""
(De)serialisation JSON des parametres du moteur, et schema JSON — pont avec
l'interface TypeScript (contrat UI -> moteur).

Tout est DERIVE des annotations de type des dataclasses (`typing.get_type_hints`),
jamais d'une liste tenue a la main : un champ ajoute a un `*Params` apparait
automatiquement dans le JSON et dans le schema, et le test
`test_every_field_in_schema` le verifie. Meme discipline que le reste du
projet : une seule source de verite.

Conventions :
- `np.ndarray` -> liste de nombres ; `tuple[...]` -> liste de longueur fixe ;
- `Enum` -> sa valeur (chaine) ;
- `dict[int, int]` -> objet JSON a cles chaines ("12": 3) ;
- `Optional[X]` -> X ou null ;
- champs `typing.Any` (ex. `MachineProfile.nozzle_model`, un maillage charge
  a l'execution) : NON serialises, marques `x-runtime-only` dans le schema.
  Le maillage de tete se designe par un chemin cote UI, charge par
  `file_io.mesh_loader.load_toolhead_mesh`.
"""
from __future__ import annotations

import dataclasses
import enum
import json
import types
import typing
from typing import Any

import numpy as np

SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"


def _is_optional(t) -> tuple[bool, Any]:
    origin = typing.get_origin(t)
    if origin in (typing.Union, types.UnionType):
        args = [a for a in typing.get_args(t) if a is not type(None)]
        if len(args) == 1 and len(typing.get_args(t)) == 2:
            return True, args[0]
    return False, t


def _hints(cls) -> dict:
    return typing.get_type_hints(cls)


def runtime_only(t) -> bool:
    return t is Any


# ---------------------------------------------------------------------------
# Dataclass -> JSON
# ---------------------------------------------------------------------------

def to_jsonable(value, t=None):
    """Valeur Python (dataclass, tableau, enum...) -> structure JSON."""
    if value is None:
        return None
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        hints = _hints(type(value))
        out = {}
        for f in dataclasses.fields(value):
            ft = hints[f.name]
            if runtime_only(ft):
                continue
            out[f.name] = to_jsonable(getattr(value, f.name), ft)
        return out
    if isinstance(value, enum.Enum):
        return value.value
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError(f"valeur non finie non representable en JSON : {value}")
    return value


def dumps(value, **kw) -> str:
    return json.dumps(to_jsonable(value), **kw)


# ---------------------------------------------------------------------------
# JSON -> dataclass
# ---------------------------------------------------------------------------

def _from(t, data, default=dataclasses.MISSING):
    opt, inner = _is_optional(t)
    if data is None:
        if opt or t is Any:
            return None
        raise ValueError(f"null interdit pour le type {t}")
    t = inner
    if t is Any:
        return data
    if dataclasses.is_dataclass(t):
        return from_jsonable(t, data)
    if isinstance(t, type) and issubclass(t, enum.Enum):
        return t(data)
    if t is np.ndarray:
        return np.asarray(data, dtype=float)
    origin = typing.get_origin(t)
    args = typing.get_args(t)
    if origin is list or t is list:
        elem = args[0] if args else None
        if elem is None:
            return [_guess_list_elem(v) for v in data]
        return [_from(elem, v) for v in data]
    if origin is tuple or t is tuple:
        if args and args[-1] is not Ellipsis:
            return tuple(_from(a, v) for a, v in zip(args, data))
        return tuple(data)
    if origin is dict or t is dict:
        kt, vt = args if args else (str, Any)
        return {_from(kt, k if kt is str else kt(k)): _from(vt, v) for k, v in data.items()}
    if t is float:
        return float(data)
    if t is int:
        if isinstance(data, float) and not data.is_integer():
            raise ValueError(f"entier attendu, recu {data}")
        return int(data)
    if t is bool:
        if not isinstance(data, bool):
            raise ValueError(f"booleen attendu, recu {data!r}")
        return data
    if t is str:
        return str(data)
    return data


def _guess_list_elem(v):
    """Liste non typee (ex. `LayerParams.support_rotation_axes: list | None`) :
    un objet {axis, angle_range_deg} est un RotationAxisSpec, le reste passe tel quel."""
    if isinstance(v, dict) and set(v) == {"axis", "angle_range_deg"}:
        from core.mesh_types import RotationAxisSpec
        return from_jsonable(RotationAxisSpec, v)
    return v


def from_jsonable(cls, data: dict, base=None):
    """Structure JSON -> instance de `cls`. Les champs absents gardent la
    valeur par defaut de la dataclass (ou celle de `base` si fournie), ce qui
    permet a l'UI de n'envoyer que ce qu'elle change. Un champ INCONNU est une
    erreur : une faute de frappe cote UI ne doit pas etre ignoree en silence."""
    if not isinstance(data, dict):
        raise ValueError(f"{cls.__name__} : objet JSON attendu, recu {type(data).__name__}")
    hints = _hints(cls)
    names = {f.name for f in dataclasses.fields(cls)}
    unknown = set(data) - names
    if unknown:
        raise ValueError(f"{cls.__name__} : champ(s) inconnu(s) {sorted(unknown)}")
    if base is None:
        required = [f.name for f in dataclasses.fields(cls)
                    if f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING]
        missing = [n for n in required if n not in data]
        if missing:
            raise ValueError(f"{cls.__name__} : champ(s) obligatoire(s) absent(s) {missing}")
        return cls(**{f.name: _from(hints[f.name], data[f.name]) for f in dataclasses.fields(cls)
                      if f.name in data})
    obj = base
    for f in dataclasses.fields(cls):
        if f.name not in data:
            continue
        ft = hints[f.name]
        cur = getattr(obj, f.name)
        if dataclasses.is_dataclass(_is_optional(ft)[1]) and isinstance(data[f.name], dict) and cur is not None:
            setattr(obj, f.name, from_jsonable(type(cur), data[f.name], base=cur))
        else:
            setattr(obj, f.name, _from(ft, data[f.name]))
    return obj


def loads(cls, text: str, base=None):
    return from_jsonable(cls, json.loads(text), base=base)


# ---------------------------------------------------------------------------
# Schema JSON
# ---------------------------------------------------------------------------

def _type_schema(t, defs: dict) -> dict:
    opt, inner = _is_optional(t)
    if opt:
        return {"anyOf": [_type_schema(inner, defs), {"type": "null"}]}
    if t is Any:
        return {"x-runtime-only": True}
    if dataclasses.is_dataclass(t):
        name = t.__name__
        if name not in defs:
            defs[name] = {}  # reserve (recursion)
            defs[name] = _class_schema(t, defs)
        return {"$ref": f"#/$defs/{name}"}
    if isinstance(t, type) and issubclass(t, enum.Enum):
        return {"type": "string", "enum": [m.value for m in t]}
    if t is np.ndarray:
        return {"type": "array", "items": {"type": "number"}}
    origin = typing.get_origin(t)
    args = typing.get_args(t)
    if origin is list or t is list:
        return {"type": "array", "items": _type_schema(args[0], defs) if args else {}}
    if origin is tuple or t is tuple:
        if args and args[-1] is not Ellipsis:
            items = [_type_schema(a, defs) for a in args]
            return {"type": "array", "prefixItems": items, "minItems": len(items), "maxItems": len(items)}
        return {"type": "array"}
    if origin is dict or t is dict:
        return {"type": "object", "additionalProperties": _type_schema(args[1], defs) if args else {}}
    return {float: {"type": "number"}, int: {"type": "integer"}, bool: {"type": "boolean"},
            str: {"type": "string"}}.get(t, {})


def _class_schema(cls, defs: dict) -> dict:
    hints = _hints(cls)
    props = {}
    try:
        inst = cls()
    except TypeError:
        inst = None
    for f in dataclasses.fields(cls):
        s = dict(_type_schema(hints[f.name], defs))
        if inst is not None and not runtime_only(hints[f.name]):
            try:
                s["default"] = to_jsonable(getattr(inst, f.name))
            except ValueError:
                pass
        props[f.name] = s
    required = [f.name for f in dataclasses.fields(cls)
                if f.default is dataclasses.MISSING and f.default_factory is dataclasses.MISSING]
    out = {"type": "object", "properties": props, "additionalProperties": False}
    if (cls.__doc__ or "").strip():
        out["description"] = " ".join((cls.__doc__ or "").split())[:500]
    if required:
        out["required"] = required
    return out


def json_schema(root_cls, title: str | None = None) -> dict:
    defs: dict = {}
    root = _class_schema(root_cls, defs)
    root["$schema"] = SCHEMA_DIALECT
    root["title"] = title or root_cls.__name__
    root["$defs"] = defs
    return root


def walk_fields(root_cls):
    """(chemin, type) de chaque champ de l'arbre de dataclasses — pour les tests
    de couverture et le generateur d'en-tete C++."""
    seen = set()

    def rec(cls, prefix):
        hints = _hints(cls)
        for f in dataclasses.fields(cls):
            t = hints[f.name]
            yield f"{prefix}{f.name}", t
            _, inner = _is_optional(t)
            cands = [inner] + list(typing.get_args(inner))
            for c in cands:
                if dataclasses.is_dataclass(c) and c not in seen:
                    seen.add(c)
                    yield from rec(c, f"{prefix}{f.name}.")
    yield from rec(root_cls, "")
