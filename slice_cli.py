"""
Protocole d'appel du tranchage (apex/integration/RULES.md, §5 ; contrats C2 et C3).

    python main.py steps
    python main.py slice --request <requete.json>

`requete.json` :
    {
      "contract": 1,                      # APEX_CONTRACT_VERSION, obligatoire
      "mesh": "chemin/modele.stl",        # obligatoire
      "output_dir": "chemin/sortie",      # obligatoire, cree si absent
      "output_name": "modele.gcode",      # facultatif (defaut : <nom du maillage>.gcode)
      "params": { ... },                  # PipelineParams, partiel (schema C1)
      "gcode_params": { ... }             # GcodeParams, partiel (schema C1)
    }
Un champ inconnu, a n'importe quel niveau, est une ERREUR (code 2).

Sorties :
  - sortie d'erreur : une ligne JSON par evenement, `{"type":"step","slug":..,"fraction":..}` ou
    `{"type":"log","level":..,"message":..}` ; la progression vient UNIQUEMENT des etapes du
    moteur, jamais d'une minuterie. Le journal libre du pipeline est aussi envoye ici ;
  - sortie standard : un seul objet JSON, derniere ligne (le resultat).
Codes de sortie : 0 reussi ; 2 requete invalide ; 3 modele refuse ; 1 plantage.
Le G-code est ecrit sous un nom temporaire puis renomme : un G-code present est complet.
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import os
import sys
import traceback

import numpy as np

from config import schema
from core import progress as _core_progress

EXIT_OK, EXIT_ERROR, EXIT_INVALID, EXIT_REFUSED = 0, 1, 2, 3

_REQUEST_FIELDS = {"contract", "mesh", "output_dir", "output_name", "params", "gcode_params"}
_REQUIRED_FIELDS = ("contract", "mesh", "output_dir")


class InvalidRequest(ValueError):
    """Requete mal formee, version du contrat differente, champ inconnu ou valeur refusee."""


ALGORITHMS = ("element-ardor", "element-aeolis")


def _steps(algorithm: str = "element-ardor") -> list[dict]:
    """Etapes de progression d'un algorithme ; Element Ardor par defaut (pipeline principal)."""
    if algorithm == "element-aeolis":
        import importlib
        aeolis = importlib.import_module("nonplanar.slicer")
        return [{"slug": slug, "label": label} for slug, label in aeolis.AEOLIS_STEPS]
    from main import PIPELINE_STEPS
    return [{"slug": slug, "label": label} for slug, label in PIPELINE_STEPS] + \
        [{"slug": "export", "label": "Ecriture du G-code"}]


def _json_default(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (set, tuple)):
        return list(value)
    return str(value)


def _dumps(obj) -> str:
    return json.dumps(obj, default=_json_default, ensure_ascii=False, allow_nan=False)


def _clean(obj):
    """Remplace les flottants non finis (inf/nan) par None : le JSON strict les refuse."""
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (float, np.floating)) and not np.isfinite(obj):
        return None
    return obj


def _emit(stream, event: dict) -> None:
    print(_dumps(_clean(event)), file=stream, flush=True)


def steps_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="main.py steps",
                                     description="Liste des etapes affichees pendant le tranchage.")
    parser.add_argument("--algorithm", choices=ALGORITHMS, default="element-ardor",
                        help="Algorithme de tranchage (defaut : element-ardor).")
    args = parser.parse_args(argv)
    print(_dumps({"contract": schema.APEX_CONTRACT_VERSION, "algorithm": args.algorithm,
                  "steps": _steps(args.algorithm)}), flush=True)
    return EXIT_OK


def parse_request(text: str) -> dict:
    """Valide la requete (structure, version, champs) et rend ses objets : chemins, `PipelineParams`
    et `GcodeParams` complets. Ne touche ni au maillage ni au disque."""
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise InvalidRequest(f"requete illisible (JSON invalide) : {exc}") from exc
    if not isinstance(data, dict):
        raise InvalidRequest("la requete doit etre un objet JSON")
    unknown = sorted(set(data) - _REQUEST_FIELDS)
    if unknown:
        raise InvalidRequest(f"champ(s) inconnu(s) dans la requete : {unknown}")
    missing = [name for name in _REQUIRED_FIELDS if name not in data]
    if missing:
        raise InvalidRequest(f"champ(s) obligatoire(s) absent(s) : {missing}")
    contract = data["contract"]
    if isinstance(contract, bool) or contract != schema.APEX_CONTRACT_VERSION:
        raise InvalidRequest(f"version du contrat {contract!r} refusee : le moteur est a la version "
                             f"{schema.APEX_CONTRACT_VERSION}")
    for name in ("mesh", "output_dir") + (("output_name",) if "output_name" in data else ()):
        if not isinstance(data[name], str) or not data[name]:
            raise InvalidRequest(f"{name} : chemin attendu (chaine non vide)")
    name = data.get("output_name")
    if name is not None and (os.path.basename(name) != name or name in (".", "..")):
        raise InvalidRequest("output_name : un nom de fichier, sans dossier")

    from file_io.gcode_export import GcodeParams
    from file_io.pivot_length import current_pivot_to_tip
    from main import PipelineParams

    try:
        params = schema.from_jsonable(PipelineParams, data.get("params", {}), base=PipelineParams())
        gcode_base = GcodeParams(pivot_to_tip=current_pivot_to_tip()[0], strict=True)
        gcode_params = schema.from_jsonable(GcodeParams, data.get("gcode_params", {}), base=gcode_base)
    except (ValueError, TypeError, KeyError) as exc:
        raise InvalidRequest(str(exc)) from exc
    return {"mesh": data["mesh"], "output_dir": data["output_dir"], "output_name": name,
            "params": params, "gcode_params": gcode_params}


def export_refusal_message(report, gcode_params) -> str:
    """Pourquoi `strict` a refuse d'ecrire le G-code, avec les chiffres du rapport (accord en nombre)."""
    causes = []
    n = int(report.tip_below_plate)
    if n > 0:
        depth = max(0.0, -float(report.tip_z_min))
        causes.append(f"la pointe de la buse passerait sous le plateau pour {n} mouvement{'s' if n > 1 else ''}"
                      f" (jusqu'à {depth:.2f} mm sous le plateau)".replace(".", ","))
    limit = gcode_params.max_y_turns
    if limit is not None and report.y_turns > limit:
        causes.append(f"la rotation cumulée de l'axe Y ({report.y_turns:.1f} tours) dépasse la limite "
                      f"({limit:.1f} tours)".replace(".", ","))
    if not causes:
        causes.append("un garde-fou de l'export a refusé l'écriture")
    return "Le G-code n'a pas été écrit : " + " ; ".join(causes) + "."


def _result(status: str, **fields) -> dict:
    out = {"status": status, "contract": schema.APEX_CONTRACT_VERSION}
    out.update(fields)
    return out


def _counts(result) -> dict:
    return {
        "layers": len(result.layers),
        "toolpath_points": len(result.adjusted_toolpath.points),
        "step8_unreachable_orientations": len(result.step8_diagnostics),
        "collisions_detected": len(result.collision_report.events),
        "collisions_unresolved": len(result.adjustment_report.unresolved_collisions),
        "adjustments": len(result.adjustment_report.entries),
        "adjustment_converged": bool(result.adjustment_report.converged),
    }


def slice_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="main.py slice", description="Tranche un modele selon une requete JSON.")
    parser.add_argument("--request", required=True, help="Fichier requete.json.")
    args = parser.parse_args(argv)
    real_stdout, log = sys.stdout, sys.stderr

    def refuse(status: str, code: int, message: str, **extra) -> int:
        _emit(log, {"type": "log", "level": "error", "code": status, "message": message})
        print(_dumps(_clean(_result(status, message=message, **extra))), file=real_stdout, flush=True)
        return code

    try:
        with open(args.request, encoding="utf-8") as handle:
            request = parse_request(handle.read())
    except OSError as exc:
        return refuse("invalid_request", EXIT_INVALID, f"requete illisible : {exc}")
    except InvalidRequest as exc:
        return refuse("invalid_request", EXIT_INVALID, str(exc))
    if not os.path.isfile(request["mesh"]):
        return refuse("invalid_request", EXIT_INVALID, f"maillage introuvable : {request['mesh']}")

    from main import export_gcode_checked, run_pipeline
    from pipeline.step1_ingest import MeshIngestError, SurfaceInconsistentError

    os.makedirs(request["output_dir"], exist_ok=True)
    name = request["output_name"] or os.path.splitext(os.path.basename(request["mesh"]))[0] + ".gcode"
    final_path = os.path.join(request["output_dir"], name)
    temp_path = final_path + ".part"
    for stale in (temp_path,):
        with contextlib.suppress(FileNotFoundError):
            os.remove(stale)

    def progress(slug: str, fraction: float) -> None:
        _emit(log, {"type": "step", "slug": slug, "fraction": float(fraction)})

    try:
        # Le journal libre du pipeline (print) va sur la sortie d'erreur : la sortie standard ne
        # porte que le resultat.
        with contextlib.redirect_stdout(log):
            result = run_pipeline(request["mesh"], request["params"], progress=progress)
            progress("export", 0.0)
            _core_progress.begin(lambda fraction: progress("export", fraction))   # R-017
            try:
                _lines, report = export_gcode_checked(result.adjusted_toolpath, request["params"], temp_path,
                                                      gcode_params=request["gcode_params"])
            finally:
                _core_progress.end()
        if not report.written:
            with contextlib.suppress(FileNotFoundError):
                os.remove(temp_path)
            return refuse("export_refused", EXIT_REFUSED,
                          export_refusal_message(report, request["gcode_params"]),
                          report=dataclasses.asdict(report), **_counts(result))
        os.replace(temp_path, final_path)
        progress("export", 1.0)
    except MeshIngestError as exc:
        code = "surface_inconsistent" if isinstance(exc, SurfaceInconsistentError) else "model_refused"
        return refuse(code, EXIT_REFUSED, str(exc), diagnostic=getattr(exc, "diagnostic", None))
    except Exception as exc:  # plantage : code 1, trace sur la sortie d'erreur, jamais un G-code incomplet
        with contextlib.suppress(FileNotFoundError):
            os.remove(temp_path)
        traceback.print_exc(file=log)
        return refuse("error", EXIT_ERROR, f"{type(exc).__name__} : {exc}")
    finally:
        _core_progress.end()   # R-017 : jamais de rappel vers `log`, qui va etre ferme
        _core_progress.set_stage_switcher(None)

    warnings = []
    if result.adjustment_report.unresolved_collisions:
        warnings.append({"code": "collisions_unresolved",
                         "count": len(result.adjustment_report.unresolved_collisions)})
    ingest = result.ingest_report
    print(_dumps(_clean(_result(
        "ok", gcode=os.path.abspath(final_path), report=dataclasses.asdict(report),
        placement_offset_z=getattr(ingest, "placement_offset_z", None),
        platform_z=getattr(ingest, "platform_z_used", None),
        warnings=warnings, **_counts(result)))), file=real_stdout, flush=True)
    return EXIT_OK
