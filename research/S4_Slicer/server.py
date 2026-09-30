"""
Serveur interactif local pour piloter la methode S4-Slicer iteration par
iteration, sans re-tetraedriser a chaque fois (le tet courant reste en
memoire process).

Aucune dependance ajoutee : http.server natif de la stdlib. L'UI vit dans
ui.html (charge tel quel), qui parle a ce serveur en JSON via fetch().

Endpoints :
  GET  /                  -> sert ui.html
  GET  /state             -> etat courant (params, iteration, en cours ou non)
  POST /iterate           -> lance UNE iteration avec les params du corps JSON, en thread
  POST /reset             -> revient au maillage initial (re-tetraedrise une fois)
  GET  /model.stl         -> le maillage courant (undeformed si iteration=0, deforme sinon)
  GET  /progress          -> etape courante de l'iteration en cours (poll cote client)

Usage : .venv/Scripts/python.exe server.py [model_name] [port] [max_tet_volume]

Chaque iteration est aussi sauvegardee (surface STL + positions .npy) sous
output_models/ui_<modele>[_a<volume>]_pass<k>_ovh<MAX_OVERHANG>, au format
que relit analyze_series.py. Reset efface la serie precedente.
"""
import glob
import io
import os
import sys
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

import bird_core

MODEL_NAME = sys.argv[1] if len(sys.argv) > 1 else "real_overhang"
PORT = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
# Densite du maillage : sans elle, real_overhang n'a que 48 tetraedres, et
# tout repli obtenu est un artefact du maillage grossier (mesure le 2026-09-26).
MAX_TET_VOLUME = float(sys.argv[3]) if len(sys.argv) > 3 else None
SERIES_PREFIX = "output_models/ui_" + MODEL_NAME.replace(" ", "_") + (
    f"_a{MAX_TET_VOLUME:g}" if MAX_TET_VOLUME else "")

DEFAULT_PARAMS = {
    "NEIGHBOUR_LOSS_WEIGHT": 20,
    "MAX_OVERHANG": 30,
    "ROTATION_MULTIPLIER": 1,
    "SET_INITIAL_ROTATION_TO_ZERO": False,
    "INITIAL_ROTATION_FIELD_SMOOTHING": 30,
    "MAX_POS_ROTATION_DEG": 360,
    "MAX_NEG_ROTATION_DEG": -360,
    "ROTATE_ITERATIONS": 100,
    "DEFORM_ITERATIONS": 1000,
    "STEEP_OVERHANG_COMPENSATION": True,
}

# ---------------------------------------------------------------------------
# Etat global du process -- un seul client attendu (usage local interactif),
# proteger juste contre un double-clic concurrent.
# ---------------------------------------------------------------------------

_lock = threading.Lock()
_state = {
    "tet": None,            # UnstructuredGrid courant (deforme ou initial)
    "model_state": None,    # bird_core.ModelState (topologie, invariant)
    "initial_tet": None,    # copie du tet initial, pour Reset
    "iteration": 0,
    "busy": False,
    "progress_step": "",    # etape courante, mise a jour PENDANT le calcul
    "progress_pct": None,   # nfev/max_nfev * 100 pour l'etape en cours, ou None si non applicable
    "last_info": None,      # dict rotate_cost/deform_cost/overhang_fraction
    "history": [],          # liste de last_info, une entree par iteration
    "error": None,
}


def _ingest_model():
    _state["progress_step"] = "ingestion + tetraedrisation"
    tet = bird_core.load_and_tetrahedralize(f"input_models/{MODEL_NAME}.stl", max_tet_volume=MAX_TET_VOLUME)
    _state["progress_step"] = "calcul topologie (voisinages, cell_to_face)"
    tet, model_state = bird_core.calculate_tet_attributes(tet)
    _state["tet"] = tet
    _state["model_state"] = model_state
    _state["initial_tet"] = tet.copy(deep=True)
    for old in glob.glob(f"{SERIES_PREFIX}_pass*"):
        os.remove(old)
    np.save(f"{SERIES_PREFIX}_cells.npy", np.asarray(tet.field_data["cells"]))
    np.save(f"{SERIES_PREFIX}_pass0.npy", np.asarray(tet.points))
    _state["iteration"] = 0
    _state["last_info"] = None
    _state["history"] = []
    _state["progress_step"] = ""


_PHASE_LABELS = {
    "rotation": "optimisation des rotations (least_squares)",
    "deformation": "recollement du maillage (least_squares)",
}


def _on_progress(phase, nfev, max_nfev, cost):
    # Appele par le callback least_squares -- une fois par iteration EXTERNE
    # de trf, jamais a l'interieur d'une evaluation de fonction/jacobienne.
    # Ecriture directe sans lock : ce sont de simples entiers/floats/strings,
    # une lecture concurrente legerement en retard n'est jamais incoherente.
    pct = round(100 * nfev / max_nfev) if max_nfev else None
    _state["progress_step"] = f"{_PHASE_LABELS.get(phase, phase)} ({nfev}/{max_nfev} evaluations)"
    _state["progress_pct"] = pct


def _do_iteration(params):
    try:
        _state["busy"] = True
        _state["error"] = None
        _state["progress_step"] = "optimisation des rotations (least_squares)"
        _state["progress_pct"] = 0
        tet = _state["tet"]
        model_state = _state["model_state"]

        t0 = time.time()
        deformed_tet, info = bird_core.run_one_iteration(tet, model_state, params, progress_callback=_on_progress)
        info["elapsed_s"] = round(time.time() - t0, 2)
        # Garder les parametres de chaque iteration : sans eux, un repli
        # obtenu dans l'interface ne peut pas etre reproduit.
        info["params"] = dict(params)

        _state["progress_step"] = "export STL"
        _state["progress_pct"] = None
        export_tet = deformed_tet.copy(deep=True)
        bird_core.export_surface_stl(export_tet, "output_models/_server_current.stl")

        _state["tet"] = deformed_tet
        _state["iteration"] += 1
        k, ovh = _state["iteration"], params["MAX_OVERHANG"]
        bird_core.export_surface_stl(deformed_tet.copy(deep=True), f"{SERIES_PREFIX}_pass{k}_ovh{ovh:g}.stl")
        np.save(f"{SERIES_PREFIX}_pass{k}.npy", np.asarray(deformed_tet.points))
        _state["last_info"] = info
        _state["history"].append(info)
    except Exception as exc:  # noqa: BLE001 - on veut remonter n'importe quelle erreur a l'UI
        _state["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        _state["progress_step"] = ""
        _state["progress_pct"] = None
        _state["busy"] = False


def _do_reset():
    try:
        _state["busy"] = True
        _state["error"] = None
        _ingest_model()
        export_tet = _state["tet"].copy(deep=True)
        bird_core.export_surface_stl(export_tet, "output_models/_server_current.stl")
    except Exception as exc:  # noqa: BLE001
        _state["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        _state["progress_step"] = ""
        _state["busy"] = False


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):  # noqa: A002 - signature imposee par BaseHTTPRequestHandler
        pass  # silence l'access log par defaut, la console reste lisible

    def _send_json(self, payload, status=200):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/" or self.path == "/ui.html":
            with open("ui.html", "rb") as f:
                body = f.read()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        if self.path == "/state":
            with _lock:
                self._send_json({
                    "model_name": MODEL_NAME,
                    "iteration": _state["iteration"],
                    "busy": _state["busy"],
                    "progress_step": _state["progress_step"],
                    "last_info": _state["last_info"],
                    "history": _state["history"],
                    "error": _state["error"],
                    "number_of_cells": _state["tet"].number_of_cells if _state["tet"] is not None else None,
                    "default_params": DEFAULT_PARAMS,
                })
            return

        if self.path == "/progress":
            with _lock:
                self._send_json({
                    "busy": _state["busy"],
                    "progress_step": _state["progress_step"],
                    "progress_pct": _state["progress_pct"],
                })
            return

        if self.path == "/model.stl":
            with _lock:
                if _state["tet"] is None:
                    self.send_response(404)
                    self.end_headers()
                    return
            try:
                with open("output_models/_server_current.stl", "rb") as f:
                    body = f.read()
                self.send_response(200)
                self.send_header("Content-Type", "model/stl")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Cache-Control", "no-store")
                self.end_headers()
                self.wfile.write(body)
            except FileNotFoundError:
                self.send_response(404)
                self.end_headers()
            return

        self.send_response(404)
        self.end_headers()

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw) if raw else {}
        except json.JSONDecodeError:
            payload = {}

        if self.path == "/iterate":
            with _lock:
                if _state["busy"]:
                    self._send_json({"ok": False, "error": "une iteration est deja en cours"}, status=409)
                    return
                if _state["tet"] is None:
                    self._send_json({"ok": False, "error": "aucun modele charge, appeler /reset d'abord"}, status=409)
                    return
                params = {**DEFAULT_PARAMS, **payload}
                _state["busy"] = True  # pose immediatement le verrou avant de lancer le thread
            threading.Thread(target=_do_iteration, args=(params,), daemon=True).start()
            self._send_json({"ok": True})
            return

        if self.path == "/reset":
            with _lock:
                if _state["busy"]:
                    self._send_json({"ok": False, "error": "un calcul est deja en cours"}, status=409)
                    return
                _state["busy"] = True
            threading.Thread(target=_do_reset, daemon=True).start()
            self._send_json({"ok": True})
            return

        self.send_response(404)
        self.end_headers()


def main():
    print(f"[server] modele: {MODEL_NAME}")
    print("[server] ingestion initiale...")
    _do_reset()
    if _state["error"]:
        print(f"[server] ERREUR a l'ingestion: {_state['error']}")
        return
    print(f"[server] pret -- {_state['tet'].number_of_cells} tetraedres")
    print(f"[server] ouvrir http://127.0.0.1:{PORT}/ dans le navigateur")

    httpd = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
