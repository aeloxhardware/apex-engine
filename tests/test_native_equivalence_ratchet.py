"""R-030 point 6 : pas de correctif natif sans test d'equivalence Python/C++ (cliquet).

Toute fonction ou classe exportee par `apex_core` doit etre NOMMEE dans au moins un test. Les 20 noms ci-dessous
(audit du 2026-11-06) ne le sont pas : certains sont couverts indirectement par leur enveloppe Python
(`step4_field_solver`, `step8_ik`, ...), ce qui n'a PAS ete verifie un a un. La liste ne peut que RETRECIR : un
nouveau nom natif sans test fait echouer ce test ; un nom nomme dans un test doit etre retire d'ici.
"""
import re
from pathlib import Path

import pytest

from core import _cpp

TESTS = Path(__file__).resolve().parent

UNNAMED_BASELINE = {
    "FemSparsity", "PrintedGrid", "SkinCheckIndex", "apply_per_point_normals", "assemble_gradient_rhs",
    "blend_tet_directions", "build_raw_points", "clamp_head_height_apply", "fast_penetrations_batch", "global_stats",
    "heat_method_distance", "interpolate_out_of_plane", "median_filter_free_axes", "reset_global_stats",
    "resolve_axis_angles_batch", "solve_linear_system", "stable_topological_order", "tet_face_topology", "tet_volumes",
    "try_local_reorientation_coupled",
}


def _unnamed() -> set:
    module = _cpp.module()
    names = {n for n in dir(module) if not n.startswith("_") and callable(getattr(module, n))}
    text = "".join(p.read_text(encoding="utf-8", errors="replace")
                   for p in TESTS.glob("*.py") if p.name != Path(__file__).name)
    return {n for n in names if not re.search(r"\b" + re.escape(n) + r"\b", text)}


@pytest.mark.skipif(not _cpp.available(), reason="noyau natif non construit")
def test_every_native_function_is_named_in_a_test_or_is_in_the_shrinking_baseline():
    fresh = _unnamed() - UNNAMED_BASELINE
    assert not fresh, f"fonctions natives sans test d'equivalence : {sorted(fresh)}"


@pytest.mark.skipif(not _cpp.available(), reason="noyau natif non construit")
def test_baseline_only_lists_names_that_are_still_untested():
    stale = UNNAMED_BASELINE - _unnamed()
    assert not stale, f"a retirer de UNNAMED_BASELINE (desormais testes) : {sorted(stale)}"
