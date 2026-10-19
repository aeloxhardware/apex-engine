"""Champ multibranche (step4b) : opt-in, valeurs finies, continuite, arbre de branches."""
import os
from dataclasses import replace

import numpy as np

from main import PipelineParams
from pipeline import step4b_multibranch as mb
from pipeline.step1_ingest import run as ingest_run
from pipeline.step2_overhang import run as overhang_run
from pipeline.step3_segmentation import run as segmentation_run
from pipeline.step4_field_solver import FieldSolverParams, run as field_run

MODEL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "res", "test_models", "real_overhang.stl")


def _prepare():
    params = PipelineParams(resolution=1.5)
    md, tm, rep = ingest_run(MODEL, params.resolution, place_on_platform=True)
    om = overhang_run(md, params.overhang)
    seg = replace(params.segmentation, platform_z=rep.platform_z_used)
    return md, tm, om, segmentation_run(md, om, seg, tm), params


def test_multibranch_is_opt_in_and_default_unchanged():
    assert FieldSolverParams().multibranch is False


def test_multibranch_field_is_finite_and_branches_form_a_tree():
    md, tm, om, rg, params = _prepare()
    fs = replace(params.field_solver, multibranch=True, enforce_upward_growth=False,
                 multibranch_direction_mode="up", multibranch_partition="flow",
                 multibranch_min_branch_tets=10)
    field = field_run(tm, md, om, rg, None, fs)
    assert field.values.shape == (tm.node_count,)
    assert np.isfinite(field.values).all() and np.isfinite(field.gradient).all()
    diag = mb.LAST_DIAGNOSTIC
    assert diag.branch_count >= 2 and diag.root_count >= 1
    assert field.branch_labels.shape == (tm.tet_count,)
    # une branche sans racine plateau a un parent : profondeur finie partout
    assert all(d >= 0 for d in field.branch_depths.values())
    # les racines sont ancrees au plateau (valeur d'ancrage = 0)
    assert abs(field.values.min()) < 1e-6 or field.values.min() >= -1e-6
