"""
Socle Vulkan (phase 1 du plan « portage C++, puis Vulkan », 2026-10-13).

- l'absence de GPU (module construit sans Vulkan, pas de pilote) n'est jamais une
  erreur : `gpu_status()` rend ok = False et une RAISON ;
- quand un GPU est la, le noyau de fumee `saxpy` (y <- a*x + y, sans contraction
  en FMA) rend un resultat IDENTIQUE AU BIT au calcul CPU, en float32 et en float64,
  y compris au-dela d'un dispatch (decoupage).
"""
import numpy as np
import pytest

from core import _cpp


def test_gpu_status_always_has_a_reason_when_unavailable():
    st = _cpp.gpu_status()
    assert "ok" in st
    if not st["ok"]:
        assert st["reason"]
    assert _cpp.gpu_enabled() is False or st["ok"]


def test_gpu_switch_is_off_by_default(monkeypatch):
    monkeypatch.delenv("APEX_USE_GPU", raising=False)
    previous = _cpp.set_gpu_enabled(None)
    try:
        assert _cpp.gpu_enabled() is False
    finally:
        _cpp.set_gpu_enabled(previous)


def _require_gpu():
    st = _cpp.gpu_status()
    if not st.get("ok"):
        pytest.skip(f"GPU indisponible : {st.get('reason')}")
    return st


@pytest.mark.parametrize("n", [1, 1000, 4096 * 256 + 17])     # le dernier depasse un dispatch
def test_saxpy_f32_bit_exact(n):
    _require_gpu()
    rng = np.random.default_rng(n)
    x = rng.normal(size=n).astype(np.float32)
    y = rng.normal(size=n).astype(np.float32)
    a = np.float32(1.7)
    out, gpu_s = _cpp.module().gpu_saxpy(float(a), x, y)
    expected = (a * x).astype(np.float32) + y
    assert out.dtype == np.float32
    assert np.array_equal(out, expected)
    assert gpu_s != 0.0


def test_saxpy_f64_bit_exact():
    st = _require_gpu()
    if not st.get("fp64"):
        pytest.skip("shaderFloat64 indisponible")
    rng = np.random.default_rng(7)
    x = rng.normal(size=50_000)
    y = rng.normal(size=50_000)
    out, _ = _cpp.module().gpu_saxpy(0.3, x, y)
    assert np.array_equal(out, 0.3 * x + y)
