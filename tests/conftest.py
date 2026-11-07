"""Position du noyau natif pour toute la suite : `pytest --native` (R-030 point 3, plus de variable d'environnement).

Sans l'option, la suite tourne sur Python ; avec, les `PipelineParams` construits pendant les tests heritent de
`RuntimeOptions(native=True)`. Les tests d'equivalence forcent eux-memes les deux positions (`_cpp.set_enabled`).
"""
import pytest

from core import runtime


def pytest_addoption(parser):
    parser.addoption("--native", action="store_true", default=False,
                     help="active le noyau natif (apex_core) pour toute la suite")


def pytest_configure(config):
    if config.getoption("--native"):
        runtime.replace_options(native=True)
