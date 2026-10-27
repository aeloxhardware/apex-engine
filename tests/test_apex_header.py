"""R-010 : le chiffrement de l'en-tete Apex rejoue les vecteurs de l'interface (C5)."""
import json
import os

import pytest

from file_io import apex_header as ah

pytestmark = pytest.mark.skipif(not ah.cipher_available(), reason="paquet `cryptography` absent")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _vectors():
    with open(os.path.join(ROOT, "tests", "vectors", "apex_cipher_vectors.json"), encoding="utf-8") as f:
        return json.load(f)["vectors"]


def test_every_stage_matches_the_ui_vectors():
    vectors = _vectors()
    assert vectors
    for v in vectors:
        trace = {}
        out = ah.apex_encrypt(v["plain"], bytes.fromhex(v["nonce1Hex"]), bytes.fromhex(v["nonce2Hex"]), trace)
        assert trace["greatCipher"] == v["greatCipher"]
        # deflate n'est pas normatif octet pour octet : on compare a partir de l'etage 2 de l'interface
        assert ah._vigenere(ah._bytes_to_letters(bytes.fromhex(v["afterAes1Hex"])), 1) == v["vigenere"]
        padded = v["vigenere"] + "X" * ((ah.BLOCK - len(v["vigenere"]) % ah.BLOCK) % ah.BLOCK)
        assert ah._length_letters(len(v["vigenere"])) + ah._fleissner_encode(padded) == v["fleissner"]
        assert ah._enigma(v["fleissner"]) == v["enigma"]
        assert ah.apex_decrypt(out) == v["plain"]


def test_decrypts_what_the_ui_produced():
    for v in _vectors():
        assert ah.apex_decrypt(v["output"]) == v["plain"]


def test_block_roundtrip_and_reproducible():
    meta = {"generator": "Apex engine", "pivot_to_tip": 46.0, "note": "é à ü  "}
    a, b = ah.header_block_lines(meta), ah.header_block_lines(meta)
    assert a == b and a[0] == ah.BEGIN and a[-1] == ah.END and all(len(ln) <= len(ah.DATA) + 100 for ln in a[1:-1])
    got = ah.read_apex_header(["G21"] + a + ["G00 X1"])
    assert got == {"format": 1, "sliced": True, "meta": meta}
    assert ah.read_apex_header(["G21"]) is None
