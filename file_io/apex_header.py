"""
Bloc d'en-tete chiffre Apex v1 du G-code (R-010, contrat C5).

Reimplementation de `apex-ui/src/utils/apexCipher.ts` et `apexHeader.ts`; format et constantes :
`apex-ui/docs/APEX_HEADER_FORMAT.md`, verifie contre `tests/vectors/apex_cipher_vectors.json` (copie par
`scripts/sync_from_ui.py`).

CE N'EST PAS UNE PROTECTION : le projet est open source, toutes les cles sont publiques. Ne jamais le
presenter comme de la securite.

Le moteur ecrit `sliced` et `meta` ; le `journal` est celui de l'interface, que le moteur ne connait pas
(`read_apex_header` le rend tel quel). AES-256-GCM vient du paquet `cryptography` ; sans lui,
`cipher_available()` est faux et l'export n'ecrit que les lignes en clair (jamais une erreur). Les nonces
sont DERIVES du contenu (SHA-256) et non tires au hasard : un meme G-code est reproductible octet pour octet.
Port natif : a faire avec le reste du moteur (R-014) ; Python reste la reference.
"""
from __future__ import annotations

import base64
import hashlib
import json
import zlib
from typing import Optional

try:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
except ImportError:  # pragma: no cover - environnement sans le paquet
    AESGCM = None

MAGIC = b"APX1"
NONCE_BYTES = 12
AES_KEY_1 = bytes.fromhex("4170657820537563726520466f726d617420312d204c6f756973205849562121")
AES_KEY_2 = bytes.fromhex("4165656c6f782048617264776172652d417065782d466f726d61742d32303236")
VIGENERE_KEY = "ROSSIGNOLGRANDCHIFFRE"
A = 65


def cipher_available() -> bool:
    return AESGCM is not None


# --- 1. grand chiffre -----------------------------------------------------------------------------------

def _digrams() -> list:
    french = ("es le de en re nt on er te ou an la et ne se qu ur it ai is me ti ra ce ue ie ar in au ll ns eu co pa pr "
              "ca ch ma po st ta tr ve sa mo lo io ro al el em ri ui ss tt om di ge ni so pe vi li si ea ec ta ux ot un us "
              "ac ag ec ic il ir ol or ul").split(" ")
    js = ['":', '",', '{"', '"}', '[[', ']]', '],', '[', ',"', '":[', '0.', '.0', '-0', 'e-', 'nu', 'll']
    digits = [f"{i:02d}" for i in range(100)]
    merged = list(dict.fromkeys(french + [s for s in js if len(s) == 2] + digits))   # Set de JS : ordre d'insertion
    return [s for s in merged if len(s) == 2]


def _mulberry32(seed: int):
    state = [seed & 0xFFFFFFFF]

    def next_float() -> float:
        state[0] = (state[0] + 0x6D2B79F5) & 0xFFFFFFFF
        t = state[0]
        t = ((t ^ (t >> 15)) * (t | 1)) & 0xFFFFFFFF
        t ^= (t + (((t ^ (t >> 7)) * (t | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0
    return next_float


def _nomenclator():
    entries = [bytes([b]) for b in range(256)] + [d.encode() for d in _digrams()]
    codes = list(range(1000))
    rnd = _mulberry32(1666)
    for i in range(len(codes) - 1, 0, -1):
        j = int(rnd() * (i + 1))
        codes[i], codes[j] = codes[j], codes[i]
    encode = {e: codes[i] for i, e in enumerate(entries)}
    decode = {codes[i]: e for i, e in enumerate(entries)}
    return encode, decode


_ENCODE, _DECODE = _nomenclator()


def _great_encode(data: bytes) -> bytes:
    out, i = [], 0
    while i < len(data):
        pair = _ENCODE.get(data[i:i + 2]) if i + 1 < len(data) else None   # digramme s'il existe, sinon octet seul
        if pair is not None:
            code, i = pair, i + 2
        else:
            code, i = _ENCODE[data[i:i + 1]], i + 1
        out.append(f"{code:03d}")
    return "".join(out).encode("ascii")


def _great_decode(digits: bytes) -> bytes:
    if len(digits) % 3:
        raise ValueError("Format Apex : grand chiffre tronque.")
    out = bytearray()
    for i in range(0, len(digits), 3):
        code = int(digits[i:i + 3])
        if code not in _DECODE:
            raise ValueError("Format Apex : code du grand chiffre inconnu.")
        out += _DECODE[code]
    return bytes(out)


# --- 3. Vigenere ----------------------------------------------------------------------------------------

def _bytes_to_letters(data: bytes) -> str:
    return "".join(chr(A + (b >> 4)) + chr(A + (b & 15)) for b in data)


def _letters_to_bytes(s: str) -> bytes:
    if len(s) % 2:
        raise ValueError("Format Apex : lettres impaires.")
    out = bytearray()
    for i in range(0, len(s), 2):
        hi, lo = ord(s[i]) - A, ord(s[i + 1]) - A
        if not (0 <= hi <= 15 and 0 <= lo <= 15):
            raise ValueError("Format Apex : lettre hors de A-P.")
        out.append((hi << 4) | lo)
    return bytes(out)


def _vigenere(s: str, direction: int) -> str:
    k = VIGENERE_KEY
    return "".join(chr(A + ((ord(c) - A + direction * (ord(k[i % len(k)]) - A)) % 26)) for i, c in enumerate(s))


# --- 4. Fleissner ---------------------------------------------------------------------------------------

GRID = 8
BLOCK = GRID * GRID


def _rotate(p: int) -> int:
    r, c = divmod(p, GRID)
    return c * GRID + (GRID - 1 - r)


def _holes():
    rnd = _mulberry32(1626)
    base = []
    for r in range(GRID // 2):
        for c in range(GRID // 2):
            p = r * GRID + c
            for _ in range(int(rnd() * 4)):
                p = _rotate(p)
            base.append(p)
    out, holes = [], base
    for _ in range(4):
        out.append(sorted(holes))
        holes = [_rotate(h) for h in holes]
    return out


_HOLES = _holes()


def _fleissner_encode(s: str) -> str:
    out = []
    for b in range(0, len(s), BLOCK):
        block, grid, i = s[b:b + BLOCK], [""] * BLOCK, 0
        for holes in _HOLES:
            for h in holes:
                grid[h] = block[i]
                i += 1
        out.append("".join(grid))
    return "".join(out)


def _fleissner_decode(s: str) -> str:
    out = []
    for b in range(0, len(s), BLOCK):
        grid = s[b:b + BLOCK]
        out.extend(grid[h] for holes in _HOLES for h in holes)
    return "".join(out)


# --- 5. Enigma I ----------------------------------------------------------------------------------------

_ROTORS = [("EKMFLGDQVZNTOWYHXUSPAIBRCJ", "Q"), ("AJDKSIRUXBLHWTMCQGZNPYFVOE", "E"), ("BDFHJLCPRTXVZNYEIWGAKMUSQO", "V")]


def _pairs_table(pairs: str) -> list:
    table = list(range(26))
    for pair in pairs.split(" "):
        a, b = ord(pair[0]) - A, ord(pair[1]) - A
        table[a], table[b] = b, a
    return table


_FWD = [[ord(c) - A for c in w] for w, _ in _ROTORS]
_BWD = [[w.index(i) for i in range(26)] for w in _FWD]
_PLUG = _pairs_table("AE BF CM DQ HU JN LX PR SZ VW")
_REFLECTOR = _pairs_table("AY BR CU DH EQ FS GL IP JX KN MO TZ VW")
_NOTCHES = [ord(n) - A for _, n in _ROTORS]


def _enigma(s: str) -> str:
    pos = [ord(c) - A for c in "APX"]

    def through(wiring, c, p):
        return (wiring[(c + p) % 26] - p + 26) % 26
    out = []
    for ch in s:
        if pos[1] == _NOTCHES[1]:
            pos[1] = (pos[1] + 1) % 26
            pos[0] = (pos[0] + 1) % 26
        elif pos[2] == _NOTCHES[2]:
            pos[1] = (pos[1] + 1) % 26
        pos[2] = (pos[2] + 1) % 26
        c = _PLUG[ord(ch) - A]
        for r in (2, 1, 0):
            c = through(_FWD[r], c, pos[r])
        c = _REFLECTOR[c]
        for r in (0, 1, 2):
            c = through(_BWD[r], c, pos[r])
        out.append(chr(A + _PLUG[c]))
    return "".join(out)


# --- assemblage -----------------------------------------------------------------------------------------

def _length_letters(n: int) -> str:
    s = ""
    for _ in range(8):
        s = chr(A + n % 26) + s
        n //= 26
    return s


def _aes_encrypt(key: bytes, plain: bytes, nonce: bytes) -> bytes:
    return nonce + AESGCM(key).encrypt(nonce, plain, None)


def _aes_decrypt(key: bytes, data: bytes) -> bytes:
    return AESGCM(key).decrypt(data[:NONCE_BYTES], data[NONCE_BYTES:], None)


def _derived_nonce(tag: bytes, content: bytes) -> bytes:
    return hashlib.sha256(tag + content).digest()[:NONCE_BYTES]


def apex_encrypt(plain: str, nonce1: Optional[bytes] = None, nonce2: Optional[bytes] = None, trace: Optional[dict] = None) -> str:
    """Chiffre un texte au format Apex v1 (base64). Sans nonces, ils sont derives du contenu."""
    if AESGCM is None:
        raise RuntimeError("le paquet `cryptography` est requis pour chiffrer l'en-tete Apex")
    digits = _great_encode(plain.encode("utf-8"))
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    deflated = comp.compress(digits) + comp.flush()
    stage2 = _aes_encrypt(AES_KEY_1, deflated, nonce1 or _derived_nonce(b"1", deflated))
    stage3 = _vigenere(_bytes_to_letters(stage2), 1)
    padded = stage3 + "X" * ((BLOCK - len(stage3) % BLOCK) % BLOCK)
    stage4 = _length_letters(len(stage3)) + _fleissner_encode(padded)
    stage5 = _enigma(stage4)
    stage6 = _aes_encrypt(AES_KEY_2, stage5.encode("ascii"), nonce2 or _derived_nonce(b"2", stage5.encode("ascii")))
    out = base64.b64encode(MAGIC + stage6).decode("ascii")
    if trace is not None:
        trace.update(greatCipher=digits.decode("ascii"), afterAes1Hex=stage2.hex(), vigenere=stage3,
                     fleissner=stage4, enigma=stage5, output=out)
    return out


def apex_decrypt(encoded: str) -> str:
    if AESGCM is None:
        raise RuntimeError("le paquet `cryptography` est requis pour dechiffrer l'en-tete Apex")
    data = base64.b64decode(encoded.strip())
    if data[:4] != MAGIC:
        raise ValueError("Format Apex : signature absente.")
    stage5 = _aes_decrypt(AES_KEY_2, data[4:]).decode("ascii")
    stage4 = _enigma(stage5)
    length = 0
    for ch in stage4[:8]:
        length = length * 26 + (ord(ch) - A)
    stage3 = _fleissner_decode(stage4[8:])[:length]
    stage2 = _letters_to_bytes(_vigenere(stage3, -1))
    deflated = _aes_decrypt(AES_KEY_1, stage2)
    return _great_decode(zlib.decompress(deflated, -15)).decode("utf-8")


# --- bloc du G-code -------------------------------------------------------------------------------------

BEGIN = "; apex_block_begin v1"
DATA = "; apex_data:"
END = "; apex_block_end"
CHUNK = 100


def header_block_lines(meta: dict) -> list:
    """Lignes du bloc d'en-tete chiffre (`sliced` et `meta`, sans journal), ou `[]` si le chiffrement est
    indisponible. Les cles de `meta` sont libres ; l'ordre est fige pour que le fichier soit reproductible."""
    if not cipher_available():
        return []
    payload = json.dumps({"format": 1, "sliced": True, "meta": meta}, separators=(",", ":"), ensure_ascii=False)
    encoded = apex_encrypt(payload)
    return [BEGIN] + [DATA + encoded[i:i + CHUNK] for i in range(0, len(encoded), CHUNK)] + [END]


def read_apex_header(lines) -> Optional[dict]:
    """Contenu dechiffre du bloc (dict, `journal` de l'interface compris), ou `None` s'il n'y en a pas."""
    lines = [ln.rstrip() for ln in lines]
    try:
        start = lines.index(BEGIN)
        end = lines.index(END, start + 1)
    except ValueError:
        return None
    data = "".join(ln[len(DATA):] for ln in lines[start + 1:end] if ln.startswith(DATA))
    return json.loads(apex_decrypt(data))
