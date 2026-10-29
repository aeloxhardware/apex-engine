"""Garde-fou memoire des calculs lourds du lissage (porte de `pipeline/bird/skin/skin_carve.py`, branche `bird_old` :
un calcul non borne avait sature 47 Go). Windows seulement (sans psutil) ; ailleurs, sans effet."""
import os

MAX_PROCESS_BYTES = 12 * 2**30  # memoire du processus au-dela de laquelle on s'arrete


def process_bytes():
    """Memoire de travail du processus (Windows, sans psutil ; 0 ailleurs)."""
    if os.name != "nt":
        return 0
    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
    c = Counters()
    c.cb = ctypes.sizeof(c)
    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
    k32.K32GetProcessMemoryInfo.restype = wintypes.BOOL
    if not k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb):
        raise OSError("GetProcessMemoryInfo a echoue : garde-fou memoire inoperant")
    return c.WorkingSetSize


def check_memory(where):
    used = process_bytes()
    if used > MAX_PROCESS_BYTES:
        raise MemoryError(f"{where} : {used / 2**30:.1f} Go utilises (> {MAX_PROCESS_BYTES / 2**30:.1f}) -- arret")
