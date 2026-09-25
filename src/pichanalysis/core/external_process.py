"""Launch external tools without leaking PyInstaller's DLL search directory."""

from __future__ import annotations

import subprocess
import sys
import threading


_DLL_LOCK = threading.RLock()


def run_external(command, *, timeout: int | float, **kwargs) -> subprocess.CompletedProcess:
    """Equivalent to captured subprocess.run, isolating a frozen Windows child.

    PyInstaller prepends its private Qt DLL directory through SetDllDirectoryW.
    Child Rscript must inherit the normal Windows loader path instead. Restore
    the parent's directory immediately after CreateProcess, not after R exits.
    """
    if not (sys.platform == "win32" and getattr(sys, "frozen", False)):
        return subprocess.run(command, timeout=timeout, check=False, **kwargs)

    import ctypes

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetDllDirectoryW.argtypes = (ctypes.c_uint32, ctypes.c_wchar_p)
    kernel.GetDllDirectoryW.restype = ctypes.c_uint32
    kernel.SetDllDirectoryW.argtypes = (ctypes.c_wchar_p,)
    kernel.SetDllDirectoryW.restype = ctypes.c_int
    buffer = ctypes.create_unicode_buffer(32768)
    with _DLL_LOCK:
        length = kernel.GetDllDirectoryW(len(buffer), buffer)
        if length >= len(buffer):
            raise OSError("Current Windows DLL directory is too long to preserve safely.")
        original = buffer.value if length else None
        if not kernel.SetDllDirectoryW(None):
            raise OSError(ctypes.get_last_error(), "Could not isolate external process DLL search path.")
        try:
            capture = kwargs.pop("capture_output", False)
            process = subprocess.Popen(command, **{
                "stdout": subprocess.PIPE if capture else None,
                "stderr": subprocess.PIPE if capture else None,
                **kwargs,
            })
        finally:
            if not kernel.SetDllDirectoryW(original):
                raise OSError(ctypes.get_last_error(), "Could not restore application DLL search path.")
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        process.communicate()
        raise
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
