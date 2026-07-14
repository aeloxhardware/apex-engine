"""Worker asynchrone exécutant le script de traitement 'exec.py' dans un thread séparé."""
import sys
import os
import json
import subprocess
from PyQt6.QtCore import QThread, pyqtSignal


class ProcessWorker(QThread):
    """Gère l'exécution asynchrone du script de traitement 'exec.py'."""
    progress_signal = pyqtSignal(dict)
    waiting_signal = pyqtSignal()
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(int)

    def __init__(self, script_path):
        super().__init__()
        self.script_path = script_path
        self.process = None

    def run(self):
        try:
            resolved_path = self.script_path
            if not os.path.isabs(resolved_path):
                resolved_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), resolved_path)

            if not os.path.isfile(resolved_path):
                self.error_signal.emit(f"SCRIPT INTROUVABLE : {resolved_path}")
                self.finished_signal.emit(-1)
                return

            env = os.environ.copy()
            env["PYTHONIOENCODING"] = "utf-8"
            
            self.process = subprocess.Popen(
                [sys.executable, resolved_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE,
                text=True, encoding="utf-8", bufsize=1, env=env 
            )
            
            for line in iter(self.process.stdout.readline, ''):
                stripped = line.strip()
                if "STATUS:WAITING_FOR_USER_INPUT" in stripped:
                    self.waiting_signal.emit()
                    continue
                try:
                    data = json.loads(stripped)
                    if isinstance(data, dict) and data.get("type") == "progress":
                        self.progress_signal.emit(data)
                except json.JSONDecodeError: pass
                    
            _, stderr = self.process.communicate()
            exit_code = self.process.returncode

            if stderr:
                filtered = [l for l in stderr.splitlines() if "UserWarning" not in l and "warnings.warn" not in l]
                if filtered: self.error_signal.emit("\n".join(filtered[-6:]))

            self.finished_signal.emit(exit_code if exit_code is not None else -1)
        except Exception as e:
            self.error_signal.emit(f"ERREUR PROCESS: {str(e)}")
            self.finished_signal.emit(-1)
