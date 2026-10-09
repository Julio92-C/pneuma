import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
DATA = Path(__file__).resolve().parents[1] / "data"


def load_bin(name):
    """Import bin/<name>.py as a module (bin/ is not a package)."""
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
