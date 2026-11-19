import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
DATA = Path(__file__).resolve().parents[1] / "data"

sys.path.insert(0, str(BIN))   # bin scripts import each other (e.g. circular_candidates -> extend_context)


def load_bin(name):
    """Import bin/<name>.py as a module (bin/ is not a package)."""
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


TEST_TAXA = [("1", "1", "no rank", "root"), ("2", "1", "superkingdom", "Bacteria"),
             ("1224", "2", "phylum", "Pseudomonadota"), ("543", "1224", "family", "Enterobacteriaceae"),
             ("561", "543", "genus", "Escherichia"), ("562", "561", "species", "Escherichia coli"),
             ("1239", "2", "phylum", "Bacillota"), ("81852", "1239", "family", "Enterococcaceae"),
             ("1350", "81852", "genus", "Enterococcus"), ("1351", "1350", "species", "Enterococcus faecalis"),
             ("33958", "1239", "family", "Lactobacillaceae"), ("1578", "33958", "genus", "Lactobacillus")]


def write_taxdump(directory, taxa):
    """Minimal NCBI taxdump: taxa = [(taxid, parent, rank, name), ...] -> directory with nodes.dmp and names.dmp."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "nodes.dmp").write_text("".join("%s\t|\t%s\t|\t%s\t|\t\t|\n" % (t, p, r) for t, p, r, _n in taxa))
    (directory / "names.dmp").write_text(
        "".join("%s\t|\t%s\t|\t\t|\tscientific name\t|\n" % (t, n) for t, _p, _r, n in taxa))
    return str(directory)
