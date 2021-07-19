import importlib
from pathlib import Path
try:
    from importlib import metadata
except ImportError: # for Python<3.8
    try:
        import importlib_metadata as metadata
    except ImportError:
        metadata = None

PACKAGE_BASE = "dpgen.generator.styles"
NOT_LOADABLE = ("__init__.py",)

for module_file in Path(__file__).parent.glob("*.py"):
    if module_file.name not in NOT_LOADABLE:
        module_name = f".{module_file.stem}"
        importlib.import_module(module_name, PACKAGE_BASE)

# https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
if metadata:
    eps = metadata.entry_points().get('dpdata.styles', [])
    for ep in eps:
        plugin = ep.load()
