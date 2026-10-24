"""Bundle the canonical public skills with binary distributions."""

from __future__ import annotations

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildPy(build_py):
    """Copy public workflow documentation without a second maintained copy."""

    def run(self):
        super().run()
        source = Path(__file__).parent / "skills"
        destination = Path(self.build_lib) / "dpgen" / "_skills"
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)


setup(cmdclass={"build_py": BuildPy})
