#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
from  dpgen import NAME,SHORT_CMD
import setuptools

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file
    readme = parse_from_file(readme_file)
except ImportError:
    with open(readme_file) as f:
        readme = f.read()

install_requires=['numpy>=1.14.3', 'dpdata>=0.1.5', 'pymatgen>=2017.9.1', 'ase', 'monty>2.0.0', 'paramiko', 'custodian']

setuptools.setup(
    name=NAME,
    version_format='{tag}.dev{commitcount}+{gitsha}',
    setup_requires=['setuptools-git-version'],
    author="Han Wang",
    author_email="wang_han@iapcm.ac.cn",
    description="DPGen: The deep potential generator",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/deepmodeling/dpgen",
    python_requires="~=3.6",
    packages=['dpgen', 
              'dpgen/generator',
              'dpgen/generator/lib',
              'dpgen/auto_test',
              'dpgen/auto_test/lib',
              'dpgen/data',
              'dpgen/data/tools',
              'dpgen/remote',
              'dpgen/dispatcher',
              'dpgen/database'
    ],
    # package_data={'example':['*.json']},
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    keywords='lammps vasp deepmd-kit',
    install_requires=install_requires,    
        entry_points={
          'console_scripts': [
              SHORT_CMD+'= dpgen.main:main']
   }
)

