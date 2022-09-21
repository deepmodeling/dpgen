#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
from  dpgen import NAME,SHORT_CMD
import setuptools, datetime

readme_file = path.join(path.dirname(path.abspath(__file__)), 'README.md')
try:
    from m2r import parse_from_file
    readme = parse_from_file(readme_file)
except ImportError:
    with open(readme_file) as f:
        readme = f.read()

today = datetime.date.today().strftime("%b-%d-%Y")
with open(path.join('dpgen', '_date.py'), 'w') as fp :
    fp.write('date = \'%s\'' % today)

install_requires=[
    'numpy>=1.14.3',
    'dpdata>=0.2.6',
    'pymatgen>=2019.1.13',
    'ase',
    'monty>2.0.0',
    'paramiko',
    'custodian',
    'GromacsWrapper>=0.8.0',
    'dpdispatcher>=0.3.11',
    'netCDF4',
    'dargs>=0.2.9',
    'pybind11>=2.4',
    'pymatgen-analysis-defects'
]

setuptools.setup(
    name=NAME,
    use_scm_version={'write_to': 'dpgen/_version.py'},
    setup_requires=['setuptools_scm'],
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
              'dpgen/database',
              'dpgen/tools',
              'dpgen/simplify',
              'dpgen/collect',
    ],
    # data_files = [('dpgen/tools/', ['dpgen/tools/update_time.sh', ])],
    # package_data={'example':['*.json']},
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
    keywords='deep potential generator active learning deepmd-kit',
    install_requires=install_requires,    
        entry_points={
          'console_scripts': [
              SHORT_CMD+'= dpgen.main:main']
   }
)

