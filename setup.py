"""
openmmdl
Script for preparation and simulation of protein-ligand complexes with OpenMM.
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

########################
__version__ = '1.0'
VERSION = __version__
########################

requirements = [
'numpy>=1.23.3',
'requests>=2.28.1',
'rdkit>=2022.03.5',
'mdtraj=1.9.7',
'pdbfixer=1.8.1',
'openff-toolkit=0.11.1',
'openmmforcefields=0.11.2',
'cudatoolkit>=11.7.0',
'cuda-python>=11.7.1',
'MDAnalysis>=2.3.0',
'pytest-shutil>=1.7.0',
'flask>=2.2.2'
]

setup(
    name='OpenmMMDL',
    version=__version__,
    description="Script for preparation and simulation of protein-ligand complexes with OpenMM",
    license="MIT",
    author="Valerij Talagayev & Yu Chen",
    author_email='v.talagayev@fu-berlin.de',
    url='https://github.com/pacificore/OpenmMMDL',
    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),
    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,
    # Entry point
    entry_points={"console_scripts": ['mmdl-setup = openmmdl_setup.openmmdlsetup:main', 'mmdl-simulation = openmmdl_setup.openmmdlsetup:main']},
    install_requires=requirements,
    keywords='openmmmdl',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
