# File: setup.py
# Date: 10-Oct-2018
#
# Update:
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
<<<<<<< HEAD
thisPackage = 'wwpdb.utils.emdb'

with open('wwpdb/utils/emdb/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError('Cannot find version information')
=======
thisPackage = "wwpdb.utils.emdb"

with open("wwpdb/utils/emdb/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")
>>>>>>> origin/develop

setup(
    name=thisPackage,
    version=version,
<<<<<<< HEAD
    description='wwPDB EMDB to XML converter',
    long_description="See:  README.md",
    author='Ezra Peisach',
    author_email='ezra.peisach@rcsb.org',
    url='https://github.com/rcsb/py-wwpdb_utils_config',
    #
    license='Apache 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        # 'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points={
        'console_scripts': ['em_emd_conversion=wwpdb.utils.config.emdb.cif_emdb_translator:main']
    },
    #
    install_requires=['lxml', 'mmcif', 'wwpdb.utils.config'],
    packages=find_packages(exclude=['wwpdb.utils.tests-emdb', 'mock-data', 'tests.*', 'wwpdb.utils.emdb/cif_emdb)translator/*test*.py']),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        '': ['*.md', '*.rst', "*.txt", "*.cfg"],
=======
    description="wwPDB EMDB to XML converter",
    long_description="See:  README.md",
    author="Ezra Peisach",
    author_email="ezra.peisach@rcsb.org",
    url="https://github.com/rcsb/py-wwpdb_utils_config",
    #
    license="Apache 2.0",
    classifiers=[
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    entry_points={"console_scripts": ["em_emd_conversion=wwpdb.utils.config.emdb.cif_emdb_translator:main"]},
    #
    install_requires=["lxml", "mmcif", "wwpdb.utils.config"],
    packages=find_packages(exclude=["wwpdb.utils.tests-emdb", "mock-data", "tests.*", "wwpdb.utils.emdb/cif_emdb)translator/*test*.py"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"],
>>>>>>> origin/develop
    },
    #
    # These basic tests require no database services -
    test_suite="wwpdb.utils.tests-emdb",
<<<<<<< HEAD
    tests_require=['tox', 'wwpdb.utils.testing'],
    #
    # Not configured ...
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },
    # Added for
    command_options={
        'build_sphinx': {
            'project': ('setup.py', thisPackage),
            'version': ('setup.py', version),
            'release': ('setup.py', version)
        }
    },
=======
    tests_require=["tox", "wwpdb.utils.testing"],
    #
    # Not configured ...
    extras_require={
        "dev": ["check-manifest"],
        "test": ["coverage"],
    },
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
>>>>>>> origin/develop
    # This setting for namespace package support -
    zip_safe=False,
)
