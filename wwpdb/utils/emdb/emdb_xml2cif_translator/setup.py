from distutils.core import setup
from setuptools import find_packages

setup(
    name='emdb_xml2cif_translator',
    version='1.0',
    packages=find_packages(),
    url='',
    license='',
    author='sanja',
    author_email='sanja@ebi.ac.ik',
    description='EMDB XML to CIF translator package',
    zip_safe=False,
    install_requires=[
            'gitpython',
            'wwpdb.utils.emdb'
        ]
)
