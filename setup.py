#!/usr/bin/env python3

from setuptools import setup, find_packages
from os import path

PKG_NAME = "Counting Peptide Inserts"
MOD_NAME = "CoPI"

DESCRIPTION  = f'''
Counting Peptide Inserts (CoPI) quantify Peptide Variants for Protein Engineering. 
Options include to Collapse variants with low counts.
'''

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as fh, open(path.join(here, "requirements.txt")) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""

exec(open("{}/_version.py".format(MOD_NAME)).read())

setup(
    name=PKG_NAME,
    version=__version__,
    author="M Irfan",
    author_email="bioirfanatics@gmail.com",
    description=f"{PKG_NAME}",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/chewlabSB2/CoPI",
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_dir={'CoPI': 'CoPI'}, 
    entry_points={
        "console_scripts": [
            "CoPI={}.analysis:main".format(MOD_NAME),
        ],
    },
    
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)
