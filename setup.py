# -*- coding: utf-8 -*-
from setuptools import setup

from build import build

package_dir = {"": "src"}

packages = ["mmseqs"]

package_data = {"": ["*"]}

install_requires = [
    "biopython>=1.78,<2.0",
    "pandas>=1.2.1,<2.0.0",
    "sqlitedict>=1.7.0,<2.0.0",
]

setup_kwargs = {
    "name": "mmseqs",
    "version": "0.1.3",
    "description": "Python bindings for UNAFold to determine hybridization energies / folding of RNA/DNA sequences.",
    "long_description": None,
    "author": "Piotr Styczyński",
    "author_email": "piotrs@radcode.co",
    "maintainer": None,
    "maintainer_email": None,
    "url": None,
    "package_dir": package_dir,
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "python_requires": ">=3.8,<4.0",
}

build(setup_kwargs)

setup(**setup_kwargs)
