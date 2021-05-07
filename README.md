# Biosnake

![PyPi version](https://img.shields.io/pypi/v/biosnake?style=flat-square)
![Build status](https://img.shields.io/github/workflow/status/covid-genomics/biosnake/Build?style=flat-square)
[![Maintainability](https://api.codeclimate.com/v1/badges/bb922f4c1b788ed1a6a4/maintainability)](https://codeclimate.com/github/covid-genomics/biosnake/maintainability)
![Supported Python versions](https://img.shields.io/pypi/pyversions/biosnake?style=flat-square)
![Python implementation](https://img.shields.io/pypi/implementation/biosnake?style=flat-square)
![Commit activity](https://img.shields.io/github/commit-activity/w/covid-genomics/biosnake?style=flat-square)

## Project description

This project aims at providing unified biological Python processing framework.
The main idea is that the bio-datascience sector lacks an efficient, modern framework that can effectively and easily process large volumes of data.
Biopython is slow, but Biotite solves that problem. Nevertheless both tools do not provide tools to common problems that can be solved only using external legacy programs.

## Project roadmap

- Create a MMSeqs2 PoC Python bindings
- Provide bindings for all Biosnake2 functionalities with pandas/numpy compatibility
- Provide visualizations that use custom components or Plotly (similar to those of Biotite)
- Incorporate [Phylip](https://evolution.genetics.washington.edu/phylip/progs.data.prot.html)
- Incorporate protein [algorithms](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6394400/) - or similar to [this](https://pypi.org/project/seqfold/)

## Current project status

Currently the project in state of heavy development. The API is not stable and can change daily. It's rather RnD state of project than actually production verified quality.
* Currently Biosnake2 bindings are implemented

## Building

You can build the project by calling:
```bash
    $ poetry install
```
You need to install [Poetry](https://python-poetry.org/docs/#installation) to execute that.

## Usage

Please see `example/` folder to see usages of the library.
