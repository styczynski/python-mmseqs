# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

#
# -- Path preparation/config loading ----------------------------------------------------
#

import sys
import toml
from pathlib import Path

pyproject = toml.load(Path(__file__).parents[1].joinpath('pyproject.toml'))
sys.path.append(str(Path(__file__).parents[1]))

import os, sys
sys.path.insert(0, os.path.abspath('../src'))
#sys.path.insert(0, os.path.abspath('../src/mmseqs_native'))

# -- Project information -----------------------------------------------------

project_poetry_conf = pyproject["tool"]["poetry"]
project = project_poetry_conf["name"]
copyright = ', '.join(project_poetry_conf["authors"])
author = ', '.join(project_poetry_conf["authors"])
full_title = project_poetry_conf["name"] + " Documentation"


# The full version, including alpha/beta/rc tags
version = project_poetry_conf["version"]

# -- General configuration ---------------------------------------------------

extensions = [
    #'numpydoc',
    'sphinx_pyreverse',
    #'sphinxcontrib.apidoc',
    'sphinx.ext.autodoc',
    'sphinx_git',
    'sphinx_autodoc_typehints',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon',  # numpy style docstrings
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx_markdown_builder',
    'sphinxcontrib.spelling',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxcontrib.mermaid',
    'sphinx_inline_tabs',
    'sphinx_click',
    'sphinx.ext.graphviz',
    'sphinx.ext.inheritance_diagram',
]

master_doc = 'index'
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

#
# -- Options for extensions ----------------------------------------------------
#

mermaid_version = "8.5.2"
doctest_global_setup = '''
try:
    import pandas as pd
except ImportError:
    pd = None
try:
    import numpy as np
except ImportError:
    np = None
'''

apidoc_toc_file = False
apidoc_module_dir = '../src'
apidoc_output_dir = '_reference'
apidoc_excluded_paths = ['tests']
apidoc_separate_modules = True

autosummary_generate = True

#
# -- Options for the theme ----------------------------------------------------
#

html_theme = "furo"
html_theme_options = {
    "style_external_links": True,
    "navigation_depth": 2
}

#
# -- Options for HTML output --------------------------------------------------
#

html_title = full_title
htmlhelp_basename = "MMseqshandbookdoc"
html_static_path = ['_static']

source_suffix = ['.rst']

#
# -- Options for spelling -----------------------------------------------------
#

spelling_lang = "en_GB"
spelling_word_list_filename = "spelling_wordlist"
spelling_ignore_pypi_package_names = True


# def setup(app):
#     print('Setup Sphinx... LOLAYY')
#     from shutil import copyfile
#     from pathlib import Path
#     Path("./_reference").mkdir(parents=True, exist_ok=True)
#     copyfile('module_index.rst', '_reference/modules.rst')

#
# -- Documentation options -----------------------------------------------------
#

add_module_names = False

# for CommonMarkParser
from recommonmark.parser import CommonMarkParser

def setup(app):
    from shutil import copyfile
    from pathlib import Path
    Path("./_sphinx_resources").mkdir(parents=True, exist_ok=True)
    copyfile('../README.md', './_sphinx_resources/README.md')
    app.add_source_suffix('.md', 'markdown')
    app.add_source_parser(CommonMarkParser)
    app.add_config_value('markdown_parser_config', {
        'auto_toc_tree_section': 'Content',
        'enable_auto_doc_ref': True,
        'enable_auto_toc_tree': True,
        'enable_eval_rst': True,
        'enable_inline_math': True,
        'enable_math': True,
    }, True)