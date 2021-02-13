# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

#
# -- Path preparation/config loading ----------------------------------------------------
#

import sys, os
import shutil
import toml
from pathlib import Path
def copytree(src, dst, symlinks=False, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                shutil.copy2(s, d)

pyproject = toml.load(Path(__file__).parents[1].joinpath('pyproject.toml'))
sys.path.append(str(Path(__file__).parents[1]))

sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('_ext'))
sys.path.insert(0, os.path.abspath('.'))
#sys.path.insert(0, os.path.abspath('../src/biosnake_native'))

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
    'add_links',
    'sphinx_pyreverse',
    'sphinx_sitemap',
    #'sphinxcontrib.apidoc',
    'sphinx_gitstamp',
    'sphinx_last_updated_by_git',
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
    'sphinxcontrib.needs',
    'sphinxcontrib.test_reports',
    'sphinxcontrib.plantuml'
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

html_baseurl = 'https://my-site.com/docs/'

gitstamp_fmt = "%b %d, %Y"

#
# -- Options for the theme ----------------------------------------------------
#

html_logo = "logo.png"
html_theme = "furo"
html_theme_options = {
    "style_external_links": True,
    "navigation_depth": 2
}

#
# -- Options for HTML output --------------------------------------------------
#

html_title = full_title
htmlhelp_basename = "Biosnakehandbookdoc"
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
    from utils.import_coverage import import_coverage_report
    import_coverage_report()
    from shutil import copyfile
    from pathlib import Path
    Path("./_sphinx_resources").mkdir(parents=True, exist_ok=True)
    copyfile('../README.md', './_sphinx_resources/README.md')
    copytree('js', './_static/js')
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


# --------

# conf.py
# srclink settings
srclink_project = 'https://github.com/styczynski/python-biosnake'
#srclink_project = 'https://bitbucket.org/westurner/sphinxcontrib-srclinks'
#srclink_project = 'hg@bitbucket.org/westurner/sphinxcontrib-srclinks'
#srclink_project = 'git@bitbucket.org/westurner/sphinxcontrib-srclinks'
srclink_src_path = '.doc/'
#srclink_src_path = ''
srclink_branch = 'master'
#srclink_branch = 'develop'

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': [
        "sidebar/scroll-start.html",
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/navigation.html",
        "sidebar/ethical-ads.html",
        "sidebar/scroll-end.html",
        'srclinks.html',
    ]
}

html_js_files = [
    'js/cleaner_navigation.js',
]

# plant uml
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    plantuml = 'java -Djava.awt.headless=true -jar /usr/share/plantuml/plantuml.jar'
else:
    cwd = os.getcwd()
    plantuml = 'java -jar %s' % os.path.join(cwd, "utils/plantuml_beta.jar")

# If we are running on windows, we need to manipulate the path,
# otherwise plantuml will have problems.
if os.name == "nt":
    plantuml = plantuml.replace("/", "\\")
    plantuml = plantuml.replace("\\", "\\\\")

plantuml_output_format = 'png'