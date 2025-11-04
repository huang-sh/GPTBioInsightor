# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys
import re
from pathlib import Path
from importlib.metadata import PackageNotFoundError, metadata

sys.path.insert(0, os.path.abspath("../../.."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = "gptbioinsightor"
copyright = "2024, huangsh"
author = "huangsh"

try:
    info = metadata("gptbioinsightor")
    release = info["Version"]
except PackageNotFoundError:
    init_path = Path(__file__).resolve().parents[3] / "gptbioinsightor" / "__init__.py"
    version_match = re.search(r'__version__\s*=\s*"([^"]+)"', init_path.read_text())
    release = version_match.group(1) if version_match else "0.0.0"

version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns = []

language = "en"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"


html_theme_options = {
    "repository_url": "https://github.com/huang-sh/GPTBioInsightor",
    "use_repository_button": True,
}

html_title = "GPTBioInsightor"

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
myst_heading_anchors = 6  # create anchors for h1-h6

nb_execution_mode = "off"

autodoc_mock_imports = [
    "anndata",
    "scanpy",
    "gseapy",
    "openai",
    "anthropic",
    "litellm",
    "perplexityai",
    "instructor",
    "frozendict",
]
