# -*- coding: utf-8 -*-


__title__ = 'DNApy'
__version__ = '0.1'
__desc__ = 'DNA plasmid editing'
__license__ = 'GNU General Public License v3 or later (GPL3+)'
__url__ = 'https://github.com/mengqvist/DNApy'

import logging
import os

try:
    # Try systemwide installation and fallback to included
    import pyperclip
except ImportError:
    from .external import pyperclip

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

DEFAULT_SETTINGS_DIR = os.path.join(ROOT_DIR, "defaults")
SETTINGS_DIR = DEFAULT_SETTINGS_DIR

RESOURCES_DIR = os.path.join(ROOT_DIR, "resources")

GUI_DIR = os.path.join(ROOT_DIR, "gui")
ICONS_DIR = os.path.join(GUI_DIR, "icon")

logging.getLogger(__name__).addHandler(logging.NullHandler())
