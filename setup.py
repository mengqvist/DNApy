#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
#
#dnapy is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#
#colcol is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Library General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software Foundation,
#Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import io
import re
import sys

from setuptools import setup, find_packages

if sys.version_info[0] == 2:
    sys.exit("Sorry, only Python 3 is supported by this package.")

# requirements
install_requirements = []

# package informations
with io.open('dnapy/__init__.py', 'r') as f:
    text = f.read()
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]$', text,
                        re.MULTILINE).group(1)
    desc = re.search(r'^__desc__\s*=\s*[\'"]([^\'"]*)[\'"]$', text,
                        re.MULTILINE).group(1)

if not version:
    raise RuntimeError('Cannot find version information')

with io.open('README.md', 'r', encoding='utf-8') as f:
    readme = f.read()

#with io.open('HISTORY.md', 'r', encoding='utf-8') as f:
#    history = f.read()


setup(name='dnapy',
        version=version,
        license='GPLv3+',
        description=desc,
    	author='Martin Engqvist',
    	author_email='martin.engqvist@chalmers.se',
        #long_description=readme + '\n\n' + history,
        keywords='DNA plasmid editing',
        packages=find_packages(exclude=['contrib', 'docs', 'tests*']), #find folders containing scripts, exclude irrelevant ones
        #entry_points={'console_scripts': ['DNApy = DNApy.gui.main_GUI:startgui'],},
        install_requires=install_requirements,
    	classifiers=[
    	# How mature is this project? Common values are
    	#   3 - Alpha
    	#   4 - Beta
    	#   5 - Production/Stable
    	'Development Status :: 3 - Alpha',

    	# Indicate who your project is intended for
    	'Intended Audience :: Science/Research',

    	# Pick your license as you wish (should match "license" above)
    	'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

    	# Specify the Python versions you support here. In particular, ensure
    	# that you indicate whether you support Python 2, Python 3 or both.
    	'Programming Language :: Python :: 3 :: Only',
    	'Programming Language :: Python :: 3',
    	'Programming Language :: Python :: 3.2',
    	'Programming Language :: Python :: 3.3',
    	'Programming Language :: Python :: 3.4',
    	'Programming Language :: Python :: 3.5',
    	'Programming Language :: Python :: 3.6'],
        python_requires='>=3' #python version
)
