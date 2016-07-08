from __future__ import print_function, absolute_import, division

import numpy as np
from os.path import join as pjoin
from setuptools import setup, find_packages, Extension

VERSION = '1.1.0'
ISRELEASED = True
__version__ = VERSION


classifiers = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: Python
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3
Programming Language :: Python :: 3.4
Programming Language :: Python :: 3.5
Operating System :: Unix
Operating System :: MacOS
Operating System :: Microsoft :: Windows
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Information Analysis"""

extensions = []
extensions.append(
    Extension('msmb_data._muller',
              sources=[pjoin('msmb_data', '_muller.pyx')],
              include_dirs=[np.get_include()]))

setup(
    name='msmb_data',
    author='Carlos Xavier Hernandez',
    author_email='cxh@stanford.edu',
    url='https://github.com/msmbuilder/msmb_data',
    classifiers=[e.strip() for e in classifiers.splitlines()],
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    version=VERSION,
    packages=find_packages(),
    zip_safe=False,
    package_data={'msmb_data': ['data/*/*.*']},
    ext_modules=extensions
)
