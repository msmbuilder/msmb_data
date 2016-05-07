from __future__ import print_function, absolute_import, division

from setuptools import setup, find_packages

VERSION = '1.0.0.dev0'
ISRELEASED = False
__version__ = VERSION


classifiers = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: Apache Software License
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
setup(
    name='msmb_data',
    author='Carlos Xavier Hernandez',
    author_email='cxh@stanford.edu',
    url='https://github.com/msmbuilder/msmb_data',
    classifiers=[e.strip() for e in classifiers.splitlines()],
    platforms=["Windows", "Linux", "Mac OS-X", "Unix"],
    license='Apache Software License',
    version=VERSION,
    packages=find_packages(),
    zip_safe=False,
    package_data={'msmb_data': ['data/*/*.*']},
)
