from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description = "Water properties base on IAPWS Industrial 1997"
long_description = """
droplet!
=======
Water properties base on IAPWS Industrial 1997
http://www.iapws.org/relguide/IF97-Rev.html

Oscar J. Delgado (oscarjdm19@gmail.com)
Feb-May 2016  - London (ON) Canada

The formulation is valid from 273.15 K to 1073.15 K at @100 MPa,
and there is a high-temperature region extending to 2273.15 K @50 MPa.
There is also a separate equation for metastable steam at
pressures up to 10 MPa.

License
=======
"DROPLET" is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2016, Oscar J. Delgado
"""

NAME = "droplet"
MAINTAINER = "Oscar J. Delgado"
MAINTAINER_EMAIL = "lijantropique@protonmail.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Oscar J. Delgado "
AUTHOR_EMAIL = "lijantropique@protonmail.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['droplet',
            'droplet.tests']
PACKAGE_DATA = {'droplet': [pjoin('data', '*')]}
REQUIRES = ["numpy"]
