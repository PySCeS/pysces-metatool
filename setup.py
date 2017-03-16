#!/usr/bin/env python

"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2017 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""

__doc__ = "PySCeS/MetaTool add-on module"
__version__ = '0.7.1'

import os

print """\nMETATOOL is a C program developed from 1998 to 2000 by Thomas Pfeiffer
(Berlin) in cooperation with Stefan Schuster and Ferdinand
Moldenhauer (Berlin) and Juan Carlos Nuno (Madrid).\n"""

try:
    from numpy.distutils.core import setup, Extension
except Exception, ex:
    print ex
    print "This requires NumPy and SciPy 0.5x+\n"
    os.sys.exit(-1)

########## From here on it's up to distutils ##########

# get the dir of setup.py
local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))
os.chdir(local_path)

myscripts = []
mydata_files = []

print '\nBuilding metatool'
# Compile the metatool binaries and setup the data file list
meta_path = os.path.join(local_path,'pysces','metatool')
os.chdir(meta_path)
if os.sys.platform == 'win32':
    if not os.path.exists(os.path.join(meta_path,'meta43_int.exe')) or not os.path.exists(os.path.join(meta_path,'meta43_double.exe')):
        os.spawnl(os.P_WAIT,'build_win32.bat','build_win32.bat')
    mydata_files = [(os.path.join('pysces','metatool'), [os.path.join(local_path, 'pysces', 'metatool','meta43_double.exe'), os.path.join(local_path, 'pysces', 'metatool','meta43_int.exe')])]
    mydata_files.append((os.path.join('pysces','metatool'), [os.path.join(local_path, 'pysces', 'metatool','readme.txt'), os.path.join(local_path, 'pysces', 'metatool','readme.txt')]))
else:
    print '\nBuilding metatool...\n(mktemp warnings do not affect the compile process)'
    #change mode to rwxrwxr--
    os.chmod('build_linux',508)
    if not os.path.exists(os.path.join(meta_path,'meta43_int')) or not os.path.exists(os.path.join(meta_path,'meta43_double')):
        os.spawnl(os.P_WAIT,os.path.join(local_path,'pysces','metatool','build_linux'), os.path.join(local_path,'pysces','metatool','build_linux'))
    mydata_files = [(os.path.join('pysces','metatool'), [os.path.join(local_path, 'pysces', 'metatool','meta43_double'), os.path.join(local_path, 'pysces', 'metatool','meta43_int')])]
    mydata_files.append((os.path.join('pysces','metatool'), [os.path.join(local_path, 'pysces', 'metatool','readme.txt'), os.path.join(local_path, 'pysces', 'metatool','readme.txt')]))
os.chdir(local_path)

# my subpackage list
mypackages= ['pysces.metatool']

os.chdir(local_path)
# Install packages and the metatool binaries as "data"
setup(name="pysces_metatool",
    version = __version__,
    description = "PySCeS MetaTool module",
    long_description = """
METATOOL is a C program developed from 1998 to 2000 by Thomas Pfeiffer
(Berlin) in cooperation with Stefan Schuster and Ferdinand
Moldenhauer (Berlin) and Juan Carlos Nuno (Madrid).
It serves to derive conclusions about the pathway
structure of metabolic networks from the stoichiometric reaction
equations and information about reversibility  and irreversibility of
enzymes.
It should preferably be compiled with the GNU compiler.
For DOS and Win32 console applications, comment out the two lines
#include<conio.h> and #include<malloc.h>.

Python compatible compile scripts - Brett G. Olivier

    """,
    author = "Thomas Pfeiffer, Stefan Schuster, Ferdinand Moldenhauer, Juan Carlos Nuno",
    maintainer = "",
    maintainer_email = "",
    url = "http://pysces.sourceforge.net",
    download_url = "http://pysces.sourceforge.net/download.html",
    license = "",
    keywords = "elementary modes" ,
    zip_safe = False,
    requires = ['pysces'],
    platforms = ["Windows", "Linux"],
    classifiers = [
    'Environment :: Console',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry'],
    packages = mypackages,
    data_files = mydata_files,
    )
