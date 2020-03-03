#!/usr/bin/python3

from distutils.core import setup, Extension

module1 = Extension('he3lib',
                    libraries = ['he3'],
                    sources = ['he3module.c'])

setup (name = 'he3lib',
       version      = '1.0',
       description  = 'He3 library',
       author       = 'Vladislav Zavjalov',
       author_email = 'slazav at altlinux.org',
       url          = 'http://slazav.github.io/he3lib/index.html',
       ext_modules = [module1])
