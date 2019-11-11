#!/usr/bin/env python3
# Time-stamp: <2019-11-11 12:48:19 taoliu>

"""Description: 

Setup script for genericCTA -- generic Count Table Analysis (dev)

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""
import sys
from setuptools import setup, Extension

def main():
    if float(sys.version[:3])<3.5:
        sys.stderr.write("CRITICAL: Python version must >= 3.5!\n")
        sys.exit(1)

    with open("README.md", "r") as fh:
        long_description = fh.read()
        
    setup(name="genericCTA",
          version="0.0.1",
          description="generic Count Table Analysis (dev)",
          long_description = long_description,
          long_description_content_type="text/markdown",
          author='Tao Liu',
          author_email='vladimir.liu@gmail.com',
          url='http://github.com/taoliu/genericCTA/',
          package_dir={'CTA' : 'CTA'},
          packages=['CTA'],
          scripts=['bin/decomp_cluster_plot', 'bin/build_sp_from_countable', 'bin/print_features_in_topics'],
          classifiers=[
              'Development Status :: 3 - Alpha',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',              
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python :: 3.5',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',              
              ],
          #install_requires=['numpy>=1.17'],          
          #ext_modules = ext_modules
          )

if __name__ == '__main__':
    main()
