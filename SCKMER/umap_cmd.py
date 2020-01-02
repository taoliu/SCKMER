# Time-stamp: <2019-12-31 15:18:48 taoliu>

"""Description: SCKMER umap cmd

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime
import tempfile

# ------------------------------------
# own python modules
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """UMAP
    
    """
    # Parse options...
    options = args
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error

    
