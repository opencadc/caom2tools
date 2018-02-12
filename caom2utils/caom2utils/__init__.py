# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
TODO
"""

from .fits2caom2 import *  # noqa
from .legacy import *  # noqa

import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
