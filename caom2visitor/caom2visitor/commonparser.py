
import logging
from argparse import ArgumentParser
import os
import sys
from __version__ import version



class CommonParser(ArgumentParser):
    """Class to hold and parse common command-line arguments for vos clients"""

    def __init__(self, *args, **kwargs):
        # call the parent constructor
        ArgumentParser.__init__(self, *args, **kwargs)

        # inherit the VOS client version
        self.version = version
        self.log_level = logging.ERROR

        # now add on the common parameters
        self.add_argument("--cert", metavar=('<CertFile>'),
                        help="location of your CADC security certificate file",
                        default=os.path.join(os.getenv("HOME", "."),
                                             ".ssl/cadcproxy.pem"))
        self.add_argument("--token", metavar=('<TokenString>'),
                        help="token string (alternative to certfile)",
                        default=None)
        self.add_argument('--version', action='version', version='%(prog)s ' +\
                          str(self.version))
        self.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Print debug level log messages")
        self.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Print verbose level log messages")
        self.add_argument("-w", "--warning", action="store_true", default=False,
                        help="Print warning level log messages")





