# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package implements a client for accessing a caom2repo Web service (WS)

The package can be used as a library as well as through the caom2-repo-client
command it installs.

1. Instantiate CAOM2RepoClient and use it to access the caom2repo WS.

The only mandatory argument that the CadcDataClient constructor takes is
a cadcutils.net.Subject that holds the user credentials. The client interacts
with the caom2repo WS through the get_observation, put_observation,
post_observation, delete_observation and visit functions of the client.

Example (appropriate credentials required):
   from caom2repo import CAOM2RepoClient
   from cadcutils import net

   client = CAOM2RepoClient(net.Subject(certificate='/path/mycert.pem'))
   print(client.get_observation('CFHT', '700000'))

2. Invoke the caom2repo entry point function. This is the function that
is used to generate the caom2-repo-client application

Example:
   from caom2repo import main_app
   import sys

   sys.argv = ['caom2-repo-client', 'read', '--netrc', '/path/mynetrc',
              --collection', 'CFHT', '700000']
   main_app()

3. Invoke the caom2-repo-client as an external command
Example:
   import os
   os.system('caom2-repo-client read -n --collection CFHT 700000')

Method 1. is the recommended method as it does not required forking external
processes and also allows trapping the exceptions and reacting according to
the type of the error. Method 2 also works but the sys.exit needs to be
trapped in order to prevent the script from quiting. Method 3, while simple,
must rely on inter processes communication to determine the result of
running the command.
"""

from .core import *  # noqa
