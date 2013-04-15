#/*+
#************************************************************************
#****  C A N A D I A N   A S T R O N O M Y   D A T A   C E N T R E  *****
#*
#* (c) 2012.                            (c)2012.
#* National Research Council            Conseil national de recherches
#* Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#* All rights reserved                  Tous droits reserves
#*
#* NRC disclaims any warranties,        Le CNRC denie toute garantie
#* expressed, implied, or statu-        enoncee, implicite ou legale,
#* tory, of any kind with respect       de quelque nature que se soit,
#* to the software, including           concernant le logiciel, y com-
#* without limitation any war-          pris sans restriction toute
#* ranty of merchantability or          garantie de valeur marchande
#* fitness for a particular pur-        ou de pertinence pour un usage
#* pose.  NRC shall not be liable       particulier.  Le CNRC ne
#* in any event for any damages,        pourra en aucun cas etre tenu
#* whether direct or indirect,          responsable de tout dommage,
#* special or general, consequen-       direct ou indirect, particul-
#* tial or incidental, arising          ier ou general, accessoire ou
#* from the use of the software.        fortuit, resultant de l'utili-
#*                                      sation du logiciel.
#*
#************************************************************************
#*
#*   Script Name:       setup.py
#*
#*   Purpose:
#*      Distutils setup script for caom2repoClient
#*
#*   Functions:
#*
#*
#*
#****  C A N A D I A N   A S T R O N O M Y   D A T A   C E N T R E  *****
#************************************************************************
#-*/

# Use "distribute"
from setuptools import setup, find_packages
import sys

if sys.version_info[0] > 2:
    print 'The caom2repoClient package is only compatible with Python version 2.n'
    sys.exit(-1)

setup(name="caom2repoClient",
      version='trunk',
      description='Client application to the CAOM2 Repository',
      long_description='The caom2repo client allows the caller to get/put/update/delete CAOM2 Observations in a central repository.',
      url='http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2/repo',
      author='Canadian Astronomy Data Centre',
      author_email='cadc@hia.nrc.ca',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      provides=['caom2repoClient'],
      install_requires = ['distribute'],
      zip_safe=False
)

