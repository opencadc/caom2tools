# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2016.                            (c) 2016.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.io import fits
from astropy.wcs import WCS as awcs
from caom2utils import FitsParser, WcsParser, main_app

from caom2 import ObservationWriter
from caom2 import Artifact, ProductType, ReleaseType
from lxml import etree

from mock import Mock, patch
from six import StringIO
from caom2utils import fits2caom2

from io import BytesIO
import os
import sys
import tempfile
import re

import pytest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')
expected_cgps_obs = os.path.join(TESTDATA_DIR, 'cgps.xml')
expected_local_cgps_obs = os.path.join(TESTDATA_DIR, 'cgps_local.xml')
sample_file_4axes = os.path.join(TESTDATA_DIR, '4axes.fits')
sample_file_time_axes = os.path.join(TESTDATA_DIR, 'time_axes.fits')

def test_fits2caom2():
    # test fits2caom2 on a known existing CGPS file
    expected = open(expected_cgps_obs).read()
    with patch('sys.stdout', new_callable=BytesIO) as stdout_mock:
        sys.argv = ('fits2caom2 --observation TEST myOBS myplane '
                    'ad:CGPS/CGPS_MA1_HI_line_image').split()
        fits2caom2.main_app()
    actual = stdout_mock.getvalue().decode('ascii')
    _cmp(expected, actual)

    # repeat the test when the observation is saved
    temp = tempfile.NamedTemporaryFile();
    sys.argv = ('fits2caom2 --observation TEST myOBS -o {} myplane '
                'ad:CGPS/CGPS_MA1_HI_line_image'.format(temp.name)).split()
    fits2caom2.main_app()

    actual = open(temp.name).read()
    _cmp(expected, actual)

    # test fits2caom2 on a known existing but now local CGPS file
    expected = open(expected_local_cgps_obs).read()
    with patch('sys.stdout', new_callable=BytesIO) as stdout_mock:
        sys.argv = ('fits2caom2 --local {} --observation TEST myOBS myplane '
                    'ad:CGPS/CGPS_MA1_HI_line_image').format(
            sample_file_4axes).split()
        fits2caom2.main_app()
    actual = stdout_mock.getvalue().decode('ascii')
    _cmp(expected, actual)

    # repeat the test when the observation is saved
    temp = tempfile.NamedTemporaryFile();
    sys.argv = ('fits2caom2 --observation TEST myOBS --local {} -o {} myplane '
                'ad:CGPS/CGPS_MA1_HI_line_image'.format(
        sample_file_4axes, temp.name)).split()
    fits2caom2.main_app()

    actual = open(temp.name).read()
    _cmp(expected, actual)


def _cmp(expected_obs_xml, actual_obs_xml):
    """
    Textual comparison of the xml representation of 2 observations ignoring
    the UUIDs.
    :param expected_obs_xml:
    :param actual_obs_xml:
    :return:
    """
    expected = re.sub(r'caom2:id=".*"', 'caom2:id=""', expected_obs_xml)
    actual = re.sub(r'caom2:id=".*"', 'caom2:id=""', actual_obs_xml)

    assert expected == actual