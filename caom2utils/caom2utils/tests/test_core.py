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

from builtins import str

from astropy.io import fits
from astropy.wcs import WCS as awcs
from caom2utils import FitsParser

from caom2 import ObservationWriter
from lxml import etree

from mock import Mock, patch
from six import StringIO

import os
import sys

import pytest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')
sample_file_4axes = os.path.join(TESTDATA_DIR, '4axes.fits')
sample_file_time_axes = os.path.join(TESTDATA_DIR, 'time_axes.fits')


class MyExitError(Exception):
    pass


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_augment_plarization(test_file):
    test_fitsparser = FitsParser(test_file)
    artifact = test_fitsparser.augment_artifact()


EXPECTED_ENERGY_XML = '''<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:energy>
    <caom2:axis>
      <caom2:axis>
        <caom2:ctype>VRAD</caom2:ctype>
        <caom2:cunit>m / s</caom2:cunit>
      </caom2:axis>
      <caom2:function>
        <caom2:naxis>1</caom2:naxis>
        <caom2:delta>-824.46001999999999</caom2:delta>
        <caom2:refCoord>
          <caom2:pix>145.0</caom2:pix>
          <caom2:val>-60000.0</caom2:val>
        </caom2:refCoord>
      </caom2:function>
    </caom2:axis>
    <caom2:specsys>LSRK</caom2:specsys>
    <caom2:restfrq>1420406000.0</caom2:restfrq>
    <caom2:restwav>0.0</caom2:restwav>
  </caom2:energy>
</caom2:import>
'''


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_augment_energy(test_file):
    test_fitsparser = FitsParser(test_file)
    artifact = test_fitsparser.augment_artifact()
    energy = artifact.parts['0'].chunks[0].energy
    energy.bandpassName = '21 cm' # user set attribute
    check_xml(ObservationWriter()._add_spectral_wcs_element, energy, EXPECTED_ENERGY_XML)


EXPECTED_POLARIZATION_XML = \
    """<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:polarization>
    <caom2:axis>
      <caom2:axis>
        <caom2:ctype>STOKES</caom2:ctype>
      </caom2:axis>
      <caom2:function>
        <caom2:naxis>1</caom2:naxis>
        <caom2:delta>1.0</caom2:delta>
        <caom2:refCoord>
          <caom2:pix>1.0</caom2:pix>
          <caom2:val>1.0</caom2:val>
        </caom2:refCoord>
      </caom2:function>
    </caom2:axis>
  </caom2:polarization>
</caom2:import>
"""


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_augment_polarization(test_file):
    test_fitsparser = FitsParser(test_file)
    artifact = test_fitsparser.augment_artifact()
    polarization = artifact.parts['0'].chunks[0].polarization
    check_xml(ObservationWriter()._add_polarization_wcs_element, polarization, EXPECTED_POLARIZATION_XML)


EXPECTED_POSITION_XML = \
    """<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:position>
    <caom2:axis>
      <caom2:axis1>
        <caom2:ctype>GLON-CAR</caom2:ctype>
        <caom2:cunit>deg</caom2:cunit>
      </caom2:axis1>
      <caom2:axis2>
        <caom2:ctype>GLAT-CAR</caom2:ctype>
        <caom2:cunit>deg</caom2:cunit>
      </caom2:axis2>
      <caom2:function>
        <caom2:dimension>
          <caom2:naxis1>1</caom2:naxis1>
          <caom2:naxis2>1</caom2:naxis2>
        </caom2:dimension>
        <caom2:refCoord>
          <caom2:coord1>
            <caom2:pix>513.0</caom2:pix>
            <caom2:val>128.74999900270001</caom2:val>
          </caom2:coord1>
          <caom2:coord2>
            <caom2:pix>513.0</caom2:pix>
            <caom2:val>-0.99999999225360003</caom2:val>
          </caom2:coord2>
        </caom2:refCoord>
        <caom2:cd11>-0.0049999989999999998</caom2:cd11>
        <caom2:cd12>0.0</caom2:cd12>
        <caom2:cd21>0.0</caom2:cd21>
        <caom2:cd22>0.0049999989999999998</caom2:cd22>
      </caom2:function>
    </caom2:axis>
  </caom2:position>
</caom2:import>
"""


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_augment_artifact(test_file):
    test_fitsparser = FitsParser(test_file)
    test_artifact = test_fitsparser.augment_artifact()
    assert test_artifact is not None
    assert test_artifact.parts is not None
    assert len(test_artifact.parts) == 1
    test_part = test_artifact.parts['0']
    test_chunk = test_part.chunks.pop()
    assert test_chunk is not None
    assert test_chunk.position is not None
    check_xml(ObservationWriter()._add_spatial_wcs_element, test_chunk.position, EXPECTED_POSITION_XML)


EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME = \
    """<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:time>
    <caom2:axis>
      <caom2:axis>
        <caom2:ctype>TIME</caom2:ctype>
        <caom2:cunit>d</caom2:cunit>
      </caom2:axis>
      <caom2:error>
        <caom2:syser>9.9999999999999995e-08</caom2:syser>
        <caom2:rnder>9.9999999999999995e-08</caom2:rnder>
      </caom2:error>
      <caom2:function>
        <caom2:naxis>1</caom2:naxis>
        <caom2:delta>2.3148100000000001e-07</caom2:delta>
        <caom2:refCoord>
          <caom2:pix>0.5</caom2:pix>
          <caom2:val>56789.4298069</caom2:val>
        </caom2:refCoord>
      </caom2:function>
    </caom2:axis>
    <caom2:timesys>UTC</caom2:timesys>
    <caom2:exposure>0.02</caom2:exposure>
    <caom2:resolution>0.02</caom2:resolution>
  </caom2:time>
</caom2:import>
"""


# @pytest.mark.parametrize('test_file, expected', [(os.path.join(TESTDATA_DIR, '1709071g.fits'), EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME_XML),
#                                                  (os.path.join(TESTDATA_DIR, '1916216i.fits'), EXPECTED_CFHT_TIME_XML),
#                                                  (os.path.join(TESTDATA_DIR, '1916216o.fits'), EXPECTED_CFHT_TIME_XML),
#                                                  (os.path.join(TESTDATA_DIR, '1916216p.fits'), EXPECTED_CFHT_TIME_XML),
#                                                  (os.path.join(TESTDATA_DIR, '2136164p.fits'), EXPECTED_CFHT_TIME_XML)])

# @pytest.mark.parametrize('test_file, expected', [(os.path.join(TESTDATA_DIR, '1709071g.fits'),
#                                                   EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME),
#                                                  (os.path.join(TESTDATA_DIR, 'oc6n01010_x1d.fits'),
#                                                    EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME),
#                                                  (os.path.join(TESTDATA_DIR, 'jcmth20100515_00048_02_rimg_pro_000.fits'),
#                                                   EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME)])

@pytest.mark.parametrize('test_file, expected', [(sample_file_time_axes,
                                                  EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME)])
def test_augment_artifact_time(test_file, expected):
    test_fitsparser = FitsParser(test_file)
    test_artifact = test_fitsparser.augment_artifact()
    assert test_artifact is not None
    assert test_artifact.parts is not None
    assert len(test_artifact.parts) == 6
    test_part = test_artifact.parts['1']
    test_chunk = test_part.chunks.pop()
    assert test_chunk is not None
    assert test_chunk.position is not None
    check_xml(ObservationWriter()._add_temporal_wcs_element, test_chunk.time, expected)


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_get_wcs_values(test_file):
    w = get_test_wcs(test_file)
    test_fitsparser = FitsParser(test_file)
    result = test_fitsparser.fix_value(w.wcs.equinox)
    assert result is None
    result = getattr(w, '_naxis1')
    assert result == 1
    assert w.wcs.has_cd() is False


def get_test_header(test_file):
    test_input = os.path.join(TESTDATA_DIR, test_file)
    hdulist = fits.open(test_input)
    hdulist.close()
    return hdulist


def get_test_wcs(test_file):
    hdu = get_test_header(test_file)
    wcs = awcs(hdu[0])
    return wcs


def check_xml(xml_func, test_wcs, expected):
    etree.register_namespace('caom2', 'http://www.opencadc.org/caom2/xml/v2.3')
    parent_element = etree.Element('{http://www.opencadc.org/caom2/xml/v2.3}import')
    xml_func(test_wcs, parent_element)
    tree = etree.ElementTree(parent_element)
    result = etree.tostring(tree, encoding='unicode', pretty_print=True)
    assert result == expected, result


@patch('sys.exit', Mock(side_effect=[MyExitError, MyExitError, MyExitError,
                                     MyExitError, MyExitError,
                                     MyExitError]))
@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_help(test_file):
    """ Tests the helper displays for commands in main"""

    # expected helper messages
    with open(os.path.join(TESTDATA_DIR, 'help.txt'), 'r') as myfile:
        usage = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'missing_observation_help.txt'), 'r') as myfile:
        missing_observation_usage = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'missing_positional_argument_help.txt'), 'r') as myfile:
        missing_positional_argument_usage = myfile.read()

    # --help
    with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "-h"]
        with pytest.raises(MyExitError):
            test_fitsparser = FitsParser(test_file)
            test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(usage == stdout_mock.getvalue())

    # missing required --observation
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "testProductID", "testpathto/testFileURI"]
        with pytest.raises(MyExitError):
            test_fitsparser = FitsParser(test_file)
            test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(missing_observation_usage == stdout_mock.getvalue())

    # missing positional argument
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--observation", "testCollection", "testObservationID", "testPathTo/testFileURI"]
        with pytest.raises(MyExitError):
            test_fitsparser = FitsParser(test_file)
            test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(missing_positional_argument_usage == stdout_mock.getvalue())


@pytest.mark.parametrize('test_file', [sample_file_4axes])
def test_valid_arguments(test_file):
    """ Tests the parser with valid commands in main"""

    # --in
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--in", "pathTo/inObsXML",
                    "productID", "pathTo/testFileURI1", "pathTo/testFileURI2"]
        test_fitsparser = FitsParser(test_file)
        test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(not stdout_mock.getvalue())

    # --in and --out
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--in", "pathTo/inObsXML", "--out", "pathTo/outObsXML",
                    "productID", "pathTo/testFileURI1", "pathTo/testFileURI2"]
        test_fitsparser = FitsParser(test_file)
        test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(not stdout_mock.getvalue())

    # --observation
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--observation", "testCollection", "testObservationID",
                    "productID", "pathTo/testFileURI1", "pathTo/testFileURI2"]
        test_fitsparser = FitsParser(test_file)
        test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(not stdout_mock.getvalue())

    # --observation and --out
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--observation", "testCollection", "testObservationID", "--out", "pathTo/outObsXML"
                    "productID", "pathTo/testFileURI1", "pathTo/testFileURI2"]
        test_fitsparser = FitsParser(test_file)
        test_fitsparser.main_app()
        help_out = stdout_mock.getvalue()
        assert(not stdout_mock.getvalue())


