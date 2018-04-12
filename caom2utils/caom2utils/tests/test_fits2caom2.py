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
from caom2utils import FitsParser, WcsParser, main_app, update_blueprint
from caom2utils import ObsBlueprint
from caom2utils.legacy import load_config

from caom2 import ObservationWriter, SimpleObservation, Algorithm
from caom2 import Artifact, ProductType, ReleaseType, ObservationIntentType
from caom2 import get_differences, obs_reader_writer, ObservationReader, Chunk
from lxml import etree

from mock import Mock, patch
from six import StringIO, BytesIO

import os
import sys

import pytest


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')
sample_file_4axes = os.path.join(TESTDATA_DIR, '4axes.fits')
sample_file_4axes_obs = os.path.join(TESTDATA_DIR, '4axes_obs.fits')
sample_file_time_axes = os.path.join(TESTDATA_DIR, 'time_axes.fits')
sample_file_4axes_uri = 'caom:CGPS/TEST/4axes_obs.fits'
java_config_file = os.path.join(TESTDATA_DIR, 'java.config')
override_file = os.path.join(TESTDATA_DIR, 'test.override')
test_override = os.path.join(TESTDATA_DIR, '4axes.override')
text_file = os.path.join(TESTDATA_DIR, 'help.txt')
text_override = os.path.join(TESTDATA_DIR, 'text.override')

# to execute only one test in the file set this var to True and comment
# out the skipif decorator of the test
single_test = False


class MyExitError(Exception):
    pass


EXPECTED_ENERGY_XML = \
    '''<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:energy>
    <caom2:axis>
      <caom2:axis>
        <caom2:ctype>VRAD</caom2:ctype>
        <caom2:cunit>m / s</caom2:cunit>
      </caom2:axis>
      <caom2:function>
        <caom2:naxis>1</caom2:naxis>
        <caom2:delta>-824.46002</caom2:delta>
        <caom2:refCoord>
          <caom2:pix>145.0</caom2:pix>
          <caom2:val>-60000.0</caom2:val>
        </caom2:refCoord>
      </caom2:function>
    </caom2:axis>
    <caom2:specsys>LSRK</caom2:specsys>
    <caom2:restfrq>1420406000.0</caom2:restfrq>
  </caom2:energy>
</caom2:import>
'''


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_energy():
    bp = ObsBlueprint(energy_axis=1)
    test_fitsparser = FitsParser(sample_file_4axes, bp)
    artifact = Artifact('ad:{}/{}'.format('TEST', sample_file_4axes),
                        ProductType.SCIENCE, ReleaseType.DATA)
    test_fitsparser.augment_artifact(artifact)
    energy = artifact.parts['0'].chunks[0].energy
    ex = _get_from_str_xml(EXPECTED_ENERGY_XML,
                           ObservationReader()._get_spectral_wcs, 'energy')
    result = get_differences(ex, energy)
    assert result is None


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact_energy_from_blueprint():
    test_blueprint = ObsBlueprint(energy_axis=1)
    test_blueprint.set('Chunk.energyAxis', 1)
    test_blueprint.set('Chunk.energy.specsys', 'LSRK')
    test_blueprint.set('Chunk.energy.axis.axis.ctype', 'VRAD')
    test_blueprint.set('Chunk.energy.axis.axis.cunit', 'm / s')
    test_blueprint.set('Chunk.energy.axis.function.refCoord.pix', '145.0')
    test_blueprint.set('Chunk.energy.axis.function.refCoord.val', '-60000.0')
    test_blueprint.set('Chunk.energy.axis.function.delta', '-824.46002')
    test_blueprint.set('Chunk.energy.axis.function.naxis', '1')
    test_fitsparser = FitsParser(sample_file_4axes, test_blueprint,
                                 uri='ad:TEST/test_blueprint')
    test_chunk = Chunk()
    test_fitsparser._try_energy_with_blueprint(test_chunk, 0)
    ex = _get_from_str_xml(EXPECTED_ENERGY_XML,
                           ObservationReader()._get_spectral_wcs,
                           'energy')
    result = get_differences(ex, test_chunk.energy)
    assert result is None


EXPECTED_POLARIZATION_XML = \
    '''<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
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
'''


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_polarization():
    test_fitsparser = FitsParser(sample_file_4axes,
                                 ObsBlueprint(polarization_axis=1))
    artifact = Artifact('ad:{}/{}'.format('TEST', sample_file_4axes),
                        ProductType.SCIENCE, ReleaseType.DATA)
    test_fitsparser.augment_artifact(artifact)
    polarization = artifact.parts['0'].chunks[0].polarization
    ex = _get_from_str_xml(EXPECTED_POLARIZATION_XML,
                           ObservationReader()._get_polarization_wcs,
                           'polarization')
    result = get_differences(ex, polarization)
    assert result is None, result


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact_polarization_from_blueprint():
    test_blueprint = ObsBlueprint(polarization_axis=1)
    test_blueprint.set('Chunk.polarizationAxis', '1')
    test_blueprint.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
    test_blueprint.set('Chunk.polarization.axis.function.refCoord.pix', '1.0')
    test_blueprint.set('Chunk.polarization.axis.function.refCoord.val', '1.0')
    test_blueprint.set('Chunk.polarization.axis.function.delta', '1.0')
    test_blueprint.set('Chunk.polarization.axis.function.naxis', '1')
    test_fitsparser = FitsParser(sample_file_4axes, test_blueprint,
                                 uri='test_parser')
    test_chunk = Chunk()
    test_fitsparser._try_polarization_with_blueprint(test_chunk, 0)
    ex = _get_from_str_xml(EXPECTED_POLARIZATION_XML,
                           ObservationReader()._get_polarization_wcs,
                           'polarization')
    result = get_differences(ex, test_chunk.polarization)
    assert result is None


EXPECTED_POSITION_XML = \
    '''<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
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
            <caom2:val>128.7499990027</caom2:val>
          </caom2:coord1>
          <caom2:coord2>
            <caom2:pix>513.0</caom2:pix>
            <caom2:val>-0.9999999922536</caom2:val>
          </caom2:coord2>
        </caom2:refCoord>
        <caom2:cd11>-0.004999999</caom2:cd11>
        <caom2:cd12>0.0</caom2:cd12>
        <caom2:cd21>0.0</caom2:cd21>
        <caom2:cd22>0.004999999</caom2:cd22>
      </caom2:function>
    </caom2:axis>
  </caom2:position>
</caom2:import>
'''


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact():
    test_blueprint = ObsBlueprint(position_axes=(1, 2))
    test_fitsparser = FitsParser(sample_file_4axes, test_blueprint)
    artifact = Artifact('ad:{}/{}'.format('TEST', sample_file_4axes),
                        ProductType.SCIENCE, ReleaseType.DATA)
    test_fitsparser.augment_artifact(artifact)
    assert artifact.parts is not None
    assert len(artifact.parts) == 1
    test_part = artifact.parts['0']
    test_chunk = test_part.chunks.pop()
    assert test_chunk is not None
    assert test_chunk.position is not None
    ex = _get_from_str_xml(EXPECTED_POSITION_XML,
                           ObservationReader()._get_spatial_wcs, 'position')
    result = get_differences(ex, test_chunk.position)
    assert result is None


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact_position_from_blueprint():
    test_blueprint = ObsBlueprint(position_axes=(1, 2))
    test_blueprint.set('Chunk.positionAxis1', '1')
    test_blueprint.set('Chunk.positionAxis2', '2')
    test_blueprint.set('Chunk.position.axis.axis1.ctype', 'GLON-CAR')
    test_blueprint.set('Chunk.position.axis.axis1.cunit', 'deg')
    test_blueprint.set('Chunk.position.axis.axis2.ctype', 'GLAT-CAR')
    test_blueprint.set('Chunk.position.axis.axis2.cunit', 'deg')
    test_blueprint.set('Chunk.position.axis.function.cd11', '-0.004999999')
    test_blueprint.set('Chunk.position.axis.function.cd12', '0.0')
    test_blueprint.set('Chunk.position.axis.function.cd21', '0.0')
    test_blueprint.set('Chunk.position.axis.function.cd22', '0.004999999')
    test_blueprint.set('Chunk.position.axis.function.dimension.naxis1', '1')
    test_blueprint.set('Chunk.position.axis.function.dimension.naxis2', '1')
    test_blueprint.set('Chunk.position.axis.range.start.coord1.pix', '513.0')
    test_blueprint.set('Chunk.position.axis.range.start.coord1.val',
                       '128.7499990027')
    test_blueprint.set('Chunk.position.axis.range.start.coord2.pix', '513.0')
    test_blueprint.set('Chunk.position.axis.range.start.coord2.val',
                       '-0.9999999922536')
    test_fitsparser = FitsParser(sample_file_4axes, test_blueprint,
                                 uri='test_parser')
    test_chunk = Chunk()
    test_fitsparser._try_position_with_blueprint(test_chunk, 0)
    ex = _get_from_str_xml(EXPECTED_POSITION_XML,
                           ObservationReader()._get_spatial_wcs, 'position')
    result = get_differences(ex, test_chunk.position)
    assert result is None


EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME = \
    '''<caom2:import xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3">
  <caom2:time>
    <caom2:axis>
      <caom2:axis>
        <caom2:ctype>TIME</caom2:ctype>
        <caom2:cunit>d</caom2:cunit>
      </caom2:axis>
      <caom2:error>
        <caom2:syser>1e-07</caom2:syser>
        <caom2:rnder>1e-07</caom2:rnder>
      </caom2:error>
      <caom2:function>
        <caom2:naxis>1</caom2:naxis>
        <caom2:delta>2.31481e-07</caom2:delta>
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
'''


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact_time():
    test_fitsparser = FitsParser(sample_file_time_axes,
                                 ObsBlueprint(time_axis=1))
    artifact = Artifact('ad:{}/{}'.format('TEST', sample_file_time_axes),
                        ProductType.SCIENCE, ReleaseType.DATA)
    test_fitsparser.augment_artifact(artifact)
    assert artifact.parts is not None
    assert len(artifact.parts) == 6
    test_part = artifact.parts['1']
    test_chunk = test_part.chunks.pop()
    assert test_chunk is not None
    assert test_chunk.time is not None
    ex = _get_from_str_xml(EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME,
                           ObservationReader()._get_temporal_wcs, 'time')
    result = get_differences(ex, test_chunk.time)
    assert result is None


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_artifact_time_from_blueprint():
    test_blueprint = ObsBlueprint(time_axis=1)
    test_blueprint.set('Chunk.timeAxis', '1')
    test_blueprint.set('Chunk.time.exposure', '0.02')
    test_blueprint.set('Chunk.time.resolution', '0.02')
    test_blueprint.set('Chunk.time.timesys', 'UTC')
    test_blueprint.set('Chunk.time.axis.axis.ctype', 'TIME')
    test_blueprint.set('Chunk.time.axis.axis.cunit', 'd')
    test_blueprint.set('Chunk.time.axis.error.syser', '1e-07')
    test_blueprint.set('Chunk.time.axis.error.rnder', '1e-07')
    test_blueprint.set('Chunk.time.axis.function.refCoord.pix', '0.5')
    test_blueprint.set('Chunk.time.axis.function.refCoord.val',
                       '56789.4298069')
    test_blueprint.set('Chunk.time.axis.function.delta', '2.31481e-07')
    test_blueprint.set('Chunk.time.axis.function.naxis', '1')
    test_fitsparser = FitsParser(sample_file_4axes, test_blueprint,
                                 uri='ad:TEST/test_blueprint')
    test_chunk = Chunk()
    test_fitsparser._try_time_with_blueprint(test_chunk, 0)
    ex = _get_from_str_xml(EXPECTED_CFHT_WIRCAM_RAW_GUIDE_CUBE_TIME,
                           ObservationReader()._get_temporal_wcs, 'time')
    result = get_differences(ex, test_chunk.time)
    assert result is None


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_get_wcs_values():
    w = get_test_wcs(sample_file_4axes)
    test_parser = WcsParser(get_test_header(sample_file_4axes)[0].header,
                            sample_file_4axes, 0)
    result = test_parser._sanitize(w.wcs.equinox)
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


def _get_from_str_xml(string_xml, get_func, element_tag):
    etree.register_namespace(
        'caom2', 'http://www.opencadc.org/caom2/xml/v2.3')
    parent_element = etree.fromstring(string_xml)
    ns = parent_element.nsmap['caom2']
    act_obj = get_func(element_tag, parent_element, ns, False)
    return act_obj


@pytest.mark.skipif(single_test, reason='Single test mode')
@patch('sys.exit', Mock(side_effect=[MyExitError, MyExitError, MyExitError,
                                     MyExitError, MyExitError,
                                     MyExitError]))
def test_help():
    """ Tests the helper displays for commands in main"""

    # expected helper messages
    with open(os.path.join(TESTDATA_DIR, 'bad_product_id.txt'), 'r') \
            as myfile:
        bad_product_id = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'missing_product_id.txt'), 'r') \
            as myfile:
        missing_product_id = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'too_few_arguments_help.txt'), 'r') \
            as myfile:
        too_few_arguments_usage = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'help.txt'), 'r') as myfile:
        usage = myfile.read()
    with open(os.path.join(TESTDATA_DIR, 'missing_observation_help.txt'), 'r')\
            as myfile:
        myfile.read()
    with open(os.path.join(TESTDATA_DIR,
                           'missing_positional_argument_help.txt'), 'r') \
            as myfile:
        myfile.read()

    # too few arguments error message when running python3
    with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2"]
        with pytest.raises(MyExitError):
            main_app()
        if stdout_mock.getvalue():
            assert (too_few_arguments_usage == stdout_mock.getvalue())

    # --help
    with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "-h"]
        with pytest.raises(MyExitError):
            main_app()
        assert (usage == stdout_mock.getvalue())

    # missing productID when plane count is wrong
    with patch('sys.stderr', new_callable=StringIO) as stderr_mock:
        bad_product_file = os.path.join(TESTDATA_DIR, 'bad_product_id.xml')
        sys.argv = ["fits2caom2", "--in", bad_product_file,
                    "ad:CGPS/CGPS_MA1_HI_line_image.fits"]
        with pytest.raises(MyExitError):
            main_app()
        assert bad_product_id == stderr_mock.getvalue()

    # missing productID when blueprint doesn't have one either
    with patch('sys.stderr', new_callable=StringIO) as stderr_mock:
        sys.argv = ["fits2caom2", "--observation", "test_collection_id",
                    "test_observation_id",
                    "ad:CGPS/CGPS_MA1_HI_line_image.fits"]
        with pytest.raises(MyExitError):
            main_app()
        assert missing_product_id == stderr_mock.getvalue()

    # missing required --observation
    """
    TODO: fix the tests
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "testProductID", "testpathto/testFileURI"]
        with pytest.raises(MyExitError):
            main_app()
        assert(missing_observation_usage == stdout_mock.getvalue())

    # missing positional argument
    with patch('sys.stderr', new_callable=StringIO) as stdout_mock:
        sys.argv = ["fits2caom2", "--observation", "testCollection",
                    "testObservationID", "testPathTo/testFileURI"]
        with pytest.raises(MyExitError):
            main_app()
        assert(missing_positional_argument_usage == stdout_mock.getvalue())
    """


EXPECTED_OBS_XML = """<?xml version='1.0' encoding='UTF-8'?>
<caom2:Observation""" + \
        """ xmlns:caom2="vos://cadc.nrc.ca!vospace/CADC/xml/CAOM/v2.0" """ + \
        """xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" """ + \
        """xsi:type="caom2:SimpleObservation" caom2:id="123">
  <caom2:collection>collection</caom2:collection>
  <caom2:observationID>MA1_DRAO-ST</caom2:observationID>
  <caom2:metaRelease>1999-01-01T00:00:00.000</caom2:metaRelease>
  <caom2:algorithm>
    <caom2:name>exposure</caom2:name>
  </caom2:algorithm>
  <caom2:intent>science</caom2:intent>
  <caom2:proposal>
    <caom2:id>HI-line</caom2:id>
  </caom2:proposal>
  <caom2:target>
    <caom2:name>CGPS Mosaic MA1</caom2:name>
    <caom2:standard>false</caom2:standard>
  </caom2:target>
  <caom2:telescope>
    <caom2:name>DRAO-ST</caom2:name>
    <caom2:geoLocationX>-2100330.87517</caom2:geoLocationX>
    <caom2:geoLocationY>-3694247.82445</caom2:geoLocationY>
    <caom2:geoLocationZ>4741018.33097</caom2:geoLocationZ>
  </caom2:telescope>
  <caom2:instrument>
    <caom2:name>DRAO-ST</caom2:name>
  </caom2:instrument>
  <caom2:planes>
    <caom2:plane caom2:id="123">
      <caom2:productID>HI-line</caom2:productID>
      <caom2:dataProductType>cube</caom2:dataProductType>
      <caom2:calibrationLevel>2</caom2:calibrationLevel>
      <caom2:provenance>
        <caom2:name>CGPS MOSAIC</caom2:name>
        <caom2:project>CGPS</caom2:project>
        <caom2:producer>CGPS Consortium</caom2:producer>
        <caom2:reference>http://dx.doi.org/10.1086/375301</caom2:reference>
        <caom2:lastExecuted>2000-10-16T00:00:00.000</caom2:lastExecuted>
      </caom2:provenance>
      <caom2:artifacts>
        <caom2:artifact caom2:id="123">
          <caom2:uri>caom:CGPS/TEST/4axes_obs.fits</caom2:uri>
          <caom2:productType>info</caom2:productType>
          <caom2:parts>
            <caom2:part caom2:id="123">
              <caom2:name>0</caom2:name>
              <caom2:chunks/>
            </caom2:part>
          </caom2:parts>
        </caom2:artifact>
      </caom2:artifacts>
    </caom2:plane>
  </caom2:planes>
</caom2:Observation>
"""


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_augment_observation():
    test_obs_blueprint = ObsBlueprint(position_axes=(1, 2))
    test_obs_blueprint.set('Observation.target.name', 'CGPS Mosaic MA1')
    test_obs_blueprint.set('Observation.target.standard', False)
    test_obs_blueprint.set('Observation.telescope.name', 'DRAO-ST')
    test_obs_blueprint.set('Observation.instrument.name', 'DRAO-ST')
    test_obs_blueprint.set('Observation.telescope.geoLocationX',
                           '-2100330.87517')
    test_obs_blueprint.set('Observation.telescope.geoLocationY',
                           '-3694247.82445')
    test_obs_blueprint.set('Observation.telescope.geoLocationZ',
                           '4741018.33097')

    test_obs_blueprint.set('Plane.dataProductType', 'cube')
    test_obs_blueprint.set('Artifact.productType', 'info')
    test_obs_blueprint.set('Artifact.releaseType', 'data')
    test_obs_blueprint.set('Plane.calibrationLevel', '2')
    test_fitsparser = FitsParser(sample_file_4axes_obs, test_obs_blueprint)
    test_fitsparser.blueprint = test_obs_blueprint
    test_obs = SimpleObservation('collection', 'MA1_DRAO-ST',
                                 Algorithm('exposure'))
    test_fitsparser.augment_observation(test_obs, sample_file_4axes_uri,
                                        product_id='HI-line')
    assert test_obs is not None
    assert test_obs.planes is not None
    assert len(test_obs.planes) == 1
    test_plane = test_obs.planes['HI-line']
    assert test_plane.artifacts is not None
    assert len(test_plane.artifacts) == 1
    test_artifact = test_plane.artifacts[sample_file_4axes_uri]
    assert test_artifact is not None
    test_part = test_artifact.parts['0']
    # remove the chunk bit, as it's part of other tests -
    # results in <caom2:chunks/> xml output
    test_part.chunks.pop()
    output = BytesIO()
    ow = ObservationWriter(False, False, "caom2",
                           obs_reader_writer.CAOM20_NAMESPACE)
    ow.write(test_obs, output)
    result = output.getvalue().decode('UTF-8')
    output.close()
    expected = _get_obs(EXPECTED_OBS_XML)
    actual = _get_obs(result)
    diff_result = get_differences(expected, actual, 'Observation')
    assert diff_result is None


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_get_from_list():
    test_fitsparser = FitsParser(sample_file_4axes)
    test_fitsparser.blueprint = ObsBlueprint()
    result = test_fitsparser._get_from_list('Observation.intent', 0,
                                            ObservationIntentType.SCIENCE)
    assert result == ObservationIntentType.SCIENCE


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_update_fits_headers():
    # The rules for the values:
    # all upper case - a FITS keyword
    # has an '{' or an '}' - a FITS keyword with an index
    #
    # The rules for the keys:
    # has a '.' - a config keyword
    # all upper case - a FITS keyword

    hdr1 = fits.Header()
    hdr2 = fits.Header()
    hdr3 = fits.Header()
    hdr4 = fits.Header()
    hdr5 = fits.Header()
    hdr6 = fits.Header()
    hdr7 = fits.Header()
    test_blueprint = ObsBlueprint()
    test_blueprint.configure_time_axis(3)

    test_parser = FitsParser(src=[hdr1, hdr2, hdr3, hdr4, hdr5, hdr6, hdr7])

    test_uri = 'ad:CFHT/1709071g.fits.gz'
    update_blueprint(test_blueprint, test_uri, config={}, defaults={},
                     overrides={})
    assert test_parser.blueprint._get('Observation.type') == \
        (['OBSTYPE'], None), 'unmodified blueprint'

    test_defaults = {'CTYPE1': 'RA---TAN',
                     'CTYPE2': 'DEC--TAN',
                     'CTYPE3': 'TIME',
                     'CTYPE4': 'WAVE',
                     'CDELT4': '1.2',
                     'CRVAL4': '32'}
    update_blueprint(test_blueprint, test_uri, config={},
                     defaults=test_defaults, overrides={})
    assert test_blueprint._get('Chunk.position.axis.axis1.ctype') == \
        (['CTYPE1'], 'RA---TAN'), 'default value assigned'
    assert test_blueprint._get('Chunk.position.axis.axis2.ctype') == \
        (['CTYPE2'], 'DEC--TAN'), 'default value assigned'
    assert test_blueprint._get('Chunk.time.axis.axis.ctype') ==  \
        (['CTYPE3'], 'TIME'), 'default value assigned, value all upper case'

    # print(test_parser.blueprint)

    test_defaults = {'CTYPE1': 'RA--TAN',
                     'CTYPE2': 'DEC--TAN',
                     'CTYPE3': 'TIME',
                     'plane.dataProductType': 'image',
                     'provenance.producer': 'CFHT',
                     'provenance.project': 'STANDARD PIPELINE'}
    test_config = load_config(java_config_file)
    test_overrides = load_config(override_file)
    update_blueprint(test_blueprint, test_uri, test_config,
                     test_defaults, test_overrides)
    assert test_blueprint._get('Plane.dataProductType') == \
        'image', 'default value assigned to configuration'
    assert test_blueprint._get('Plane.provenance.producer') == \
        (['ORIGIN'], 'CFHT'), \
        'default value assigned to configuration, all upper-case'
    assert test_blueprint._get('Plane.provenance.project') == \
        (['ADC_ARCH'], 'STANDARD PIPELINE'), \
        'default value assigned to configuration, with white-space'
    assert test_blueprint._get('Observation.type') == 'OBJECT', \
        'default value over-ridden, value all upper case'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val',
        0) == '210.551666667', 'override HDU 0'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val', 1) == \
        '210.551666667',         'override HDU 1'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val',
        2) == '210.508333333',         'override HDU 2'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val',
        3) == '210.898333333',         'override HDU 3'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val',
        4) == '210.942083333',         'override HDU 4'
    assert test_blueprint._get(
        'Chunk.position.axis.function.refCoord.coord1.val',
        5) == '0.000000000', 'override HDU 5'

    test_parser.blueprint = test_blueprint
    assert test_parser._headers[0][
               'CRVAL1'] == 210.551666667, 'override HDU 0'
    assert test_parser._headers[1][
               'CRVAL1'] == 210.551666667, 'override HDU 1'
    assert test_parser._headers[2][
               'CRVAL1'] == 210.508333333, 'override HDU 2'
    assert test_parser._headers[3][
               'CRVAL1'] == 210.898333333, 'override HDU 3'
    assert test_parser._headers[4][
               'CRVAL1'] == 210.942083333, 'override HDU 4'
    assert test_parser._headers[5]['CRVAL1'] == 0.000000000, 'override HDU 5'
    assert test_parser._headers[0][
               'CRVAL3'] == 56789.429806900000, 'override HDU 0'
    # this will fail because of CompositeObservation.members errors
    assert len(test_parser._errors) == 0, test_parser._errors


TEST_OVERRIDES = \
    {'obs.sequenceNumber': '1709071',
     'obs.intent': 'science',
     'obs.type': 'OBJECT',
     'target.standard': 'false',
     'proposal.project': '',
     'proposal.pi': 'Jean-Gabriel Cuby',
     'proposal.title': '',
     'plane.calibrationLevel': '1',
     'plane.dataRelease': '2015-02-01T00:00:00',
     'obs.metaRelease': '2014-05-12T10:18:55',
     'plane.metaRelease': '2014-05-12T10:18:55',
     'filtername': 'H.WC8201',
     'CRVAL4': '16310.000',
     'CDELT4': '2890.000',
     'resolvingPower': '5.64',
     'artifacts': {
         'ad:CFHT/1709071o.fits.fz': {
             0: {'artifact.productType': 'science',
                 'artifact.contentChecksum':
                     'md5:88bfd03471053a916067a4e6f80d332d',
                 'CRPIX3': '0.500000000000',
                 'CRVAL3': '56789.429806900000',
                 'CDELT3': '0.000173611111',
                 'time.resolution': '15.000000000000',
                 'time.exposure': '15.000000000000',
                 'NAXIS3': '3'}},
         'ad:CFHT/1709071g.fits.gz': {
             0: {'artifact.productType': 'auxiliary',
                 'artifact.contentChecksum':
                     'md5:47cdd15371f82893ed384dec96240ae2',
                 'CD1_1': '-0.000083333333',
                 'CD1_2': '0.000000000000',
                 'CD2_1': '0.000000000000',
                 'CD2_2': '0.000083333333',
                 'CRPIX3': '0.500000000000',
                 'CRVAL3': '56789.429806900000',
                 'CDELT3': '0.000000231481',
                 'time.resolution': '0.020000000000',
                 'time.exposure': '0.020000000000',
                 'NAXIS3': '1964',
                 'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '210.551666667',
                 'CRVAL2': '54.526222222'},
             1: {'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '210.551666667',
                 'CRVAL2': '54.526222222'},
             2: {'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '210.508333333',
                 'CRVAL2': '54.345555556'},
             3: {'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '210.898333333',
                 'CRVAL2': '54.341916667'},
             4: {'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '210.942083333',
                 'CRVAL2': '54.446805556'},
             5: {'CRPIX1': '7.00000000',
                 'CRPIX2': '7.00000000',
                 'CRVAL1': '0.000000000',
                 'CRVAL2': '0.000000000'},
             6: {'BITPIX': '0'}}
     }}


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_load_config_overrides():
    # cool override file content
    result = load_config(override_file)
    assert result == TEST_OVERRIDES


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_chunk_naxis():
    hdr1 = fits.Header()
    test_blueprint = ObsBlueprint()
    test_blueprint.configure_time_axis(3)
    test_uri = 'ad:CFHT/1709071g.fits.gz'
    test_defaults = {'CTYPE3': 'TIME'}
    test_config = {'Chunk.naxis': 'chunk.naxis'}
    test_overrides = {'chunk.naxis': '1'}
    update_blueprint(test_blueprint, test_uri, config=test_config,
                     defaults=test_defaults, overrides=test_overrides)
    assert test_blueprint._get('Chunk.naxis') == '1', 'default value assigned'
    FitsParser([hdr1], test_blueprint)
    assert hdr1['NAXIS'] == 1
    assert hdr1['ZNAXIS'] == 1


EXPECTED_FILE_SCHEME_XML = """<?xml version='1.0' encoding='UTF-8'?>
<caom2:Observation""" + \
    """ xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3" """ + \
    """xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" """ + \
    """xsi:type="caom2:SimpleObservation" """ + \
    """caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
  <caom2:collection>test_collection_id</caom2:collection>
  <caom2:observationID>test_observation_id</caom2:observationID>
  <caom2:metaRelease>1999-01-01T00:00:00.000</caom2:metaRelease>
  <caom2:algorithm>
    <caom2:name>exposure</caom2:name>
  </caom2:algorithm>
  <caom2:intent>science</caom2:intent>
  <caom2:target>
    <caom2:name>CGPS Mosaic MA1</caom2:name>
  </caom2:target>
  <caom2:instrument>
    <caom2:name>DRAO ST</caom2:name>
  </caom2:instrument>
  <caom2:planes>
    <caom2:plane caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
      <caom2:productID>test_product_id</caom2:productID>
      <caom2:dataProductType>cube</caom2:dataProductType>
      <caom2:calibrationLevel>2</caom2:calibrationLevel>
      <caom2:artifacts>
        <caom2:artifact caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
          <caom2:uri>file://""" + sample_file_4axes + """</caom2:uri>
          <caom2:productType>science</caom2:productType>
          <caom2:releaseType>data</caom2:releaseType>
          <caom2:contentType>application/octet-stream</caom2:contentType>
          <caom2:contentLength>11520</caom2:contentLength>
          <caom2:contentChecksum>md5:e6c08f3b8309f05a5a3330e27e3b44eb</caom2:contentChecksum>
          <caom2:parts>
            <caom2:part caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
              <caom2:name>0</caom2:name>
              <caom2:chunks>
                <caom2:chunk caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
                  <caom2:naxis>4</caom2:naxis>
                  <caom2:positionAxis1>1</caom2:positionAxis1>
                  <caom2:positionAxis2>2</caom2:positionAxis2>
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
                            <caom2:val>128.7499990027</caom2:val>
                          </caom2:coord1>
                          <caom2:coord2>
                            <caom2:pix>513.0</caom2:pix>
                            <caom2:val>-0.9999999922536</caom2:val>
                          </caom2:coord2>
                        </caom2:refCoord>
                        <caom2:cd11>-0.004999999</caom2:cd11>
                        <caom2:cd12>0.0</caom2:cd12>
                        <caom2:cd21>0.0</caom2:cd21>
                        <caom2:cd22>0.004999999</caom2:cd22>
                      </caom2:function>
                    </caom2:axis>
                  </caom2:position>
                </caom2:chunk>
              </caom2:chunks>
            </caom2:part>
          </caom2:parts>
        </caom2:artifact>
      </caom2:artifacts>
    </caom2:plane>
  </caom2:planes>
</caom2:Observation>
"""


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_file_scheme_uris():
    """ Tests that local files as URIs will be accepted and processed."""

    fname = 'file://{}'.format(sample_file_4axes)
    with patch('sys.stdout', new_callable=BytesIO) as stdout_mock:
        sys.argv = ['fits2caom2', '--observation', 'test_collection_id',
                    'test_observation_id', '--productID', 'test_product_id',
                    '--config', java_config_file, '--override', test_override,
                    fname]
        main_app()
        if stdout_mock.getvalue():
            expected = _get_obs(EXPECTED_FILE_SCHEME_XML)
            actual = _get_obs(stdout_mock.getvalue().decode('ascii'))
            result = get_differences(expected, actual, 'Observation')
            assert result is None


def _get_obs(from_xml_string):
    etree.parse = Mock(return_value=etree.ElementTree(
        etree.fromstring(from_xml_string.encode('ascii'))))
    obsr = obs_reader_writer.ObservationReader()
    obs = obsr.read(None)
    return obs


EXPECTED_GENERIC_PARSER_FILE_SCHEME_XML = """<?xml version='1.0' encoding='UTF-8'?>
<caom2:Observation""" + \
    """ xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3" """ + \
    """xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" """ + \
    """xsi:type="caom2:SimpleObservation" """ + \
    """caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
      <caom2:collection>test_collection_id</caom2:collection>
  <caom2:observationID>test_observation_id</caom2:observationID>
  <caom2:algorithm>
    <caom2:name>exposure</caom2:name>
  </caom2:algorithm>
  <caom2:planes>
    <caom2:plane caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
      <caom2:productID>test_product_id</caom2:productID>
      <caom2:artifacts>
        <caom2:artifact caom2:id="00000000-0000-0000-8ffa-fc0a5a8759df">
          <caom2:productType>thumbnail</caom2:productType>
          <caom2:releaseType>data</caom2:releaseType>
          <caom2:contentType>text/plain</caom2:contentType>
          <caom2:contentLength>2832</caom2:contentLength>
          <caom2:contentChecksum>md5:e6c08f3b8309f05a5a3330e27e3b44eb</caom2:contentChecksum>
          <caom2:uri>file://""" + text_file + """</caom2:uri>
        </caom2:artifact>
      </caom2:artifacts>
    </caom2:plane>
  </caom2:planes>
</caom2:Observation>
"""


@pytest.mark.skipif(single_test, reason='Single test mode')
def test_generic_parser():
    """ Tests that GenericParser will be created."""

    fname = 'file://{}'.format(text_file)
    with patch('sys.stdout', new_callable=BytesIO) as stdout_mock:
        sys.argv = ['fits2caom2', '--local', fname,
                    '--observation', 'test_collection_id',
                    'test_observation_id', '--productID', 'test_product_id',
                    '--config', java_config_file, '--override', text_override,
                    fname]
        main_app()
        if stdout_mock.getvalue():
            expected = _get_obs(EXPECTED_GENERIC_PARSER_FILE_SCHEME_XML)
            actual = _get_obs(stdout_mock.getvalue().decode('ascii'))
            result = get_differences(expected, actual, 'Observation')
            assert result is None
