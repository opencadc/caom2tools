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
from caom2utils import augment_artifact
from caom2utils import avoid_nan

from caom2 import ObservationWriter
from lxml import etree

import os
import sys

import pytest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')

# note - added significant digits for the values of coord1 val, coord2 val cd11, cd22
#
EXPECTED_CGPS_POSITION_XML = \
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
          <caom2:naxis1>1024</caom2:naxis1>
          <caom2:naxis2>1024</caom2:naxis2>
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


@pytest.mark.parametrize('test_input', ['CGPS_MA1_HI_line_image.fits'])
def test_augment_artifact(test_input):
    test_file = os.path.join(TESTDATA_DIR, test_input)
    test_artifact = augment_artifact(None, test_file)
    assert test_artifact is not None
    assert test_artifact.parts is not None
    assert len(test_artifact.parts) == 1
    test_part = test_artifact.parts['0']
    test_chunk = test_part.chunks.pop()
    assert test_chunk is not None
    assert test_chunk.position is not None

    etree.register_namespace('caom2', 'http://www.opencadc.org/caom2/xml/v2.3')
    parent_element = etree.Element('{http://www.opencadc.org/caom2/xml/v2.3}import')
    ow = ObservationWriter()
    ow._add_spatial_wcs_element(test_chunk.position, parent_element)
    tree = etree.ElementTree(parent_element)
    result = etree.tostring(tree, encoding='unicode', pretty_print=True)
    assert result == EXPECTED_CGPS_POSITION_XML, result


@pytest.mark.parametrize('test_file', ['CGPS_MA1_HI_line_image.fits'])
def test_get_wcs_values(test_file):
    w = get_test_wcs(test_file)
    result = avoid_nan(w.wcs, 'equinox')
    assert result is None
    result = getattr(w, '_naxis1')
    assert result == 1024
    assert w.wcs.has_cd() is False


def get_test_header(test_file):
    test_input = os.path.join(TESTDATA_DIR, test_file)
    hdulist = fits.open(test_input)
    # print(repr(hdulist[0].header))
    hdulist.close();
    return hdulist


def get_test_wcs(test_file):
    hdu = get_test_header(test_file)
    wcs = awcs(hdu[0])
    return wcs