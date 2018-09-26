# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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

import math
import pytest
import sys

from astropy.io import fits

from caom2pipe import astro_composable as ac


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_convert_time():
    hdr1 = fits.Header()
    mjd_start, mjd_end = ac.find_time_bounds([hdr1])
    assert mjd_start is None
    assert mjd_end is None

    hdr1['DATE-OBS'] = '2012-09-03T01:04:44'
    hdr1['TEXP'] = 20.000
    mjd_start, mjd_end = ac.find_time_bounds([hdr1])
    assert mjd_start is not None
    assert mjd_end is not None
    assert math.isclose(mjd_start, 56173.044953703706), mjd_start
    assert math.isclose(mjd_end, 56173.04518518518), mjd_end


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_datetime():
    result = ac.get_datetime('2006-12-12T12:12:12')
    assert result is not None
    assert result == '2006-12-12 12:12:12.000'

    result = ac.get_datetime('2006-12-12 12:12:12.001')
    assert result is not None
    assert result == '2006-12-12 12:12:12.001'

    result = ac.get_datetime('2006-12-12')
    assert result is not None
    assert result == '2006-12-12 00:00:00.000'

    # a format that is not understood
    result = ac.get_datetime('16-Dec-12T01:23:45')
    assert result is None

    result = ac.get_datetime(None)
    assert result is None


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_location():
    x, y, z = ac.get_location(21.0, -32.0, 12)
    assert x == 5051887.288718968, x
    assert y == -3156769.536020791, y
    assert z == 2271399.319625149, z


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_build_plane_time():
    start = ac.get_datetime('2012-09-03T01:04:44')
    end = ac.get_datetime('2012-09-03T03:04:44')
    exposure = end - start
    result = ac.build_plane_time(start, end, exposure)
    assert result is not None, 'expected a value'
    assert result.bounds is not None, 'expected a bounds value'
    assert result.exposure == 7199.999999999994, 'wrong exposure value'


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_time_delta_in_s():
    result = ac.get_timedelta_in_s('0:06:41')
    assert result is not None
    assert result == 401, 'wrong value returned'
