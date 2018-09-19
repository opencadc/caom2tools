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

import logging

from astropy.io import fits
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta

from datetime import timedelta as dt_timedelta
from datetime import datetime as dt_datetime
from time import strptime as dt_strptime

from caom2 import Interval as caom_Interval
from caom2 import Time as caom_Time
from caom2 import shape as caom_shape

from caom2pipe import manage_composable as mc


__all__ = ['convert_time', 'get_datetime', 'build_plane_time', 'get_location',
           'get_timedelta_in_s', 'make_headers_from_string']


def find_time_bounds(headers):
    """Given an observation date, and time exposure length, calculate the
    mjd_start and mjd_end time for those values."""
    logging.debug('Begin find_time_bounds.')
    date = headers[0].get('DATE-OBS')
    exposure = headers[0].get('TEXP')
    return convert_time(date, exposure)


def convert_time(start_time, exposure):
    """Convert a start time and exposure length into an mjd_start and mjd_end
    time."""
    logging.debug('Begin convert_time.')
    if start_time is not None and exposure is not None:
        logging.debug(
            'Use date {} and exposure {} to convert time.'.format(start_time,
                                                                  exposure))
        if type(start_time) is float:
            t_start = Time(start_time, format='mjd')
        else:
            t_start = Time(start_time)
        dt = TimeDelta(exposure, format='sec')
        t_end = t_start + dt
        t_start.format = 'mjd'
        t_end.format = 'mjd'
        mjd_start = t_start.value
        mjd_end = t_end.value
        logging.debug('End convert_time mjd start {} mjd end {} .'.format(
            mjd_start, mjd_end))
        return mjd_start, mjd_end
    return None, None


def get_datetime(from_value):
    """
    Ensure datetime values are in MJD. This is meant to handle any odd formats
    that telescopes have for datetime values.

    Relies on astropy, until astropy fails.

    :param from_value:
    :return: datetime instance
    """
    if from_value is not None:
        try:
            result = Time(from_value)
            result.format = 'mjd'
            return result
        except ValueError:
            try:
                # VLASS has a format astropy fails to understand
                # from datetime import datetime
                result = Time(
                    dt_datetime.strptime(from_value, '%H:%M:%S'))
                result.format = 'mjd'
                return result
            except ValueError:
                logging.error('Cannot parse datetime {}'.format(from_value))
                return None
    else:
        return None


def get_location(latitude, longitude, elevation):
    """The CAOM model expects the telescope location to be in geocentric
    coordinates. Rely on astropy to do the conversion."""
    result = EarthLocation.from_geodetic(
        longitude, latitude, elevation, 'WGS84')
    return result.x.value, result.y.value, result.z.value


def build_plane_time(start_date, end_date, exposure_time):
    """Calculate the plane-level bounding box for time."""
    start_date.format = 'mjd'
    end_date.format = 'mjd'
    time_bounds = caom_Interval(mc.to_float(start_date.value),
                                mc.to_float(end_date.value),
                                samples=[
                                    caom_shape.SubInterval(
                                        mc.to_float(start_date.value),
                                        mc.to_float(end_date.value))])
    return caom_Time(bounds=time_bounds,
                     dimension=1,
                     resolution=exposure_time.to('second').value,
                     sample_size=exposure_time.to('day').value,
                     exposure=exposure_time.to('second').value)


def get_timedelta_in_s(from_value):
    """
    :param from_value: a string representing time in H:M:S
    :return: the value as a timedelta, in seconds
    """
    temp = dt_strptime(from_value, '%H:%M:%S')
    td = dt_timedelta(
        hours=temp.tm_hour, minutes=temp.tm_min, seconds=temp.tm_sec)
    return td.seconds


def make_headers_from_string(fits_header):
    """Create a list of fits.Header instances from a string.
    ":param fits_header a string of keyword/value pairs"""
    delim = '\nEND'
    extensions = \
        [e + delim for e in fits_header.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers
