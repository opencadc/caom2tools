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

"""
WCS Validation Utilities

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from six.moves import range
from caom2 import shape

__all__ = ['TimeUtil', 'EnergyUtil', 'ORIGIN']

# TODO both these are very bad, implement more sensibly
TARGET_TIMESYS = "UTC"
TARGET_CTYPE = "TIME"
TARGET_CUNIT = "d"
ORIGIN = 0


# CoordFunction1D, double/float
def pix2val(function, pix):
    ref_pix = float(function.ref_coord.pix)
    return function.ref_coord.val + function.delta * (pix - ref_pix)


class TimeUtil:
    def __init__(self):
        pass

    @staticmethod
    def range1d_to_interval(temporal_wcs, range_1d):
        TimeUtil.validate_wcs(temporal_wcs)

        # TODO: (comment pulled from Java code):
        # if mjdref has a value then the units of axis values could be any time
        # units, like days, hours, minutes, seconds, and smaller
        # since they are offsets from mjdref
        a = range_1d.start.val
        b = range_1d.end.val
        if b < a:
            raise ValueError("range.end not >= range.start in Temporal WCS")

        if temporal_wcs.mjdref is not None:
            a += float(temporal_wcs.mjdref)
            b += float(temporal_wcs.mjdref)

        return shape.SubInterval(min(a, b), max(a, b))

    @staticmethod
    def function1d_to_interval(temporal_wcs, function_1d):
        try:
            TimeUtil.validate_wcs(temporal_wcs)
            # // TODO: (comment pulled from Java code):
            # if mjdref has a value then the units of axis values could be
            # any time
            # // units, like days, hours, minutes, seconds, and smaller
            # // since they are offsets from mjdref

            p1 = float(0.5)
            p2 = float(function_1d.naxis + 0.5)
            a = pix2val(function_1d, p1)
            b = pix2val(function_1d, p2)
            if function_1d.delta < 0.0:
                raise ValueError(
                    '{} delta must be greater than 0.0'.format(function_1d))

            if temporal_wcs.mjdref is not None:
                a += float(temporal_wcs.mjdref)
                b += float(temporal_wcs.mjdref)

            return shape.SubInterval(min(a, b), max(a, b))

        except Exception as ex:
            raise ValueError(
                "Invalid function in Temporal WCS: {}".format(repr(ex)))

    @staticmethod
    def validate_wcs(temporal_wcs):
        ctype = temporal_wcs.axis.axis.ctype
        sb = ""
        if ctype == TARGET_CTYPE \
                and (temporal_wcs.timesys is None
                     or temporal_wcs.timesys == TARGET_TIMESYS):
            pass
        elif ctype == TARGET_TIMESYS and temporal_wcs.timesys is None:
            pass
        else:
            sb = "unexpected TIMESYS, CTYPE: {},{}".format(
                temporal_wcs.timesys, ctype)

        cunit = temporal_wcs.axis.axis.cunit
        if TARGET_CUNIT != cunit:
            sb = sb + "unexpected CUNIT: {}".format(cunit)

        if len(sb) > 0:
            raise ValueError(sb)


class EnergyUtil:
    def __init__(self):
        pass

    @staticmethod
    def range1d_to_interval(range_1d):
        a = float(range_1d.start.val)
        b = float(range_1d.end.val)
        #  The energy converter work done in the Java code is skipped here.
        #  Doing it here introduced some precision errors which lead to false
        #  invalids. Ignoring the units for validation sounds like it's ok,
        #  as long as the same values come out of the p2s, s2p calculations
        #  in the main validator code. It's assumed that doing the
        #  conversions in the native units is sufficient, as long as the
        #  values are the same after p2s -> s2p is done.

        return shape.SubInterval(min(a, b), max(a, b))

    @staticmethod
    def function1d_to_interval(temporal_wcs):
        naxis = temporal_wcs.axis.function.naxis
        p1 = 0.5
        p2 = naxis + 0.5
        return shape.SubInterval(p1, p2)


class PolarizationWcsUtil:
    def _get_range(self, from_range):
        if from_range is not None:
            lb = int(round(from_range.start.val))
            ub = int(round(from_range.end.val))
            return range(lb, ub+1)
        return None

    @staticmethod
    def get_keys(range):
        """
        Examines the lower bound (lb) and upper bound (ub) of
        PolarizationWCS.axis.range and returns range(lb, ub+1) if a range is
        defined, else returns None. Since Python range iterates over [lb,ub),
        the returned range is ub+1 to ensure that ub is included in the
        range iteration.
        """
        return PolarizationWcsUtil()._get_range(range)

    @staticmethod
    def get_ranges_from_bounds(bounds):
        """
        Examines the ranges in PolarizationWCS.axis.bounds and returns the
        list of ranges in the bounds if the bounds is defined, else returns
        an empty list. The upper bound of each range is incremented by 1
        (refer to comments for get_range above)
        """
        ranges = []
        if bounds is not None:
            samples = bounds.samples
            if samples is not None:
                for sample in samples:
                    ranges.append(PolarizationWcsUtil()._get_range(sample))
        return ranges

    @staticmethod
    def get_range_from_function(function):
        """
        Examines the ranges in PolarizationWCS.axis.function and returns a
        range from 1 to Naxis+1 if the function is defined, else returns
        None. The upper bound of the range is incremented by 1 (refer to
        comments for get_range above)
        """
        if function is not None:
            if function.naxis >= 1:
                return range(1, function.naxis + 1)
            else:
                raise ValueError(
                    'Invalid naxis value: {}'.format(function.naxis))
        return None
