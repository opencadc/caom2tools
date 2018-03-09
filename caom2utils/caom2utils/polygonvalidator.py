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

from __future__ import (absolute_import, print_function, unicode_literals)

import six

import numpy as np
from spherical_geometry import polygon


__all__ = ['validate_polygon']


def validate_polygon(poly):
    """
    Performs validation on the polygon provided. The points in the polygon
    will be validated for closure and self segment intersection. The samples
    in the polygon will be validated using validate_multipolygon().

    for a closed polygon, we must have:
        points[0].cval1 == points[-1].cval1
        points[0].cval2 == points[-1].cval2

    An AssertionError is thrown if the object does not represent a polygon

    :param poly: Polygon to be validated
    """
    points = poly.points
    if points:
        if len(points) < 3:
            # points in a polygon is not required to form a closed polygon,
            # hence min 3 points
            raise AssertionError('invalid polygon: {} points (min 3)'.format(
                len(points)))

        cval1s = []
        cval2s = []
        for i in range(len(points)):
            cval1s.append(points[i].cval1)
            cval2s.append(points[i].cval2)
        if (cval1s[-1] == 0) and (cval2s[-1] == 0):
            # (0,0) represents the last point in a closed polygon
            # SphericalPolygon requires point[0]==point[-1] for closed polygons
            cval1s[-1] = cval1s[0]
            cval2s[-1] = cval2s[0]
        elif (cval1s[-1] != cval1s[0]) or (cval2s[-1] != cval2s[0]):
            cval1s.append(points[0].cval1)
            cval2s.append(points[0].cval2)

        # validate the polygons in the multipolygon
        _validate_self_intersection_and_direction(cval1s, cval2s)

    # validate the samples
    if poly.samples is not None:
        validate_multipolygon(poly.samples)


# Stub for now
def validate_multipolygon(mpoly):
    pass


def _validate_self_intersection(spolygon):
    """
    Verifies that the polygon does not contain self-intersecting segments.
    Note: The (new version) spherical-geometry package has been updated to
    detect and fix self-segment intersection. However at time of writing
    of this comment, this new feature has not been released. The
    (current version) currently released version does not detect
    self-segment intersection.

    An AssertionError is thrown if the polygon contains self-intersecting
    segments.
    """
    if hasattr(spolygon, '__len__'):
        # working with new version of spherical-geometry
        if len(spolygon) > 1:
            raise AssertionError(
                'Polygon contains self intersecting segments')
    else:
        # working with current version of spherical-geometry
        for p in spolygon.iter_polygons_flat():
            for i in range(len(p._points) - 3):
                A = p._points[i]
                B = p._points[i + 1]
                C = p._points[2:-2] if i == 0 else p._points[i + 2:-1]
                D = p._points[3:-1] if i == 0 else p._points[i + 3:]
                if np.any(polygon.great_circle_arc.intersects(A, B, C, D)):
                    raise AssertionError(
                        'Polygon contains self intersecting segments')


def _validate_is_clockwise(orig_lon, lon):
    """
    Verifies that the polygon is contructed from points in a clockwise
    direction.

    Note: The SphericalPolygon fixes a polygon with points not in a
    clockwise direction. We only need to compare a point in our
    polygon/multipolygon with the corresponding one in the SphericalPolygon
    object. We cannot use an endpoint since they are the same for a closed
    polygon.

    An AssertionError is thrown if the points are not in a clockwise
    direction
    """
    if not np.isclose(lon[1], orig_lon[1]):
        if not np.isclose(lon[1] - 360, orig_lon[1]):
            rlon = lon[::-1]
            if np.isclose(rlon[1], orig_lon[1]) or \
                    np.isclose(rlon[1] - 360, orig_lon[1]):
                raise AssertionError(
                    'invalid polygon: vertices not in clockwise direction')
            else:
                raise AssertionError(
                    'software error: compared wrong values')


def _validate_self_intersection_and_direction(ras, decs):
    """
    Verifies that the polygon does not contain self-intersecting segments
    and that the points are clockwise.

    An AssertionError is thrown if the polygon contains self-intersecting
    segments.
    """
    # use SphericalPolygon to validate our polygon
    spolygon = polygon.SphericalPolygon.from_radec(ras, decs)
    _validate_self_intersection(spolygon)
    lon, lat = six.next(spolygon.to_lonlat())
    _validate_is_clockwise(ras, lon)
