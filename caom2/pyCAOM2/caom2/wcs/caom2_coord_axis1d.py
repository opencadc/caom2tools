#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2010.                            (c) 2010.
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
#***********************************************************************
#

"""Definition of the CoordAxis1D class"""

from caom2_axis import Axis
from caom2_coord_error import CoordError
from caom2_coord_range1d import CoordRange1D
from caom2_coord_bounds1d import CoordBounds1D
from caom2_coord_function1d import CoordFunction1D
from caom2.util import caom2_util as util
from caom2.caom2_object import Caom2Object


class CoordAxis1D(Caom2Object):
    """Holds the metadata needed to transform a 1D pixel value into a
    World Coordinate value.

    """

    def __init__(self, axis, error=None,
                 _range=None, bounds=None, function=None):

        self.axis = axis
        self.error = error
        self.range = _range
        self.bounds = bounds
        self.function = function

    @property
    def axis(self):
        """An axis object which describes the type and units of this
        coordinate axis.

        eg. Axis(ctype1,cunit1)

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        util.typeCheck(value, Axis, "Axis", override=False)
        self._axis = value

    @property
    def error(self):
        """A CoordError object that describes the uncertainty in the.

        eg.  CoordError(syser=0.1, rnder=0.1)
        unit: cunit1 [of axis]

        """
        return self._error

    @error.setter
    def error(self, value):
        util.typeCheck(value, CoordError, 'error')
        self._error = value

    @property
    def range(self):
        """A range that defines a coordinate transformation.

        the transform is a linear interpolation over the range
        given which is a specified as a set of two pix/val reference pair.
        eg.  CoordRange1D(start=RefCoord(pix1,val1),end=RefCoord(pix2,val2))
        unit: same as the axis you are defining.

        """
        return self._range

    @range.setter
    def range(self, value):
        util.typeCheck(value, CoordRange1D, 'range')
        self._range = value

    @property
    def bounds(self):
        """A polygon that defines the boundary of this axis, in 1D.

        eg. CoordBounds1D(ListOfRanges())
        The ranges are like those given for the range attribute.

        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        util.typeCheck(value, CoordBounds1D, "bounds")
        self._bounds = value

    @property
    def function(self):
        """A linear function that describes the tranformation between pixel
        and world coordinate value.

        Since this is a 1D object and linear, the function is
        y = m*x + b.
        eg. CoordFunction1D(naxis, delta, RefCoord)

        """
        return self._function

    @function.setter
    def function(self, value):
        util.typeCheck(value, CoordFunction1D, 'function')
        self._function = value
