# -*- coding: latin-1 -*-
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

""" defines the CoordAxis2D class"""

from caom2_axis import Axis
from caom2_coord_error import CoordError
from caom2_coord_range2d import CoordRange2D
from caom2_coord_circle2d import CoordCircle2D
from caom2_coord_polygon2d import CoordPolygon2D
from caom2_coord_function2d import CoordFunction2D
from caom2.caom2_object import Caom2Object
from caom2.util import caom2_util as util


class CoordAxis2D(object):
    """This object hold the metadata need to transform a 2D pixel 
    array (say an image) into a World position, say RA/DEC

    """

    def __init__(self, axis1, axis2,
                 error1=None, error2=None,
                 range=None, bounds=None,
                 function=None):

        self.axis1 = axis1
        self.axis2 = axis2
        self.error1 = error1
        self.error2 = error2
        self.range = range
        self.bounds = bounds
        self.function = function


    @property
    def axis1(self):
        """An axis object that desciribes the first dimension of this 2d system.

        eg. axis1=Axis("RA","deg")

        """
        return self._axis1

    @axis1.setter
    def axis1(self, value):
        util.typeCheck(value, Axis, "axis1", override=False)
        self._axis1 = value



    @property
    def axis2(self):
        """An axis objet that describes the 2nd dimensiotn of this 2d coord
        system.

        eg. axis2=Axis("DEG","deg")

        """
        return self._axis2

    @axis2.setter
    def axis2(self,value):
        util.typeCheck(value, Axis, "axis2", override=False)
        self._axis2 = value

    @property
    def error1(self):
        """An object that descibes the uncertainty in the pix/world transform.

        eg. CoordError()
        type: CoordError

        """
        return self._error1

    @error1.setter
    def error1(self, value):
        util.typeCheck(value, CoordError, 'error1')
        self._error1 = value

    @property
    def error2(self):
        """An object that describes the uncertainty in the pix/world transform
        for the 2nd axis

        type: CoordError

        """
        return self._error2

    @error2.setter
    def error2(self, value):
        util.typeCheck(value, CoordError, "error2")
        self._error2 = value

    @property
    def range(self):
        """Coordinate range defined by this CoordAxis2d object.
        
        type: CoordRange2D
        """
        return self._range

    @range.setter
    def range(self, value):
        util.typeCheck(value, CoordRange2D, 'range')
        self._range = value

    @property
    def bounds(self):
        """The Coordinate boundary mapped by this CoordAxis2D object.

        ag. CoordPolygon2d((234,10),(234,11),(233,11),(234,11),(234,10))
        type:  CoordPolygon2D or CoordCircle2D
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        util.typeCheck(value, (CoordCircle2D, CoordPolygon2D), 'bounds')
        self._bounds = value

    @property
    def function(self):
        """A function object that describes the relation between pixels and wcs.

        ag.  CoordFunction2D (see the help for that puppy)
        type: CoordFunction2D
        """
        return self._function

    @function.setter
    def function(self, value):
        util.typeCheck(value, CoordFunction2D, 'function')
        self._function = value
