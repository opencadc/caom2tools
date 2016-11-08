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

"""defines a 'CoordFunction2D' class

"""

from caom2_dimension2d import Dimension2D
from caom2_coord2d import Coord2D
from caom2.caom2_object import Caom2Object
from caom2.util import caom2_util as util


class CoordFunction2D(Caom2Object):
    """Describes the parameters needed for the standard CD matrix.

    defines a linear translation between pixel and WCS.

    """

    def __init__(self, dimension, ref_coord, cd11, cd12, cd21, cd22):
        self.dimension = dimension
        self.ref_coord = ref_coord
        self.cd11 = cd11
        self.cd12 = cd12
        self.cd21 = cd21
        self.cd22 = cd22

    @property
    def dimension(self):
        """A Dimension2D object that holds the lengths of the axis

        eg.  Diemnsion2D(naxis1=1024,naxis2=2048)
        type: Dimension2D

        """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        util.typeCheck(value, Dimension2D, 'dimension', override=False)
        self._dimension = value

    @property
    def ref_coord(self):
        """A Coord2D object that holds the reference pixel location

        eg. Coord2D((crpix1,crval1),(crpix2,crval2))
        type: Coord2D

        """
        return self._ref_coord

    @ref_coord.setter
    def ref_coord(self, value):
        util.typeCheck(value, Coord2D, 'ref_coord', override=False)
        self._ref_coord = value

    @property
    def cd11(self):
        """The CD1_1 value (depenence of RA scale on x-pixel value)

        eg. cd11 = 5E-5
        unit: deg/pix
        type: float

        """
        return self._cd11

    @cd11.setter
    def cd11(self, value):
        util.typeCheck(value, float, 'cd11', override=False)
        self._cd11 = value

    @property
    def cd12(self):
        """The CD1_2 value (depenence of RA scale on y-pixel value)

        eg. cd12 = 5E-10
        unit: deg/pix
        type: float

        """
        return self._cd12

    @cd12.setter
    def cd12(self, value):
        util.typeCheck(value, float, 'cd12', override=False)
        self._cd12 = value

    @property
    def cd21(self):
        """The CD1_1 value (depenence of DEC scale on x-pixel value)

        eg. cd11 = 5E-10
        unit: deg/pix
        type: float

        """
        return self._cd21

    @cd21.setter
    def cd21(self, value):
        util.typeCheck(value, float, 'cd21', override=False)
        self._cd21 = value

    @property
    def cd22(self):
        """The CD2_2 value (depenence of DEC scale on y-pixel value)

        eg. cd12 = 5E-5
        unit: deg/pix
        type: float

        """
        return self._cd22

    @cd22.setter
    def cd22(self, value):
        util.typeCheck(value, float, 'cd22', override=False)
        self._cd22 = value
