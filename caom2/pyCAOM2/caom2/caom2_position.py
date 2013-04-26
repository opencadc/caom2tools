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

"""defines the caom2.Position class"""

from types.caom2_box import Box
from types.caom2_circle import Circle
from types.caom2_location import Location
from types.caom2_polygon import Polygon
from wcs.caom2_dimension2d import Dimension2D
from caom2_object import Caom2Object
import caom2.util as util

class Position(Caom2Object):
    """ Position """

    def __init__(self, location=None, 
                 bounds=None, 
                 dimension=None,
                 resolution=None,
                 sample_size=None,
                 time_dependent=None
                 ):
        """
        Initialize a Time instance.

        Arguments:
        None
        """
        self._location = location
        self._bounds = bounds
        self._dimension = dimension
        self._resolution = resolution
        self._sample_size = sample_size
        self._time_dependent = time_dependent

    # Properties

    @property
    def location(self):
        """ Location """
        return self._location

    @location.setter
    def location(self, value):
        if value is not None:
            assert isinstance(value, Location), (
                    "location is not a Location: {0}".format(value))
        self._location = value

    @property
    def bounds(self):
        """ Bounds """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            assert isinstance(value, (Box, Circle, Location, Polygon)), (
                    "bounds is not a Shape: {0}".format(value))
        self._bounds = value

    @property
    def dimension(self):
        """ Dimension """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            assert isinstance(value, Dimension2D), (
                    "dimension is not a Dimension2D: {0}".format(value))
        self._dimension = value

    @property
    def resolution(self):
        """ Resolution """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        if value is not None:
            assert isinstance(value, float), (
                    "resolution is not a float: {0}".format(value))
        self._resolution = value

    @property
    def sample_size(self):
        """ Sample size """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        if value is not None:
            assert isinstance(value, float), (
                    "sample size is not a float: {0}".format(value))
        self._sample_size = value

    @property
    def time_dependent(self):
        """ Time dependent """
        return self._time_dependent

    @time_dependent.setter
    def time_dependent(self, value):
        if value is not None:
            assert isinstance(value, bool), (
                    "time dependent is not a bool: {0}".format(value))
        self._time_dependent = value
