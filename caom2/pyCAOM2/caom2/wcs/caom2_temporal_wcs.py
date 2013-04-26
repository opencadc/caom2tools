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

"""The definition of the TemporalWCS class"""


from caom2_coord_axis1d import CoordAxis1D
from caom2.caom2_object import Caom2Object
from caom2.util import caom2_util as util

class TemporalWCS(Caom2Object):
    """Describes the Time variation within the data.

    In the case of a single exposure, define the center of the first
    pixel (0.5) as the MJD of the exposure and size of the pixel 
    as the exposure time."""

    def __init__(self, 
                 axis,
                 timesys=None,
                 trefpos=None,
                 mjdref=None,
                 exposure=None,
                 resolution=None
                 ):

        self.axis = axis
        self.timesys = timesys
        self.trefpos = trefpos
        self.mjdref = mjdref
        self.exposure = exposure
        self.resolution = resolution
        
    @property
    def axis(self):
        """A CoordAxis1D object that describes the TemporalWCS transform.

        eg. CoordAxis1D(Axis) """
        return self._axis

    @axis.setter
    def axis(self, value):
        util.typeCheck(value, CoordAxis1D, 'axis', override=False)
        self._axis = value

    @property
    def timesys(self):
        """The time scale that you are using, almost alwasy UTC.
        
        eg.  timesys = "UTC" 
        type: str
        """
        return self._timesys

    @timesys.setter
    def timesys(self, value):
        util.typeCheck(value, str, 'timesys')
        self._timesys = value

    @property
    def trefpos(self):
        """ specifies the spatial location at which the time is valid, either
        where the observation was made or the point in space for which
        light-time corrections have been applied.

        eg. trefpos = "TOPOCENTER"
        type: str
        """
        return self._trefpos

    @trefpos.setter
    def trefpos(self, value):
        util.typeCheck(value, str, 'trefpos')
        self._trefpos = value

    @property
    def mjdref(self):
        """The Modified Julian Date of the at the reference location of the
        location of the TimeWCS (aka. pixel 0.5).  Nominally this the start
        of the exposure.  

        Why 0.5? FITS: the middle of the first pixel
        is defined as pixel value 1.0  so, the start of that pixel
        is location 0.5 

        eg. mjdref = 567643.1234
        unit: d
        type: float
        """
        return self._mjdref

    @mjdref.setter
    def mjdref(self, value):
        util.typeCheck(value, float, 'mjdref')
        ### set the limits to be after 1800 but before year 2050
        util.valueCheck(value, -22000, 70000, 'mjdref')
        self._mjdref = value

    @property
    def exposure(self):
        """The median exposure time per pixel.

        The exposure time if this not a time cube you are describing.

        eg. exposure = 100.0
        unit: s
        type: float
        """
        return self._exposure

    @exposure.setter
    def exposure(self, value):
        util.typeCheck(value, float, "exposure")
        util.valueCheck(value, 0, 30*24*3600.0, "exposure")
        self._exposure = value

    @property
    def resolution(self):
        """the resolution of the time sampling available.

        Normally this is going to be the same as the exposure above,
        but a stack of exposures taken over many months has a very 
        large value for resolution while the exposure value is just
        the sum of the individual exposure times. 


        eg. resolution = 100.0s
        unit: s
        type: float
        """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        util.typeCheck(value, float, 'resolution')
        util.valueCheck(value, 0, 100*365*24*3600.0, "resolution")
        self._resolution = value
