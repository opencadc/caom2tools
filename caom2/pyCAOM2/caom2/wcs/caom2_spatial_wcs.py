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

"""defines the SpatialWCS class

"""


from caom2_coord_axis2d import CoordAxis2D
from caom2.util import caom2_util as util
from caom2.caom2_object import Caom2Object

class SpatialWCS(Caom2Object):
    """this object contains the WCS information needed to convert an
    astronomical spatial location (ie. RA/DEC) into a pixel location
    in the image.  

    During ingestion a variety of extra information is  created.

    """

    def __init__(self, 
                 axis,
                 coordsys=None,
                 equinox=None,
                 resolution=None):


        self.axis = axis
        self.coordsys = coordsys
        self.equinox = equinox
        self.resolution = resolution

    @property
    def axis(self):
        """A CoordAxis2D object that contains the actual WCS values (crpix etc.)

        type: CoordAxis2D

        """
        return self._axis

    @axis.setter
    def axis(self,value):
        util.typeCheck(value, CoordAxis2D, 'axis', override=False)
        self._axis = value

    @property
    def coordsys(self):
        """The Coordinate system of the transformation, likely ICRS or FK5.

        eg.  SpatialWCS.coordsys="ICRS"
        
        type: str

        """
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        util.typeCheck(value, str, 'coordsys')
        self._coordsys = value

    @property
    def equinox(self):
        """The Equinox of the coordinate system.  

        You might think J2000, but must be expressed as a float, so in years

        unit: years
        type: float

        """
        return self._equinox

    @equinox.setter
    def equinox(self, value):
        util.typeCheck(value, float, 'equinox')
        util.valueCheck(value, 1800, 2500, 'equinox')
        self._equinox = value

    @property
    def resolution(self):
        """The spatial resolution of the image data (account for seeing/beem).

        unit: arcsec
        type: float

        """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        util.typeCheck(value, float, 'resolution')
        util.valueCheck(value, 0, 360*3600.0, 'resolution')
        self._resolution = value
