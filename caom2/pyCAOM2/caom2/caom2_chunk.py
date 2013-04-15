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

"""Defines caom2.Chunk class.

"""

from caom2_entity import AbstractCaom2Entity
from caom2_enums import ProductType
from wcs.caom2_observable_axis import ObservableAxis
from wcs.caom2_spatial_wcs import SpatialWCS
from wcs.caom2_spectral_wcs import SpectralWCS
from wcs.caom2_polarization_wcs import PolarizationWCS
from wcs.caom2_temporal_wcs import TemporalWCS
import util.caom2_util as util

class Chunk(AbstractCaom2Entity):
    """A caom2.Chunk object.  A chunk is a peice of file part.

    eg.  a column in a Table extension of a FITS file.

    The chunk is characterised by world coordinate system (WCS)
    metadata plus an extra axis to describe different observables (the
    measured values) stored within the data. Different chunks can be
    defined which vary only on the range of coordinate values they
    include. For example, if a single data array contains different
    observable quantities then one can define a chunk (perhaps
    representing different slices through a stored array) with each
    slice having a different product type. 

    Chunks can also be used to define arbitrary tiles in a large data
    array; this is useful if there is no WCS solution to describe the
    mapping of sky to pixel coordinates but one still wants to be able
    to extract smaller sections of the data (e.g. one chunk).

    """

    def __init__(self, product_type=None,
                       naxis=None,
		       position_axis_1=None,
		       position_axis_2=None,
		       position=None,
		       energy_axis=None,
		       energy=None,
		       time_axis=None,
		       time=None,
		       polarization_axis=None,
		       polarization=None,
		       observable_axis=None,
		       observable=None,
		       ):

        super(Chunk, self).__init__()
        self.product_type = product_type
        self.naxis = naxis
        self.position_axis_1 = position_axis_1
        self.position_axis_2 = position_axis_2
        self.energy_axis = energy_axis
        self.time_axis = time_axis
        self.polarization_axis = polarization_axis
        self.observable_axis = observable_axis
        self.observable = observable
        self.position = position
        self.energy = energy
        self.time = time
        self.polarization = polarization

    @property
    def product_type(self):
        """A word that describes the content of the chunk.

        eg.  Chunk.product_type = ProductType('SCIENCE')

        Allowed values:
        """+str(ProductType.names())+"""

        """

        return self._product_type

    @product_type.setter
    def product_type(self, value):
        if isinstance(value, str) and value in ProductType.names():
            ## be helpful
            value = ProductType('value')
        util.typeCheck(value, ProductType, 'product_type')
        self._product_type = value

    @property
    def naxis(self):
        """There number of dimensions in this chunk.

        type: int
        eg: 2

        """
        return self._naxis

    @naxis.setter
    def naxis(self, value):
        util.typeCheck(value, int, 'naxis')
        util.valueCheck(value, 0, 5, 'naxis')
        self._naxis = value

    @property
    def position_axis_1(self):
        """The first spatial axis (nominally NAXIS1).

        This is the spatial axis whose WCS is connected to CRPIX1, CD1_1, CD2_1
        
        eg: position_axis_1 = 1
        type: int

        """
        return self._position_axis_1

    @position_axis_1.setter
    def position_axis_1(self, value):
        util.typeCheck(value, int, 'position_axis_1')
        util.valueCheck(value, 0, self.naxis, 'position_axis_1')
        self._position_axis_1 = value

    @property
    def position_axis_2(self):
        """The second spatial axis (nominally NAXIS2).

        This is the spatial axis whose WCS is connected to CRPIX2,
        CD2_2, CD1_2 

        eg: position_axis_2 = 2 
        type: int

        """
        return self._position_axis_2

    @position_axis_2.setter
    def position_axis_2(self, value):
        util.typeCheck(value, int, 'position_axis_2')
        util.valueCheck(value, 0, self.naxis, 'position_axis_2')
        self._position_axis_2 = value

    @property
    def energy_axis(self):
        """The axis in the file that is in the energy direction.

        This should be None if the data does not contain an
        energy axis.  In this case the energy WCS maps to a
        single pixel.

        eg: energy_axis = 3
        type: int

        """
        return self._energy_axis

    @energy_axis.setter
    def energy_axis(self, value):
        util.typeCheck(value, int, 'energy_axis')
        util.valueCheck(value, 0, self.naxis, 'energy_axis')
        self._energy_axis = value

    @property
    def time_axis(self):
        """The axis in the data chunk that is in the time direction.

        Can and should be None if no time sampling axis exist.
        
        eg. time_axis = None
        type: int

        """
        return self._time_axis

    @time_axis.setter
    def time_axis(self, value):
        util.typeCheck(value, int, 'polarization_axis')
        util.valueCheck(value, 0, self._naxis, 'polarization_axis')
        self._time_axis = value

    @property
    def polarization_axis(self):
        """The axis in the data chunk that is in the polarization direction.

        Likely None... 

        eg. polarization_axis = None
        type: int
        
        """
        return self._polarization_axis

    @polarization_axis.setter
    def polarization_axis(self, value):
        util.typeCheck(value, int, 'polarization_axis')
        util.valueCheck(value, 0, self._naxis, 'polariztion_axis')
        self._polarization_axis = value

    @property
    def observable_axis(self):
        """Used when on of the dimensions of the file contains?? ?

        type: int

        """
        return self._observable_axis

    @observable_axis.setter
    def observable_axis(self, value):
        util.typeCheck(value, int, 'obserable_axis')
        util.valueCheck(value, 0, 1E10, 'observable_axis')
        self._observable_axis = value

    @property
    def observable(self):
        """An obserable that is contained in the chunk. 
        
        Observables are quantities that are recorded directly??

        """
        return self._observable

    @observable.setter
    def observable(self, value):
        util.typeCheck(value,ObservableAxis,'observable_axis')
        self._observable = value

    @property
    def position(self):
        """A SpatialWCS object associated with this chunk.

        The spatialWCS describes the relation between the position_axis
        values and the world coordinate.

        type: SpatialWCS.
        
        """
        return self._position

    @position.setter
    def position(self, value):
        util.typeCheck(value,SpatialWCS,'position')
        self._position = value

    @property
    def energy(self):
        """A SpectralWCS object associated with this chunk.

        Even if energy_axis is None an SpectralWCS should still
        be defined.  The SpectalWCS in this case will be one pixel
        in dimension.

        type: SpectralWCS
        
        """
        return self._energy

    @energy.setter
    def energy(self, value):
        util.typeCheck(value, SpectralWCS, 'energy')
        self._energy = value

    @property
    def time(self):
        """The TemporalWCS object associated with this chunk.

        Even if time_axis is None you should define the TimeWCS
        to convey when you observation was taken.
        
        type: TemporalWCS
        
        """
        return self._time

    @time.setter
    def time(self, value):
        util.typeCheck(value, TemporalWCS, 'time')
        self._time = value

    @property
    def polarization(self):
        """The PolarizationWCS of the observation.

        ususally None
        
        type: PolarizationWCS

        """
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        util.typeCheck(value, PolarizationWCS, 'polarization')
        self._polarization = value
