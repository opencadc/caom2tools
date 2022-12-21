# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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

"""Defines caom2.Chunk class.

"""


from builtins import str

from caom2.caom_util import int_32
from . import caom_util
from . import wcs
from .common import AbstractCaomEntity
from .common import CaomObject, OrderedEnum


class ProductType(OrderedEnum):
    """
    SCIENCE: "science"
    CALIBRATION: "calibration"
    PREVIEW: "preview"
    INFO: "info"
    NOISE: "noise"
    WEIGHT: "weight"
    AUXILIARY: "auxiliary"
    THUMBNAIL: "thumbnail"
    BIAS: "bias"
    DARK: "dark"
    FLAT: "flat"
    WAVECAL: "wavecal"
    """
    SCIENCE = "science"
    CALIBRATION = "calibration"
    PREVIEW = "preview"
    INFO = "info"
    NOISE = "noise"
    WEIGHT = "weight"
    AUXILIARY = "auxiliary"
    THUMBNAIL = "thumbnail"
    BIAS = "bias"
    DARK = "dark"
    FLAT = "flat"
    WAVECAL = "wavecal"


__all__ = ['ProductType', 'Chunk', 'ObservableAxis', 'SpatialWCS',
           'SpectralWCS', 'TemporalWCS', 'PolarizationWCS', 'CustomWCS']


class Chunk(AbstractCaomEntity):
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
                 custom_axis=None,
                 custom=None,
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
        self.custom_axis = custom_axis
        self.custom = custom
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

        eg.  Chunk.product_type = ProductType.SCIENCE

        Allowed values:
        """ + str(list(ProductType)) + """

        """

        return self._product_type

    @product_type.setter
    def product_type(self, value):
        if isinstance(value, str) and value in ProductType.names():
            # be helpful
            value = ProductType('value')
        caom_util.type_check(value, ProductType, 'product_type')
        self._product_type = value

    @property
    def naxis(self):
        """There number of dimensions in this chunk.

        type: int_32
        eg: 2

        """
        return self._naxis

    @naxis.setter
    def naxis(self, value):
        caom_util.type_check(value, int_32, 'naxis')
        caom_util.value_check(value, 0, 6, 'naxis')
        self._naxis = int_32(value) if value is not None else None

    @property
    def position_axis_1(self):
        """The first spatial axis (nominally NAXIS1).

        This is the spatial axis whose WCS is connected to CRPIX1, CD1_1, CD2_1

        eg: position_axis_1 = 1
        type: int_32

        """
        return self._position_axis_1

    @position_axis_1.setter
    def position_axis_1(self, value):
        caom_util.type_check(value, int_32, 'position_axis_1')
        #         util.valueCheck(value, 0, self.naxis, 'position_axis_1')
        self._position_axis_1 = int_32(value) if value is not None else None

    @property
    def position_axis_2(self):
        """The second spatial axis (nominally NAXIS2).

        This is the spatial axis whose WCS is connected to CRPIX2,
        CD2_2, CD1_2

        eg: position_axis_2 = 2
        type: int_32

        """
        return self._position_axis_2

    @position_axis_2.setter
    def position_axis_2(self, value):
        caom_util.type_check(value, int_32, 'position_axis_2')
        #         util.valueCheck(value, 0, self.naxis, 'position_axis_2')
        self._position_axis_2 = int_32(value) if value is not None else None

    @property
    def energy_axis(self):
        """The axis in the file that is in the energy direction.

        This should be None if the data does not contain an
        energy axis.  In this case the energy WCS maps to a
        single pixel.

        eg: energy_axis = 3
        type: int_32

        """
        return self._energy_axis

    @energy_axis.setter
    def energy_axis(self, value):
        caom_util.type_check(value, int_32, 'energy_axis')
        #         util.valueCheck(value, 0, self.naxis, 'energy_axis')
        self._energy_axis = int_32(value) if value is not None else None

    @property
    def time_axis(self):
        """The axis in the data chunk that is in the time direction.

        Can and should be None if no time sampling axis exist.

        eg. time_axis = None
        type: int_32

        """
        return self._time_axis

    @time_axis.setter
    def time_axis(self, value):
        caom_util.type_check(value, int_32, 'time_axis')
        self._time_axis = int_32(value) if value is not None else None

    @property
    def custom_axis(self):
        """The axis in the data chunk that is in a custom direction.
        """
        return self._custom_axis

    @custom_axis.setter
    def custom_axis(self, value):
        caom_util.type_check(value, int_32, 'custom_axis')
        self._custom_axis = int_32(value) if value is not None else None

    @property
    def polarization_axis(self):
        """The axis in the data chunk that is in the polarization direction.

        Likely None...

        eg. polarization_axis = None
        type: int_32

        """
        return self._polarization_axis

    @polarization_axis.setter
    def polarization_axis(self, value):
        caom_util.type_check(value, int_32, 'polarization_axis')
        caom_util.value_check(value, 0, 2 ** 32, 'polariztion_axis')
        self._polarization_axis = int_32(value) if value is not None else None

    @property
    def observable_axis(self):
        """Used when on of the dimensions of the file contains?? ?

        type: int_32

        """
        return self._observable_axis

    @observable_axis.setter
    def observable_axis(self, value):
        caom_util.type_check(value, int_32, 'obserable_axis')
        #         util.valueCheck(value, 0, 1E10, 'observable_axis')
        self._observable_axis = int_32(value) if value is not None else None

    @property
    def observable(self):
        """An obserable that is contained in the chunk.

        Observables are quantities that are recorded directly??

        """
        return self._observable

    @observable.setter
    def observable(self, value):
        caom_util.type_check(value, ObservableAxis, 'observable_axis')
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
        caom_util.type_check(value, SpatialWCS, 'position')
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
        caom_util.type_check(value, SpectralWCS, 'energy')
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
        caom_util.type_check(value, TemporalWCS, 'time')
        self._time = value

    @property
    def polarization(self):
        """The PolarizationWCS of the observation.

        usually None

        type: PolarizationWCS

        """
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        caom_util.type_check(value, PolarizationWCS, 'polarization')
        self._polarization = value

    @property
    def custom(self):
        """The CustomWCS of the observation.

        usually None

        type: CustomWCS

        """
        return self._custom

    @custom.setter
    def custom(self, value):
        caom_util.type_check(value, CustomWCS, 'custom')
        self._custom = value


class ObservableAxis(CaomObject):
    """The slice through the data structure that provides the thing being
    described by this Axis.

    this data structure is used when a file contains data that has set
    of measured values (the dependent variable) and might also have
    the coordinate at which those values are measured in another Slice
    of the file.

    The Slice refers to a column in a FITS image. The bin is the
    column index and ctype/cunit for the axis describe what the column
    contains.

    eg.

    NAXIS=2
    NAXIS1=3
    NAXIS2=N
    l1 f1 s1
    l2 f2 s2
    l3 f3 s3
    .
    .
    lN fN sN

    where l? is the wavelength at which a measure of flux f? has been
    made and l is the first column of a FITS data structure that is
    3,N in size.  s? is a third slice that would be used to define
    another observable. When defining the s? observable the independent
    variable must be defined for that ObservableAxis too.

    The l/f observable would be recorded as

    dependent=Slice(Axis('wave','nm'),bin=1)
    independent=Slice(Axis('flux','Jy'),bin=2)
    Chunk.observable_axis=ObservableAxis(dependent, independent)


    """

    def __init__(self, dependent, independent=None):
        self.dependent = dependent
        self.independent = independent

    @property
    def dependent(self):
        """The dependent (y) variable slice.

        A slice provides the bin and the type/unit of the observable axis
        """

        return self._dependent

    @dependent.setter
    def dependent(self, value):
        caom_util.type_check(value, wcs.Slice, 'dependent', override=False)
        self._dependent = value

    @property
    def independent(self):
        """The dependent (y) variable slice.

        A slice that provides the pixel value and the type/unit for
        conversion

        """
        return self._independent

    @independent.setter
    def independent(self, value):
        caom_util.type_check(value, wcs.Slice, "independent")
        self._independent = value


class SpatialWCS(CaomObject):
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
        """A CoordAxis2D object that contains
        the actual WCS values (crpix etc.)

        type: CoordAxis2D

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        caom_util.type_check(value, wcs.CoordAxis2D, 'axis', override=False)
        self._axis = value

    @property
    def coordsys(self):
        """The Coordinate system of the transformation, likely ICRS or FK5.

        eg.  SpatialWCS.coordsys="ICRS"

        type: unicode string

        """
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        caom_util.type_check(value, str, 'coordsys')
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
        caom_util.type_check(value, float, 'equinox')
        caom_util.value_check(value, 1800, 2500, 'equinox')
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
        caom_util.type_check(value, float, 'resolution')
        caom_util.value_check(value, 0, 360 * 3600.0, 'resolution')
        self._resolution = value


class SpectralWCS(CaomObject):
    """A transformation that maps pixel coordinates to spectral ones.

    Note that a 2D image has implicit 'spectral' (and temporal)
    dimension that is one pixel in size.  Basically pixels are really
    multidimensional voxels.

    Due to FITS standards this pixel starts at 0.5 and runs to 1.5
    with the centre being at 1.0

    """

    def __init__(self,
                 axis,
                 specsys,
                 ssysobs=None,
                 ssyssrc=None,
                 restfrq=None,
                 restwav=None,
                 velosys=None,
                 zsource=None,
                 velang=None,
                 bandpass_name=None,
                 transition=None,
                 resolving_power=None
                 ):
        """The SpectralWCS can be defined in a number of different ways, to
        define one you must provide a CoordAxis1D object that maps the
        pixels to WCS values and a reference specsys.  After that the
        user can add what ever parts seam useful.  More info is more
        helpful for searching.

        """

        self.axis = axis
        self.specsys = specsys
        self.ssysobs = ssysobs
        self.ssyssrc = ssyssrc
        self.restfrq = restfrq
        self.restwav = restwav
        self.velosys = velosys
        self.zsource = zsource
        self.velang = velang
        self.bandpass_name = bandpass_name
        self.transition = transition
        self.resolving_power = resolving_power

    @property
    def axis(self):
        """A 1D coordinate axis object that contains the pix/wcs
        transformation values.

        eg.  CoordAxis1D(Axis('wave','flux'),...)

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        caom_util.type_check(value, wcs.CoordAxis1D, 'axis', override=False)
        self._axis = value

    @property
    def specsys(self):
        """describes the reference frame in use for the spectral-axis
        coordinate(s).

        eg. BARYCENT

        type: unicode string
        """
        return self._specsys

    @specsys.setter
    def specsys(self, value):
        caom_util.type_check(value, str, 'specsys', override=False)
        self._specsys = value

    @property
    def ssysobs(self):
        """describes the spectral reference frame that is constant over the
        range of the non-spectral world coordinates

        For example, for a large image the the wavelength at the edges
        is different from the centres.  This reference frame is one where they
        are not different.

        Nominally 'TOPOCENT'

        type: unicode string
        """
        return self._ssysobs

    @ssysobs.setter
    def ssysobs(self, value):
        caom_util.type_check(value, str, 'ssysobs')
        self._ssysobs = value

    @property
    def ssyssrc(self):
        """The reference frame in which zsource is expressed.

        eg. BARYCENT
        type: string
        """
        return self._ssyssrc

    @ssyssrc.setter
    def ssyssrc(self, value):
        caom_util.type_check(value, str, 'ssyssrc')
        self._ssyssrc = value

    @property
    def restfrq(self):
        """The frequency of the spectal feature being observed.

        unit: Hz
        type: float
        """
        return self._restfrq

    @restfrq.setter
    def restfrq(self, value):
        caom_util.type_check(value, float, 'restfrq')
        self._restfrq = value

    @property
    def restwav(self):
        """The wavelength of spectral feature being observed,
        not the wavelength observed but the wavelength of the
        feature when at rest..

        unit: m
        type: float
        """
        return self._restwav

    @restwav.setter
    def restwav(self, value):
        caom_util.type_check(value, float, 'restwav')
        self._restwav = value

    @property
    def velosys(self):
        """Relative radial velocity between the observer and the selected
        standard of rest in the direction of the celestial reference
        coordinate.

        eg. 26000 m/s


        unit: m/s
        type: float
        """
        return self._velosys

    @velosys.setter
    def velosys(self, value):
        caom_util.type_check(value, float, 'velosys')
        self._velosys = value

    @property
    def zsource(self):
        """The redshift of the source emitting the photons.

        almost always None

        unit: z
        type: float
        """
        return self._zsource

    @zsource.setter
    def zsource(self, value):
        caom_util.type_check(value, float, 'zsource')
        caom_util.value_check(value, -0.5, 1200, 'zsource')
        self._zsource = value

    @property
    def velang(self):
        """I don't know what this is... angle of the velocity ??? """
        return self._velang

    @velang.setter
    def velang(self, value):
        caom_util.type_check(value, float, 'velang')
        self._velang = value

    @property
    def bandpass_name(self):
        """string the represent the bandpass of the observation.

        eg. r'
        type: unicode string
        """
        return self._bandpass_name

    @bandpass_name.setter
    def bandpass_name(self, value):
        caom_util.type_check(value, str, 'bandpass_name')
        self._bandpass_name = value

    @property
    def transition(self):
        """which molecular transition has been observed.

        type: EnergyTransition object (see caom2.EnergyTransition for help)
        """
        return self._transition

    @transition.setter
    def transition(self, value):
        caom_util.type_check(value, wcs.EnergyTransition, "transition")
        self._transition = value

    @property
    def resolving_power(self):
        """The R value of the spectal coverage.

        Normally this is something like dlamda/lamda

        unit:  RATIO
        type: float
        """
        return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        caom_util.type_check(value, float, 'resolving_power')
        caom_util.value_check(value, 0, 1E8, 'resolving_power')
        self._resolving_power = value


class TemporalWCS(CaomObject):
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
        caom_util.type_check(value, wcs.CoordAxis1D, 'axis', override=False)
        self._axis = value

    @property
    def timesys(self):
        """The time scale that you are using, almost alwasy UTC.

        eg.  timesys = "UTC"
        type: unicode string
        """
        return self._timesys

    @timesys.setter
    def timesys(self, value):
        caom_util.type_check(value, str, 'timesys')
        self._timesys = value

    @property
    def trefpos(self):
        """ specifies the spatial location at which the time is valid, either
        where the observation was made or the point in space for which
        light-time corrections have been applied.

        eg. trefpos = "TOPOCENTER"
        type: unicode string
        """
        return self._trefpos

    @trefpos.setter
    def trefpos(self, value):
        caom_util.type_check(value, str, 'trefpos')
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
        caom_util.type_check(value, float, 'mjdref')
        # set the limits to be after 1800 but before year 2050
        caom_util.value_check(value, -22000, 70000, 'mjdref')
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
        caom_util.type_check(value, float, "exposure")
        caom_util.value_check(value, 0, float("inf"), "exposure")
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
        caom_util.type_check(value, float, 'resolution')
        caom_util.value_check(value, 0, 100 * 365 * 24 * 3600.0, "resolution")
        self._resolution = value


class PolarizationWCS(CaomObject):
    """A WCS structure that describes the relation ship between a pixel
    location and the polarization value.

    """

    def __init__(self, axis):
        """Set up a CoordAxis1D object to represent the Polarization.

        """

        self.axis = axis

    @property
    def axis(self):
        """A CoordAxis1D object that describes the pixel/value relation ship
        for polarization of the data.

        type: CoordAxis1D

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        caom_util.type_check(value, wcs.CoordAxis1D, 'axis', override=False)
        if value.axis.ctype != 'STOKES':
            raise ValueError('CTYPE must be STOKES')
        self._axis = value


class CustomWCS(CaomObject):
    """A WCS structure that describes the relation ship between a pixel
    location and a custom value.

    """

    def __init__(self, axis):
        """Set up a CoordAxis1D object to represent the Custom Axis.

        """
        caom_util.type_check(axis, wcs.CoordAxis1D, 'axis', override=False)
        self._axis = axis

    @property
    def axis(self):
        """A CoordAxis1D object that describes the pixel/value relation ship
        for customized axis of the data.

        type: CoordAxis1D

        """
        return self._axis
