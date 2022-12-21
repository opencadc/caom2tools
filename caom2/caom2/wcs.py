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

from builtins import str, int

from . import caom_util
from . import common

__all__ = ['Axis', 'Coord2D', 'CoordAxis1D', 'CoordAxis2D', 'CoordBounds1D',
           'CoordBounds2D', 'CoordCircle2D', 'CoordError', 'CoordFunction1D',
           'CoordFunction2D', 'CoordPolygon2D', 'CoordRange1D', 'CoordRange2D',
           'Dimension2D', 'EnergyTransition', 'RefCoord', 'Slice',
           'ValueCoord2D']


class Axis(common.CaomObject):
    """the Axis class holds the definition of the axis type and units"""

    def __init__(self, ctype, cunit=None):
        self.ctype = ctype
        self.cunit = cunit

    @property
    def ctype(self):
        """The Coordinate Type value for this axis.

        eg. DEC
        type: unicode string

        """
        return self._ctype

    @ctype.setter
    def ctype(self, value):
        caom_util.type_check(value, str, 'ctype', override=False)
        self._ctype = value

    @property
    def cunit(self):
        """The unit of the coordinate that results after transform.

        eg. deg
        type: unicode string

        """
        return self._cunit

    @cunit.setter
    def cunit(self, value):
        caom_util.type_check(value, str, 'cunit')
        self._cunit = value


class Coord2D(common.CaomObject):
    """Represents the reference point.

    eg:  Coord2D(RefCoord(crpix1,crval1),RefCoord(crpix2,crval2))
    """

    def __init__(self, coord1, coord2):
        self.coord1 = coord1
        self.coord2 = coord2

    @property
    def coord1(self):
        """A coordinate axis1 coordinate/pix pair, crpix1, crval1.

        eg.  RefCoord(crpix1, crval1)
        """
        return self._coord1

    @coord1.setter
    def coord1(self, value):
        caom_util.type_check(value, RefCoord, 'coord1', override=False)
        self._coord1 = value

    @property
    def coord2(self):
        """The axis2 coordinate reference pair (ei. crpix2/crval2.

        eg. RefCorrd(crpix2, crval2)
        """
        return self._coord2

    @coord2.setter
    def coord2(self, value):
        caom_util.type_check(value, RefCoord, 'coord2', override=False)
        self._coord2 = value


class CoordAxis1D(common.CaomObject):
    """Holds the metadata needed to transform a 1D pixel value into a
    World Coordinate value.

    """

    def __init__(self, axis, error=None, range=None,
                 bounds=None, function=None):
        self.axis = axis
        self.error = error
        self.range = range
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
        caom_util.type_check(value, Axis, "Axis", override=False)
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
        caom_util.type_check(value, CoordError, 'error')
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
        caom_util.type_check(value, CoordRange1D, 'range')
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
        caom_util.type_check(value, CoordBounds1D, "bounds")
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
        caom_util.type_check(value, CoordFunction1D, 'function')
        self._function = value


class CoordAxis2D(common.CaomObject):
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
        """An axis object that desciribes the first dimension
        of this 2d system.

        eg. axis1=Axis("RA","deg")

        """
        return self._axis1

    @axis1.setter
    def axis1(self, value):
        caom_util.type_check(value, Axis, "axis1", override=False)
        self._axis1 = value

    @property
    def axis2(self):
        """An axis objet that describes the 2nd dimensiotn of this 2d coord
        system.

        eg. axis2=Axis("DEG","deg")

        """
        return self._axis2

    @axis2.setter
    def axis2(self, value):
        caom_util.type_check(value, Axis, "axis2", override=False)
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
        caom_util.type_check(value, CoordError, 'error1')
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
        caom_util.type_check(value, CoordError, "error2")
        self._error2 = value

    @property
    def range(self):
        """Coordinate range defined by this CoordAxis2d object.

        type: CoordRange2D
        """
        return self._range

    @range.setter
    def range(self, value):
        caom_util.type_check(value, CoordRange2D, 'range')
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
        caom_util.type_check(value, (CoordCircle2D, CoordPolygon2D), 'bounds')
        self._bounds = value

    @property
    def function(self):
        """A function object that describes the relation
        between pixels and wcs.

        ag.  CoordFunction2D (see the help for that puppy)
        type: CoordFunction2D
        """
        return self._function

    @function.setter
    def function(self, value):
        caom_util.type_check(value, CoordFunction2D, 'function')
        self._function = value


class CoordBounds1D(common.CaomObject):
    """Contains the bounds for a 1D axis, a list of ranges

    """

    def __init__(self, samples=None):
        if samples is None:
            samples = caom_util.TypedList(CoordRange1D, )
        self.samples = samples

    @property
    def samples(self):
        """A list of CoordRange1D objects that define the
        boundary of a 1D axis.

        see also caom2.util.TypedList and caom2.wcs.CoordRange1D

        eg.
        samples.add(CoordRange1D(RefCoord(pix,val),RefCoord(pix,val)))

        """
        return self._samples

    @samples.setter
    def samples(self, value):
        caom_util.type_check(value, caom_util.TypedList, 'samples',
                             override=False)
        self._samples = value


class CoordBounds2D(common.CaomObject):
    """Contains the bounds for a 2D axis

    """

    def __init__(self, bounds):
        if (isinstance(bounds, CoordCircle2D) or
                isinstance(bounds, CoordPolygon2D)):
            self.bounds = bounds
        else:
            raise TypeError(
                "Expected CoordCircle2D or CoordPolygon2D, received {}"
                .format(type(bounds)))

    @property
    def bounds(self):
        """The bounds expressed as a circle or polygon.

        eg CoordBounds2D(CoordCircle2D())
        type: CoordCircle2D or CoordPolygon2D
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        self._bounds = value


class CoordCircle2D(common.CaomObject):
    """A circle expressed in both pixel and WCS value coordinates.

    These objects are used to map out the bounds of spatial WCS.

    Currently radius is only given in WCS units.. should be both.

    """

    def __init__(self, center, radius):
        self.center = center
        self.radius = radius

    @property
    def center(self):
        """The pixel/world coordinate location of the centre.

        eg ValueCoord2D(coord1, coord2)
        type: ValueCoord2D
        """
        return self._center

    @center.setter
    def center(self, value):
        caom_util.type_check(value, ValueCoord2D, 'center', override=False)
        self._center = value

    @property
    def radius(self):
        """The radius of the circle.

        NOTE::: This should likely be a RefCoord too...

        unit: same as centre which is pix/cunit
        type: float
        """
        return self._radius

    @radius.setter
    def radius(self, value):
        caom_util.type_check(value, float, 'radius', override=False)
        caom_util.value_check(value, 0, 1E10, 'radius')
        self._radius = value


class CoordError(common.CaomObject):
    """Holds the systematic (syser) and random (rnder) error on a
    coordinate.

    The concept is that these values are functions of the coordinate
    transformation but likely be expressed with sufficient precision
    (for a query model) using a systematic and random component.  Both
    are taken to be symetric.

    Likely either of these could be 'None', but best set to '0' if
    that's what you really want to express.

    """

    def __init__(self, syser, rnder):
        self.syser = syser
        self.rnder = rnder

    # Properties
    @property
    def syser(self):
        """the systematic uncertainty in a coordinate transform.

        units: should be the same as the CoordAxis (s, arcsec, deg?)

        """
        return self._syser

    @syser.setter
    def syser(self, value):
        caom_util.type_check(value, float, "syser", override=False)
        self._syser = value

    @property
    def rnder(self):
        """the random uncertainty in a coordinate transform.

        units: should be the same as the CoordAxis transform.

        """
        return self._rnder

    @rnder.setter
    def rnder(self, value):
        caom_util.type_check(value, float, "rnder", override=False)
        self._rnder = value


class CoordFunction1D(common.CaomObject):
    """Defines a linear function that transforms from pixel to WCS
    values.

    """

    def __init__(self, naxis, delta, ref_coord):
        """
        Need to define the length of the axis, the slope of the
        conversion and a reference coordinate.  All are needed for a
        valid 1D function.

        """

        self.naxis = naxis
        self.delta = delta
        self.ref_coord = ref_coord

    @property
    def naxis(self):
        """The length of the axis.

        unit: pix
        type: int

        """
        return self._naxis

    @naxis.setter
    def naxis(self, value):
        caom_util.type_check(value, int, 'naxis', override=False)
        self._naxis = value

    @property
    def delta(self):
        """The step in WCS between pixels.

        unit: WCS/pix  (days if this is a timeWCS)
        type: float

        """
        return self._delta

    @delta.setter
    def delta(self, value):
        caom_util.type_check(value, float, 'delta', override=False)
        self._delta = value

    @property
    def ref_coord(self):
        """the (pix,val) reference for this transformtion.

        eg. ref_coord=RefCoord(pix,val)
        type: RefCoord

        """
        return self._ref_coord

    @ref_coord.setter
    def ref_coord(self, value):
        caom_util.type_check(value, RefCoord, 'ref_coord', override=False)
        self._ref_coord = value


class CoordFunction2D(common.CaomObject):
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
        caom_util.type_check(value, Dimension2D, 'dimension', override=False)
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
        caom_util.type_check(value, Coord2D, 'ref_coord', override=False)
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
        caom_util.type_check(value, float, 'cd11', override=False)
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
        caom_util.type_check(value, float, 'cd12', override=False)
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
        caom_util.type_check(value, float, 'cd21', override=False)
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
        caom_util.type_check(value, float, 'cd22', override=False)
        self._cd22 = value


class CoordPolygon2D(common.CaomObject):
    """A object to contain a TypeList ValueCoord2D vertices that are a
    polygon.  The vertices are given as ValueCoord2D objects, which are
    coordinate pairs.

    eg. vertices.add(ValueCoord2D(coord1,coord2))

    """

    def __init__(self, vertices=None):
        if vertices is None:
            vertices = caom_util.TypedList(ValueCoord2D, )
        self.vertices = vertices

    @property
    def vertices(self):
        """A TypedList of ValueCoord2D objects that layout the vertices of a
        polygon.

        A vertices can be added using the 'add' method..
        eg: vertices.add(ValueCoord2D())

        see the caom2.wcs.ValueCoord2D help for details on making a
        coordinate pair.

        type: TypedList((ValueCoord2D),)

        """
        return self._vertices

    @vertices.setter
    def vertices(self, value):
        caom_util.type_check(value, caom_util.TypedList, 'vertices',
                             override=False)
        self._vertices = value


class CoordRange1D(common.CaomObject):
    """a CoordRange1D object contains the start and end of
     a range of values, expressed in both pixel and WCS units.

     """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    @property
    def start(self):
        """The pixel and world coordinate of the start of a range.

        eg.  RefCoord(pix,val)
        type: RefCoord

        """
        return self._start

    @start.setter
    def start(self, value):
        caom_util.type_check(value, RefCoord, "start", override=False)
        self._start = value

    @property
    def end(self):
        """The pixel and world coordinate of the end of a range.

        eg. RefCoord(pix,val)
        type: RefCoord

        """
        return self._end

    @end.setter
    def end(self, value):
        caom_util.type_check(value, RefCoord, "end", override=False)
        self._end = value


class CoordRange2D(common.CaomObject):
    """A range (x1,y1) to (x2,y2) in two dimenstions.

    The range object should know the coordinate in both
    pixels and WCS units.

    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    @property
    def start(self):
        """The starting coordinate pair  of a range, in wcs units.

        eg: Coord2D(RefCoord(crpix1,crval1), RefCoord(crpix2,crval2))
        """
        return self._start

    @start.setter
    def start(self, value):
        caom_util.type_check(value, Coord2D, 'start', override=False)
        self._start = value

    @property
    def end(self):
        """The reference for the ending coordinate of a range.

        eg: Coord2D(RefCoord(crpix1,crval1), RefCoord(crpix2,crval2))
        """
        return self._end

    @end.setter
    def end(self, value):
        caom_util.type_check(value, Coord2D, 'end', override=False)
        self._end = value


class Dimension2D(common.CaomObject):
    """Hey, how big is this thing? What are its dimensions.  That's what
    Dimension2D will tell you.

    """

    def __init__(self, naxis1, naxis2):
        self.naxis1 = naxis1
        self.naxis2 = naxis2

    @property
    def naxis1(self):
        """The length of the first (x) dimension.

        eg.  naxis1=1024
        unit: pix
        type: int

        """
        return self._naxis1

    @naxis1.setter
    def naxis1(self, value):
        caom_util.type_check(value, int, 'naxis1', override=False)
        caom_util.value_check(value, 0, 1E10, 'naxis1', override=False)
        self._naxis1 = value

    @property
    def naxis2(self):
        """The length of the second (y) dimension.

        eg.  naxis2=2048
        unit: pix
        type: int

        """
        return self._naxis2

    @naxis2.setter
    def naxis2(self, value):
        caom_util.type_check(value, int, 'naxis2', override=False)
        caom_util.value_check(value, 0, 1E10, 'naxis2', override=False)
        self._naxis2 = value


class EnergyTransition(common.CaomObject):
    """ EnergyTransition """

    def __init__(self, species, transition):
        """
        Construct an EnergyTransition instance

        Arguments:
        species
        transition
        """
        caom_util.type_check(species, str, "species", override=False)
        caom_util.type_check(transition, str, "transition", override=False)
        self._species = species
        self._transition = transition

    @property
    def species(self):
        """ Species """
        return self._species

    @property
    def transition(self):
        """ Transition """
        return self._transition


class RefCoord(common.CaomObject):
    """A refernce coordinate object, maps pixel value to wcs value

    """

    def __init__(self, pix, val):
        """maps a pixel location to a wcs value, as a reference spot.

        eg.  RefCoord(crpix1, crval1)

        """
        self.pix = pix
        self.val = val

    @property
    def pix(self):
        """The pixel location of a reference position.

        units: pix
        type: float

        """
        return self._pix

    @pix.setter
    def pix(self, value):
        caom_util.type_check(value, float, 'pix', override=False)
        self._pix = value

    @property
    def val(self):
        """The WCS value at the reference position.

        units: CUNIT
        type: float

        """
        return self._val

    @val.setter
    def val(self, value):
        caom_util.type_check(value, float, 'val', override=False)
        self._val = value


class Slice(common.CaomObject):
    """defines a slice in a set of data contains values.

    The slice keeps track of the type/unit and values.

    ctype and cunit are stored in the axis variable.
    values are stored in the bins

    """

    def __init__(self, axis, bin):
        self.axis = axis
        self.bin = bin

    @property
    def axis(self):
        """A ctype/cunit pair for this slice of data.

        type: Axis(ctype,cunit)

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        caom_util.type_check(value, Axis, 'axis', override=False)
        self._axis = value

    @property
    def bin(self):
        """The pixel value on the axis.

        This value is use to transform to the WCS for this axis.
        unit: pixel
        type: int
        """
        return self._bin

    @bin.setter
    def bin(self, value):
        caom_util.type_check(value, int, 'int', override=False)
        self._bin = value


class ValueCoord2D(common.CaomObject):
    """Represents the reference point."""

    def __init__(self, coord1, coord2):
        self.coord1 = coord1
        self.coord2 = coord2

    @property
    def coord1(self):
        """Coordinate 1"""
        return self._coord1

    @coord1.setter
    def coord1(self, value):
        caom_util.type_check(value, float, 'coord1', override=False)
        self._coord1 = value

    @property
    def coord2(self):
        """Coordinate 2"""
        return self._coord2

    @coord2.setter
    def coord2(self, value):
        caom_util.type_check(value, float, 'coord2', override=False)
        self._coord2 = value
