# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2016.                            (c) 2016.
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

"""defines the caom2.Plane class"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from datetime import datetime
from six.moves.urllib.parse import SplitResult, urlsplit
from builtins import str, int
from enum import Enum

from . import caom_util
from . import shape
from . import wcs
from .artifact import Artifact
from .common import AbstractCaomEntity
from .common import CaomObject
from .common import ObservationURI
from caom2.caom_util import int_32

__all__ = ['CalibrationLevel', 'DataProductType', 'EnergyBand', 
           'VocabularyTerm', 'PolarizationState', 'Quality', 'Plane', 
           'PlaneURI', 'DataQuality', 'Metrics', 'Provenance', 'Position', 
           'Energy', 'Polarization', 'Time']


class CalibrationLevel(Enum):
    """
    PLANNED: -1
    RAW_INSTRUMENT: 0
    RAW_STANDARD: 1
    CALIBRATED: 2
    PRODUCT: 3
    """
    PLANNED = int_32(-1)
    RAW_INSTRUMENT = int_32(0)
    RAW_STANDARD = int_32(1)
    CALIBRATED = int_32(2)
    PRODUCT = int_32(3)
    ANALYSIS_PRODUCT = int_32(4)


class VocabularyTerm(object):
    """ VocabularyTerm """

    def __init__(self, namespace, term, base=False):
        """
        Construct a VocabularyTerm instance. This creates a term in the
        specified vocabulary namespace. If the value of base is False, 
        the string value (from getvalue()) will just be the namespace URI
        plus the term added as a fragment. If the value of base is True, 
        this is a term in a base vocabulary and the value will just be the
        term (without the namespace).

        Arguments:
        namespace : namespace of the vocabulary
        term : a term in the base vocabulary
        base : if True, getValue() returns term, otherwise getvalue() returns
               namespace URI plus term
        """
        self.namespace = namespace
        self.term = term
        self.base = base

    def get_value(self):
        """ get_value """
        if self.base:
            return self._term
        else:
            return self._namespace + "#" + self._term

    def __str__(self):
        """ __str__ """
        return self.get_value()

    # Properties
    @property
    def namespace(self):
        """ namespace """
        return self._namespace

    @namespace.setter
    def namespace(self, value):
        caom_util.type_check(value, str, "namespace")
        tmp = urlsplit(value)
        assert tmp.geturl() == value, "Invalid URI: " + value
        self._namespace = value

    @property
    def term(self):
        """ term """
        return self._term

    @term.setter
    def term(self, value):
        caom_util.type_check(value, str, "term")
        self._term = value

    @property
    def base(self):
        """ base """
        return self._base

    @base.setter
    def base(self, value):
        caom_util.type_check(value, bool, "base")
        self._base = value


class DataProductType(Enum):
    """ DataproductType """

    _OBSCORE = "http://www.ivoa.net/std/ObsCore"
    _CAOM = "http://www.opencadc.org/caom2/DataProductType"

    """ 
    def __init__(self, value, namespace=None):
        Initialize a DataProductType instance 

        Arguments:
        value : fragment to be appended to namespace
        namespace : namespace the data product type belongs to,
                    defaults to OBSCORE
        if namespace == None:
            super(DataProductType, self).__init__(DataProductType._OBSCORE, value, True)
        else:
            super(DataProductType, self).__init__(namespace, value)
    """

    """
    ObsCore-1.0
    IMAGE: "image"
    CATALOG: "catalog"
    CUBE: "cube"
    EVENTLIST: "eventlist"
    SPECTRUM: "spectrum"
    TIMESERIES: "timeseries"
    VISIBILITY: "visibility"

    ObsCore-1.1
    MEASUREMENTS: "measurements"

    ObsCore-2.3
    CATALOG: "http://www.opencadc.org/caom2#catalog"
    """

    IMAGE = VocabularyTerm(_OBSCORE, "image", True).get_value()
    CUBE = VocabularyTerm(_OBSCORE, "cube", True).get_value()
    EVENTLIST = VocabularyTerm(_OBSCORE, "eventlist", True).get_value()
    SPECTRUM = VocabularyTerm(_OBSCORE, "spectrum", True).get_value()
    TIMESERIES = VocabularyTerm(_OBSCORE, "timeseries", True).get_value()
    VISIBILITY = VocabularyTerm(_OBSCORE, "visibility", True).get_value()
    MEASUREMENTS = VocabularyTerm(_OBSCORE, "measurements", True).get_value()
    CATALOG = VocabularyTerm(_CAOM, "catalog").get_value()


class EnergyBand(Enum):
    """
    GAMMARAY: "Gamma-ray"
    INFRARED: "Infrared"
    MILLIMETER: "Millimeter"
    OPTICAL: "Optical"
    RADIO: "Radio"
    UV: "UV"
    XRAY: "X-ray"
    """
    EUV = "EUV"
    GAMMARAY = "Gamma-ray"
    INFRARED = "Infrared"
    MILLIMETER = "Millimeter"
    OPTICAL = "Optical"
    RADIO = "Radio"
    UV = "UV"
    XRAY = "X-ray"


class PolarizationState(Enum):
    """
    I: "I"
    Q: "Q"
    U: "U"
    V: "V"
    LL: "LL"
    LR: "LR"
    RL: "RL"
    RR: "RR"
    XX: "XX"
    XY: "XY"
    YX: "YX"
    YY: "YY"
    """
    I = "I"
    Q = "Q"
    U = "U"
    V = "V"
    LL = "LL"
    LR = "LR"
    RL = "RL"
    RR = "RR"
    XX = "XX"
    XY = "XY"
    YX = "YX"
    YY = "YY"


class Quality(Enum):
    """
    JUNK: junk
    """
    JUNK = "junk"


class Plane(AbstractCaomEntity):
    """ Plane class """

    def __init__(self, product_id,
                 creator_id=None,
                 artifacts=None,
                 meta_release=None,
                 data_release=None,
                 data_product_type=None,
                 calibration_level=None,
                 provenance=None,
                 metrics=None,
                 quality=None):
        """
        Initialize a Plane instance

        Arguments:
        product_id : product ID
        """
        super(Plane, self).__init__()
        self.product_id = product_id
        if artifacts is None:
            artifacts = caom_util.TypedOrderedDict(Artifact, )
        self.creator_id = creator_id
        self.artifacts = artifacts

        self.meta_release = meta_release
        self.data_release = data_release
        self.data_product_type = data_product_type
        self.calibration_level = calibration_level
        self.provenance = provenance
        self.metrics = metrics
        self.quality = quality

        # computed fields
        # aggregated from the Chunks during ingestion
        self._position = None
        self._energy = None
        self._time = None
        self._polarization = None

    def _key(self):
        return self.product_id

    def __hash__(self):
        return hash(self._key())

    # Properties
    @property
    def product_id(self):
        """A string that identifies the data product, within a given
        observation, that is stored in this plane.

        eg: '1234567p'
        type: unicode string
        """
        return self._product_id

    @product_id.setter
    def product_id(self, value):
        caom_util.type_check(value, str, 'product_id', override=False)
        self._product_id = value

    @property
    def creator_id(self):
        """A URI that identifies the creator of this plane.

        eg: ivo://cadc.nrc.ca/users?tester
        type: URI 
        """
        return self._creator_id

    @creator_id.setter
    def creator_id(self, value):
        caom_util.type_check(value, str, 'creator_id')
        if value is not None:
            tmp = urlsplit(value)
            assert tmp.geturl() == value, "Invalid URI: " + value
        self._creator_id = value

    @property
    def artifacts(self):
        """A TypeList of artifacts that are part of this plane.

        individual artifacts are constructed and then added to the plane.

        eg. Plane.artifacts.add(Artifact('ad:CFHT/1234567p')), see the
        Arifact help for more info on making an Aritfact object return
        """
        return self._artifacts

    @artifacts.setter
    def artifacts(self, value):
        caom_util.type_check(value, caom_util.TypedOrderedDict, 'artifacts', override=False)
        self._artifacts = value

    @property
    def meta_release(self):
        """The date when the metadata describing this observation become
        public.

        eg.  Plane.meta_releaes=datetime.datetime(2012,1,1,0,0,0)
        indicates that the metadata become public at 00:00:00 on
        January 1st 2012.

        if there is no meta_release period, set value to None

        unit: calendar date
        type: datetime
        """
        return self._meta_release

    @meta_release.setter
    def meta_release(self, value):
        caom_util.type_check(value, datetime, 'meta_release')
        caom_util.value_check(value,
                              datetime(1800, 1, 1, 0, 0, 0),
                              datetime(2050, 1, 1, 0, 0, 0),
                              'meta_release')
        self._meta_release = value

    @property
    def data_release(self):
        """The date when the data contained in this plane become public.

        eg.  Plane.data_releaes=datetime.datetime(2012,1,1,0,0,0)
        indicates that ar 00:00:00 on January 1st 2012 the data
        assocaited with this Plane will become publicly accessible.

        If there is no data_release period, set value to None.

        unit: calendar date
        type: datetime
        """
        return self._data_release

    @data_release.setter
    def data_release(self, value):
        caom_util.type_check(value, datetime, 'data_release')
        caom_util.value_check(value,
                              datetime(1800, 1, 1, 0, 0, 0),
                              datetime(2050, 1, 1, 0, 0, 0),
                              'data_release')
        self._data_release = value

    @property
    def data_product_type(self):
        """The type of file structure that this plane contains.

        eg.
        Plane.data_product_type = DataProductType['EVENTLIST']

        see DataProductType for allowed values

        """
        return self._data_product_type

    @data_product_type.setter
    def data_product_type(self, value):
        caom_util.type_check(value, DataProductType, 'data_product_type')
        self._data_product_type = value

    @property
    def calibration_level(self):
        """a string that represents the level of calibration (aka processing)
        the data contained in this plane have received.  The string
        is converted to an integer during storage.

        eg. Plane.calibration_level = CalibrationLevel['RAW_STANDARD']
        type: unicode string

        Must be one of CalibrationLevel

        """
        return self._calibration_level

    @calibration_level.setter
    def calibration_level(self, value):
        caom_util.type_check(value, CalibrationLevel, "calibration_level")
        self._calibration_level = value

    @property
    def provenance(self):
        """The process that created the data referred to by this Plane.

        eg.  Plane.provenance=caom2.Provenance("Elixir")
        """
        return self._provenance

    @provenance.setter
    def provenance(self, value):
        caom_util.type_check(value, Provenance, "provenance")
        self._provenance = value

    @property
    def metrics(self):
        """reference to an object that contains metrics of this plane.

        eg. Plane.metrics = caom2.Metrics()
        """
        return self._metrics

    @metrics.setter
    def metrics(self, value):
        caom_util.type_check(value, Metrics, 'metrics')
        self._metrics = value

    @property
    def quality(self):
        """reference to an object that describes the quality of the data of this plane.

        eg. Plane.data_quality = caom2.DataQuality()
        """
        return self._quality

    @quality.setter
    def quality(self, value):
        caom_util.type_check(value, DataQuality, 'quality')
        self._quality = value

    # @property
    # def observable(self):
    #    """ """
    #    return self._observable

    @property
    def position(self):
        """A caom2 Position object that is developed from
        the agregation of the Chunks that are children of
        the Plane.

        agregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._position
    
    @position.setter
    def position(self, value):
        caom_util.type_check(value, Position, "position")
        self._position = value

    @property
    def energy(self):
        """A caom2 Energy object that is developed from
        the agregation of the Chunks that are children of
        the Plane.

        agregation happens during ingest and is not part
        of the python module at this time.
        """
        """ Energy """
        return self._energy
    
    @energy.setter
    def energy(self, value):
        caom_util.type_check(value, Energy, "energy")
        self._energy = value

    @property
    def time(self):
        """A caom2 Time object that is developed from
        the agregation of the Chunks that are children of
        the Plane.

        agregation happens during ingest and is not part
        of the python module at this time.
        """
        """ Time """
        return self._time
    
    @time.setter
    def time(self, value):
        caom_util.type_check(value, Time, "time")
        self._time = value
    
    @property
    def polarization(self):
        """A caom2 Polarization object that is developed from
        the agregation of the Chunks that are children of
        the Plane.

        agregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._polarization
    
    @polarization.setter
    def polarization(self, value):
        caom_util.type_check(value, Polarization, "polarization")
        self._polarization = value

    # Compute derived fields

    def compute_position(self):
        raise NotImplementedError(
            "Aggregation of position has not been implemented in this module")

    def compute_energy(self):
        raise NotImplementedError(
            "Aggregation of energy has not been implemented in this module")

    def compute_time(self):
        raise NotImplementedError(
            "Aggregation of time has not been implemented in this module")

    def compute_polarization(self):
        raise NotImplementedError(
            "Aggregation of polarization " +
            "has not been implemented in this module")


class PlaneURI(CaomObject):
    """ Plane URI """

    def __init__(self, uri):
        """
        Initializes an Plane instance

        Arguments:
        uri : URI corresponding to the plane

        Throws:
        TypeError  : if uri is not a string
        ValueError : if uri is invalid
        ValueError : if the uri is valid but does not contain the expected
        fields (collection, observation_id and product_id)
        """

        self.uri = uri

    def _key(self):
        return self.uri

    def __hash__(self):
        return hash(self._key())

    def __lt__(self, other):
        if not isinstance(other, PlaneURI):
            raise ValueError('Canot compare PlaneURI with {}'.format(type(other)))
        return self.uri < other.uri

    def __eq__(self, other):
        if not isinstance(other, PlaneURI):
            raise ValueError('Canot compare PlaneURI with {}'.format(type(other)))
        return self.uri == other.uri

    @classmethod
    def get_plane_uri(cls, observation_uri, product_id):
        """
        Initializes an Plane URI instance

        Arguments:
        observation_uri : the uri of the observation
        product_id : ID of the product
        """
        caom_util.type_check(observation_uri, ObservationURI, "observation_uri",
                             override=False)
        caom_util.type_check(product_id, str, "observation_uri", override=False)
        caom_util.validate_path_component(cls, "product_id", product_id)

        path = urlsplit(observation_uri.uri).path
        uri = SplitResult(ObservationURI._SCHEME, "", path + "/" +
                          product_id, "", "").geturl()
        return cls(uri)

    # Properties
    @property
    def uri(self):
        """A uri that locates the plane object inside caom"""
        return self._uri

    @uri.setter
    def uri(self, value):

        caom_util.type_check(value, str, "uri", override=False)
        tmp = urlsplit(value)

        if tmp.scheme != ObservationURI._SCHEME:
            raise ValueError("{} doesn't have an allowed scheme".format(value))
        if tmp.geturl() != value:
            raise ValueError("Failed to parse uri correctly: {}".format(value))

        (collection, observation_id, product_id) = tmp.path.split("/")

        if product_id is None:
            raise ValueError("Faield to get product ID from uri: {}"
                             .format(value))

        self._product_id = product_id
        self._observation_uri = \
            ObservationURI.get_observation_uri(collection, observation_id)
        self._uri = value

    def get_product_id(self):
        """return the product_id associated with this plane"""
        return self._product_id

    def get_observation_uri(self):
        """Return the uri that can be used to find the caom2 observation object that
        this plane belongs to"""
        return self._observation_uri


class DataQuality(CaomObject):
    """ DataQuality """

    def __init__(self, flag):
        """
        Construct an DataQuality instance

        Arguments:
        flag
        """
        self.flag = flag

    @property
    def flag(self):
        """ flag """
        return self._flag

    @flag.setter
    def flag(self, value):
        caom_util.type_check(value, Quality, "flag")
        self._flag = value


class Metrics(CaomObject):
    """ Metrics """

    def __init__(self):
        """
        Initializes a Metrics instance

        Arguments:
        None
        """
        self._source_number_density = None
        self._background = None
        self._background_std_dev = None
        self._flux_density_limit = None
        self._mag_limit = None

    # Properties
    @property
    def source_number_density(self):
        """The number of sources brighter than mag_limit (flux_density_limit)

        unit: ct/deg2
        type: float
        """
        return self._source_number_density

    @source_number_density.setter
    def source_number_density(self, value):
        caom_util.type_check(value, float, "source_number_density")
        caom_util.value_check(value, 0, 1E10, "source_number_density")
        self._source_number_density = value

    @property
    def background(self):
        """The flux in the sky (background).

        units: Jy/pix
        type: float
        """
        return self._background

    @background.setter
    def background(self, value):
        caom_util.type_check(value, float, "background")
        caom_util.value_check(value, 0, 1E10, "background")
        self._background = value

    @property
    def background_std_dev(self):
        """the standard deviation (per pixel) in background flux.

        Likely this only makes sense to define if background is also defined.
        units: Jy/pix
        type: float
        """
        return self._background_std_dev

    @background_std_dev.setter
    def background_std_dev(self, value):
        caom_util.type_check(value, float, "background_std_dev")
        caom_util.value_check(value, 0, 1E10, "background")
        self._background_std_dev = value

    @property
    def flux_density_limit(self):
        """flux density where S:N=5 for point source.

        this is intended to provide a measure of the limit of detection.

        units: Jy
        type: float
        """
        return self._flux_density_limit

    @flux_density_limit.setter
    def flux_density_limit(self, value):
        caom_util.type_check(value, float, "flux_denisty_limit")
        caom_util.value_check(value, 0, 1E10, "flux_density_limit")
        self._flux_density_limit = value

    @property
    def mag_limit(self):
        """AB magnitude limit where S:N=5 for point source.

        Likely specify just mag_limit or flux_density_limit, not both?

        units: AB mag
        type: float
        """
        return self._mag_limit

    @mag_limit.setter
    def mag_limit(self, value):
        caom_util.type_check(value, float, 'mag_limit')
        caom_util.value_check(value, 0, 40, 'mag_limit')
        self._mag_limit = value


class Provenance(CaomObject):
    """ Provenance """

    def __init__(self, name,
                 version=None,
                 project=None,
                 producer=None,
                 run_id=None,
                 reference=None,
                 last_executed=None):
        """
        Initializes a Provenance instance

        Arguments:
        name - name of the provenance
        """

        assert name is not None, "No name provided"
        assert isinstance(name, str), "name is not a unicode string: {0}".format(name)
        self._name = name

        self.version = version
        self.project = project
        self.producer = producer
        self.run_id = run_id
        self.reference = reference
        self.last_executed = last_executed

        self._keywords = set()
        self._inputs = caom_util.TypedSet(PlaneURI, )

    # Properties

    @property
    def name(self):
        """ Name """
        return self._name

    @property
    def version(self):
        """ Version """
        return self._version

    @version.setter
    def version(self, value):
        caom_util.type_check(value, str, 'version')
        self._version = value

    @property
    def project(self):
        """ Project """
        return self._project

    @project.setter
    def project(self, value):
        caom_util.type_check(value, str, 'project')
        self._project = value

    @property
    def producer(self):
        """ Producer """
        return self._producer

    @producer.setter
    def producer(self, value):
        caom_util.type_check(value, str, 'producer')
        self._producer = value

    @property
    def run_id(self):
        """ Run ID """
        return self._run_id

    @run_id.setter
    def run_id(self, value):
        caom_util.type_check(value, str, 'run_id')
        self._run_id = value

    @property
    def reference(self):
        """ Reference """
        return self._reference

    @reference.setter
    def reference(self, value):
        caom_util.type_check(value, str, 'version')
        if value is not None:
            tmp = urlsplit(value)
            assert tmp.geturl() == value, "Invalid URI: " + value
        self._reference = value

    @property
    def last_executed(self):
        """ Version """
        return self._last_executed

    @last_executed.setter
    def last_executed(self, value):
        caom_util.type_check(value, datetime, 'last_executed')
        self._last_executed = value

    @property
    def keywords(self):
        """ Set of keywords as unicode string"""
        return self._keywords

    @property
    def inputs(self):
        """ Set of inputs as PlaneURI"""
        return self._inputs


class Position(CaomObject):
    """ Position """

    def __init__(self, bounds=None,
                 dimension=None,
                 resolution=None,
                 sample_size=None,
                 time_dependent=None
                 ):
        """
        Initialize a Position instance.

        Arguments:
        None
        """
        self.bounds = bounds
        self.dimension = dimension
        self.resolution = resolution
        self.sample_size = sample_size
        self.time_dependent = time_dependent

    # Properties

    @property
    def bounds(self):
        """ Bounds """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            assert isinstance(value, (shape.Box, shape.Circle, shape.Polygon)), (
                "bounds is not a Shape: {0}".format(value))
        self._bounds = value

    @property
    def dimension(self):
        """ Dimension """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            assert isinstance(value, wcs.Dimension2D), (
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


class Energy(CaomObject):
    """ Energy """

    def __init__(self):
        """
        Initialize an Energy instance.

        Arguments:
        None
        """
        self._bounds = None
        self._dimension = None
        self._resolving_power = None
        self._sample_size = None
        self._bandpass_name = None
        self._em_band = None
        self._transition = None
        self._restwav = None

    # Properties

    @property
    def bounds(self):
        """ Energy bounds """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            assert isinstance(value, shape.Interval), (
                "energy bounds is not an Interval: {0}".format(value))
        self._bounds = value

    @property
    def dimension(self):
        """DIMENSION (NUMBER OF PIXELS) ALONG ENERGY AXIS."""
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            assert isinstance(value, int), (
                "energy dimension is not an int: {0}".format(value))
        self._dimension = value

    @property
    def resolving_power(self):
        """ Resolving power """
        return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        if value is not None:
            assert isinstance(value, float), (
                "resolving power is not a float: {0}".format(value))
        self._resolving_power = value

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
    def bandpass_name(self):
        """ Bandpass name """
        return self._bandpass_name

    @bandpass_name.setter
    def bandpass_name(self, value):
        if value is not None:
            assert isinstance(value, str), (
                "bandpass name is not unicode string: {0}".format(value))
        self._bandpass_name = value

    @property
    def em_band(self):
        """ EM Band """
        return self._em_band

    @em_band.setter
    def em_band(self, value):
        if value is not None:
            assert isinstance(value, EnergyBand), (
                "em_Band is not an EnergyBand: {0}".format(value))
        self._em_band = value

    @property
    def transition(self):
        """ Energy transition """
        return self._transition

    @transition.setter
    def transition(self, value):
        if value is not None:
            assert isinstance(value, wcs.EnergyTransition), (
                "transition is not an EnergyTransition: {0}".format(value))
        self._transition = value

    @property
    def restwav(self):
        """ rest wavelength of the target energy transition """
        return self._restwav

    @restwav.setter
    def restwav(self, value):
        if value is not None:
            assert isinstance(value, float), (
                "restwav is not a float: {0}".format(value))
        self._restwav = value


class Polarization(CaomObject):
    """ Polarization """

    def __init__(self,
                 dimension=None,
                 polarization_states=None):
        """
        Initialize a Polarization instance.

        Arguments:
        None
        """
        self.dimension = dimension
        self.polarization_states = polarization_states

    # Properties
    @property
    def dimension(self):
        """number of samples (pixels) along polarization axis.

        unit: pix
        type: int
        """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        caom_util.type_check(value, int, 'dimension')
        caom_util.value_check(value, 0, 1E10, 'dimension')
        self._dimension = value
        
    @property
    def polarization_states(self):
        """
        type: list
        """
        return self._polarization_states

    @polarization_states.setter
    def polarization_states(self, value):
        if value is not None:
            caom_util.type_check(value, list, 'polarization_states', override=False)
        self._polarization_states = value


class Time(CaomObject):
    """ Time """

    def __init__(self,
                 bounds=None,
                 dimension=None,
                 resolution=None,
                 sample_size=None,
                 exposure=None):
        """
        Initialize a Time instance.

        Arguments:
        None
        """
        self.bounds = bounds
        self.dimension = dimension
        self.resolution = resolution
        self.sample_size = sample_size
        self.exposure = exposure

    # Properties

    @property
    def bounds(self):
        """an interval object that gives start and end of an time interval

        Not actually implemented yet.  Likely you want a TemporalWCS
        instead

        type: Interval(lower_mjd, upper_mjd)
        unit: mjd

        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        caom_util.type_check(value, shape.Interval, 'bounds')
        self._bounds = value

    @property
    def dimension(self):
        """Number of pixel in the time direction, normally 1.

        eg 1
        type: int

        """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        caom_util.type_check(value, int, 'dimension')
        self._dimension = value

    @property
    def resolution(self):
        """Time resolution of the samples, in seconds.

        normally this is the same as the exposure time,
        but in a stack you might have a larger resolution value than
        exposure time.

        eg. 1000
        unit: s
        type: float

        """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        caom_util.type_check(value, float, 'resolution')
        self._resolution = value

    @property
    def sample_size(self):
        """nominally the exposure time, in seconds.

        """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        caom_util.type_check(value, float, 'sample_size')
        self._sample_size = value

    @property
    def exposure(self):
        """Duration of the exposure, in seconds"""
        return self._exposure

    @exposure.setter
    def exposure(self, value):
        caom_util.type_check(value, float, 'exposure')
        self._exposure = value
