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

"""defines the caom2.Plane class"""

from datetime import datetime

from builtins import str, int
from urllib.parse import SplitResult, urlsplit
from deprecated import deprecated

from caom2.caom_util import int_32
from . import caom_util
from . import shape
from . import wcs
from .artifact import Artifact
from .common import AbstractCaomEntity, CaomObject, ObservationURI,\
    VocabularyTerm, OrderedEnum
from .common import _CAOM_VOCAB_NS, _OBSCORE_VOCAB_NS
import warnings
from enum import Enum
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aenum import Enum, extend_enum

__all__ = ['CalibrationLevel', 'DataProductType', 'EnergyBand',
           'PolarizationState', 'Quality', 'Plane',
           'PlaneURI', 'DataQuality', 'Metrics', 'Provenance', 'Position',
           'Energy', 'Polarization', 'Time', 'Observable']


class CalibrationLevel(Enum):
    """
    PLANNED: -1
    RAW_INSTRUMENTAL: 0
    RAW_STANDARD: 1
    CALIBRATED: 2
    PRODUCT: 3
    """
    PLANNED = int_32(-1)
    RAW_INSTRUMENTAL = int_32(0)
    RAW_STANDARD = int_32(1)
    CALIBRATED = int_32(2)
    PRODUCT = int_32(3)
    ANALYSIS_PRODUCT = int_32(4)


class DataProductType(OrderedEnum):
    """ DataproductType - enum of data product types"""

    IMAGE = VocabularyTerm(_OBSCORE_VOCAB_NS, "image", True).get_value()
    CUBE = VocabularyTerm(_OBSCORE_VOCAB_NS, "cube", True).get_value()
    EVENTLIST = VocabularyTerm(_OBSCORE_VOCAB_NS, "eventlist",
                               True).get_value()
    SPECTRUM = VocabularyTerm(_OBSCORE_VOCAB_NS, "spectrum", True).get_value()
    TIMESERIES = VocabularyTerm(_OBSCORE_VOCAB_NS, "timeseries",
                                True).get_value()
    VISIBILITY = VocabularyTerm(_OBSCORE_VOCAB_NS, "visibility",
                                True).get_value()
    MEASUREMENTS = VocabularyTerm(_OBSCORE_VOCAB_NS, "measurements",
                                  True).get_value()
    CATALOG = VocabularyTerm(_CAOM_VOCAB_NS, "catalog").get_value()
    EVENT = VocabularyTerm(_CAOM_VOCAB_NS, "event", True).get_value()
    SED = VocabularyTerm(_CAOM_VOCAB_NS, "sed", True).get_value()

    @staticmethod
    def extend(namespace, name):
        """
        Extends the DataProductType with a new, user-defined, entry
        :param namespace: Namespace for the new data product type
        :param name: Name of the new data product type
        """
        extend_enum(DataProductType, name.upper(),
                    VocabularyTerm(namespace, name).get_value())


class EnergyBand(OrderedEnum):
    """
    RADIO: "Radio"
    MILLIMETER: "Millimeter"
    INFRARED: "Infrared"
    OPTICAL: "Optical"
    UV: "UV"
    EUV: "EUV"
    XRAY: "X-ray"
    GAMMARAY: "Gamma-ray"
    """
    RADIO = "Radio"
    MILLIMETER = "Millimeter"
    INFRARED = "Infrared"
    OPTICAL = "Optical"
    UV = "UV"
    EUV = "EUV"
    XRAY = "X-ray"
    GAMMARAY = "Gamma-ray"


class PolarizationState(OrderedEnum):
    """
    I: "I"
    Q: "Q"
    U: "U"
    V: "V"
    POLI: "POLI"
    FPOLI: "FPOLI"
    POLA: "POLA"
    EPOLI: "EPOLI"
    CPOLI: "CPOLI"
    NPOLI: "NPOLI"
    RR: "RR"
    LL: "LL"
    RL: "RL"
    LR: "LR"
    XX: "XX"
    YY: "YY"
    XY: "XY"
    YX: "YX"
    """
    I = "I"  # noqa
    Q = "Q"
    U = "U"
    V = "V"
    RR = "RR"
    LL = "LL"
    RL = "RL"
    LR = "LR"
    XX = "XX"
    YY = "YY"
    XY = "XY"
    YX = "YX"
    POLI = "POLI"
    FPOLI = "FPOLI"
    POLA = "POLA"
    EPOLI = "EPOLI"
    CPOLI = "CPOLI"
    NPOLI = "NPOLI"


class Quality(Enum):
    """
    JUNK: junk
    """
    JUNK = VocabularyTerm(_CAOM_VOCAB_NS, "junk", True).get_value()


class Observable():
    """ Observable class"""

    def __init__(self, ucd):
        self.ucd = ucd

    @property
    def ucd(self):
        return self._ucd

    @ucd.setter
    def ucd(self, value):
        caom_util.type_check(value, str, 'ucd', override=False)
        self._ucd = value


class Plane(AbstractCaomEntity):
    """ Plane class """

    def __init__(self, product_id,
                 creator_id=None,
                 artifacts=None,
                 meta_release=None,
                 data_release=None,
                 meta_read_groups=None,
                 data_read_groups=None,
                 data_product_type=None,
                 calibration_level=None,
                 provenance=None,
                 metrics=None,
                 quality=None,
                 observable=None):
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
        self.meta_read_groups = meta_read_groups
        self.data_read_groups = data_read_groups
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
        self._custom = None
        self.observable = observable

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
            if tmp.geturl() != value:
                raise ValueError("Invalid URI: " + value)
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
        caom_util.type_check(value, caom_util.TypedOrderedDict, 'artifacts',
                             override=False)
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
        caom_util.value_check(value, caom_util.MIN_DATETIME,
                              caom_util.MAX_DATETIME, "meta_release",)
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
        caom_util.value_check(value, caom_util.MIN_DATETIME,
                              caom_util.MAX_DATETIME, 'data_release')
        self._data_release = value

    @property
    def meta_read_groups(self):
        return self._meta_read_groups

    @meta_read_groups.setter
    def meta_read_groups(self, value):
        """
            value is a caom_util.URISet
        """
        if value is None:
            self._meta_read_groups = caom_util.URISet()
        else:
            caom_util.type_check(value, caom_util.URISet,
                                 'meta_read_groups')
            self._meta_read_groups = value

    @property
    def data_read_groups(self):
        return self._data_read_groups

    @data_read_groups.setter
    def data_read_groups(self, value):
        """
            value is a caom_util.URISet
        """
        if value is None:
            self._data_read_groups = caom_util.URISet()
        else:
            caom_util.type_check(value, caom_util.URISet,
                                 'data_read_groups')
            self._data_read_groups = value

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

    @property
    def observable(self):
        return self._observable

    @observable.setter
    def observable(self, value):
        caom_util.type_check(value, str, 'observable')
        self._observable = value

    @property
    def position(self):
        """A caom2 Position object that is developed from
        the aggregation of the Chunks that are children of
        the Plane.

        aggregation happens during ingest and is not part
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
        the aggregation of the Chunks that are children of
        the Plane.

        aggregation happens during ingest and is not part
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
        the aggregation of the Chunks that are children of
        the Plane.

        aggregation happens during ingest and is not part
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
        the aggregation of the Chunks that are children of
        the Plane.

        aggregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        caom_util.type_check(value, Polarization, "polarization")
        self._polarization = value

    @property
    def custom(self):
        """A caom2 Custom Axis object that is developed from
        the aggregation of the Chunks that are children of the Plane.

        aggregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._custom

    @custom.setter
    def custom(self, value):
        caom_util.type_check(value, CustomAxis, "custom")
        self._custom = value

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
            raise ValueError(
                'Cannot compare PlaneURI with {}'.format(type(other)))
        return self.uri < other.uri

    def __eq__(self, other):
        if not isinstance(other, PlaneURI):
            raise ValueError(
                'Cannot compare PlaneURI with {}'.format(type(other)))
        return self.uri == other.uri

    @classmethod
    def get_plane_uri(cls, observation_uri, product_id):
        """
        Initializes an Plane URI instance

        Arguments:
        observation_uri : the uri of the observation
        product_id : ID of the product
        """
        caom_util.type_check(observation_uri, ObservationURI,
                             "observation_uri",
                             override=False)
        caom_util.type_check(product_id, str, "observation_uri",
                             override=False)
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
        self._sample_snr = None

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

    @property
    def sample_snr(self):
        """
        TBD
        """
        return self._sample_snr

    @sample_snr.setter
    def sample_snr(self, value):
        caom_util.type_check(value, float, 'sample_snr')
        caom_util.value_check(value, 0, 1E10, 'mag_limit')  # TODO limits?
        self._sample_snr = value


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

        if not name:
            raise AttributeError("No name provided")
        caom_util.type_check(name, str, 'name', override=False)
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
            if tmp.geturl() != value:
                raise ValueError("Invalid URI: " + value)
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
                 resolution_bounds=None,
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
        self.resolution_bounds = resolution_bounds
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
            caom_util.type_check(value,
                                 (shape.Box, shape.Circle, shape.Polygon),
                                 'bounds', override=False)
        self._bounds = value

    @property
    def dimension(self):
        """ Dimension """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            caom_util.type_check(value, wcs.Dimension2D,
                                 'dimension', override=False)
        self._dimension = value

    @property
    def resolution(self):
        """ Resolution """
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        if value is not None:
            caom_util.type_check(value, float, 'resolution')
        self._resolution = value

    @property
    def resolution_bounds(self):
        """ Resolution bounds"""
        return self._resolution_bounds

    @resolution_bounds.setter
    def resolution_bounds(self, value):
        if value is not None:
            caom_util.type_check(value, shape.Interval, 'resolution bounds')
        self._resolution_bounds = value

    @property
    def sample_size(self):
        """ Sample size """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        if value is not None:
            caom_util.type_check(value, float, 'sample_size')
        self._sample_size = value

    @property
    def time_dependent(self):
        """ Time dependent """
        return self._time_dependent

    @time_dependent.setter
    def time_dependent(self, value):
        if value is not None:
            caom_util.type_check(value, bool, 'time_dependent')
        self._time_dependent = value


class Energy(CaomObject):
    """ Energy """

    def __init__(self, bounds=None, dimension=None, resolving_power=None,
                 resolving_power_bounds=None, energy_bands=None,
                 sample_size=None, bandpass_name=None, em_band=None,
                 transition=None, restwav=None):
        """
        Initialize an Energy instance.

        Arguments:
        None
        """
        self.bounds = bounds
        self.dimension = dimension
        self.resolving_power = resolving_power
        self.resolving_power_bounds = resolving_power_bounds
        self.sample_size = sample_size
        self.bandpass_name = bandpass_name
        self.energy_bands = energy_bands
        self.em_band = em_band
        if em_band is not None:
            self.energy_bands.add(em_band)
        self.transition = transition
        self.restwav = restwav

    # Properties

    @property
    def bounds(self):
        """ Energy bounds """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            caom_util.type_check(value, shape.Interval, 'bounds')
        self._bounds = value

    @property
    def dimension(self):
        """DIMENSION (NUMBER OF PIXELS) ALONG ENERGY AXIS."""
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            caom_util.type_check(value, int, 'dimension')
        self._dimension = value

    @property
    def resolving_power(self):
        """ Resolving power """
        return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        if value is not None:
            caom_util.type_check(value, float, 'resolving_power')
        self._resolving_power = value

    @property
    def resolving_power_bounds(self):
        """ Resolving power bounds"""
        return self._resolving_power_bounds

    @resolving_power_bounds.setter
    def resolving_power_bounds(self, value):
        if value is not None:
            caom_util.type_check(value, shape.Interval,
                                 'resolving power bounds')
        self._resolving_power_bounds = value

    @property
    def sample_size(self):
        """ Sample size """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        if value is not None:
            caom_util.type_check(value, float, 'sample_size')
        self._sample_size = value

    @property
    def energy_bands(self):
        """Energy bands"""
        return self._energy_bands

    @energy_bands.setter
    def energy_bands(self, value):
        if value is not None:
            if not isinstance(value, caom_util.TypedSet) and \
                            value.oktypes != EnergyBand:
                raise TypeError('Energy bands must be of type '
                                'caom_util.TypedSet(EnergyBand)')
            self._energy_bands = value
        else:
            self._energy_bands = caom_util.TypedSet(EnergyBand)

    @property
    def bandpass_name(self):
        """ Bandpass name """
        return self._bandpass_name

    @bandpass_name.setter
    def bandpass_name(self, value):
        if value is not None:
            caom_util.type_check(value, str, 'bandpass_name')
        self._bandpass_name = value

    @property
    @deprecated(version='CAOM2.4',
                reason='Replaced by energy_bands, gone in 2.5')
    def em_band(self):
        """ EM Band """
        return None

    @em_band.setter
    @deprecated(version='CAOM2.4',
                reason='Replaced by energy_bands, gone in 2.5')
    def em_band(self, value):
        pass

    @property
    def transition(self):
        """ Energy transition """
        return self._transition

    @transition.setter
    def transition(self, value):
        if value is not None:
            caom_util.type_check(value, wcs.EnergyTransition, 'transition')
        self._transition = value

    @property
    def restwav(self):
        """ rest wavelength of the target energy transition """
        return self._restwav

    @restwav.setter
    def restwav(self, value):
        if value is not None:
            caom_util.type_check(value, float, 'restwav')
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
            caom_util.type_check(value, list, 'polarization_states',
                                 override=False)
        self._polarization_states = value


class Time(CaomObject):
    """ Time """

    def __init__(self,
                 bounds=None,
                 dimension=None,
                 resolution=None,
                 resolution_bounds=None,
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
        self.resolution_bounds = resolution_bounds
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
    def resolution_bounds(self):
        """ Resolution bounds"""
        return self._resolution_bounds

    @resolution_bounds.setter
    def resolution_bounds(self, value):
        if value is not None:
            caom_util.type_check(value, shape.Interval, 'resolution bounds')
        self._resolution_bounds = value

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


class CustomAxis(CaomObject):
    """
    Custom Axis == NonStandard Axis, where Standard Axis is either
     Position, Energy, Polarization or Time
     """

    def __init__(self,
                 ctype,
                 bounds=None,
                 dimension=None):
        """
        Initialize a Custom Axis instance.
        """
        if ctype is None:
            raise AttributeError('ctype of CustomAxis cannot be None')
        self._ctype = ctype
        self.bounds = bounds
        self.dimension = dimension

    # Properties

    @property
    def ctype(self):
        return self._ctype

    @property
    def bounds(self):
        """
        An interval object that gives start and end of a custom interval
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        caom_util.type_check(value, shape.Interval, 'bounds')
        self._bounds = value

    @property
    def dimension(self):
        """
        Number of pixel in the custom direction, normally 1.
        """
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        caom_util.type_check(value, int, 'dimension')
        self._dimension = value
