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

"""definition of the  caom2.Observation object."""

from datetime import datetime

from builtins import str
import warnings
from deprecated import deprecated

from . import caom_util
from .caom_util import int_32
from .common import AbstractCaomEntity, CaomObject, ObservationURI, \
    VocabularyTerm, OrderedEnum
from .common import _CAOM_VOCAB_NS
from .plane import Plane
from .shape import Point
from urllib.parse import urlsplit
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aenum import Enum

__all__ = ['ObservationIntentType', 'Status', 'TargetType',
           'Observation', 'ObservationURI', 'Algorithm', 'SimpleObservation',
           'DerivedObservation', 'Environment', 'Instrument', 'Proposal',
           'Requirements', 'Target', 'TargetPosition', 'Telescope',
           'CompositeObservation']


class ObservationIntentType(OrderedEnum):
    """
    CALIBRATION: "calibration"
    SCIENCE: "science"
    """
    SCIENCE = "science"
    CALIBRATION = "calibration"


class Status(Enum):
    """
    FAIL: "fail"
    """
    FAIL = VocabularyTerm(_CAOM_VOCAB_NS, "fail", True).get_value()


class TargetType(Enum):
    """
    FIELD: "field",
    OBJECT: "object"
    """
    FIELD = VocabularyTerm(_CAOM_VOCAB_NS, "field", True).get_value()
    OBJECT = VocabularyTerm(_CAOM_VOCAB_NS, "object", True).get_value()


class Observation(AbstractCaomEntity):
    """
    Observation

    An observation describes a set of empirical data.

    The CAOM2 observation is described at:
    http://www.cadc.hia.nrc.gc.ca/caom2

    The Observation object is the top level container for the
    meta-data assocaited with an observation.  The variaous attributes
    described here provide ways of storing the complete set of
    associated meta-data as well as links to the actual observation
    data.

    Observation -> Target
                -> Plane(s)
                -> Instrument
                -> Telescope
                -> Proposal
                -> Environment
                -> TargetPosition

    Plane -> Artifact(s)

    Artifact -> Part(s)

    Part -> Chunk(s)

    The actual 'chunks' of observational data are reference by Chunk
    objects.  Information about the Spatial/Frequency/Time aspects of
    an Observation are expressed at the Chunk level.

    The Chunk contains references to caom2 objects that fully describe
    the circumstances of that chunk of observation.  Often a 'Chunk'
    is a single extension in a FITS image. But can also be a
    particular column of values from a FITS Table or some other data
    object.

    Chunk -> SpatialWCS
          -> TemporalWCS
          -> SpectralWCS
          -> PolarizationWCS
          -> (Observable)


    """

    def __init__(self,
                 collection,
                 observation_id,
                 algorithm,
                 sequence_number=None,
                 intent=None,
                 type=None,
                 proposal=None,
                 telescope=None,
                 instrument=None,
                 target=None,
                 meta_release=None,
                 meta_read_groups=None,
                 planes=None,
                 environment=None,
                 target_position=None,
                 requirements=None
                 ):
        """
        Initializes an Observation instance

        Arguments: collection : where the observation is from
        (eg. 'HST')

        observation_id : a unique identifier within that collection
                         (eg. '111')

        algorithm : the algorithm used to create the observation.  For
                    a telescope observation this is always 'exposure'
        """
        super(Observation, self).__init__()

        self.collection = collection
        self.observation_id = observation_id
        if not algorithm:
            raise AttributeError('Algorithm required')
        self.algorithm = algorithm

        self._uri = ObservationURI.get_observation_uri(collection,
                                                       observation_id)

        self.sequence_number = sequence_number
        self.intent = intent
        self.type = type
        self.proposal = proposal
        self.telescope = telescope
        self.instrument = instrument
        self.target = target
        self.environment = environment
        self.target_position = target_position
        self.requirements = requirements
        self.meta_release = meta_release
        self.meta_read_groups = meta_read_groups
        self.planes = planes

    # Properties
    @property
    def collection(self):
        """The name of the collection of observations, normally a telescope
        name.

        type: unicode string
        """
        return self._collection

    @collection.setter
    def collection(self, value):
        caom_util.type_check(value, str, 'collection', override=False)
        self._collection = value

    @property
    def observation_id(self):
        """A string that uniquely identifies this obseravtion within the given
        collection.

        type: unicode string
        """
        return self._observation_id

    @observation_id.setter
    def observation_id(self, value):
        caom_util.type_check(value, str, 'observation_id', override=False)
        self._observation_id = value

    def get_uri(self):
        """A URI for this observation referenced in the caom system.

        This attribute is auto geneqrated from the other metadata.
        type: unicode string
        """
        return self._uri

    @property
    def planes(self):
        """A typed ordered dictionary containing plane objects associated with
        this observations.

        see caom2.Plane for details about a creating a plane.

        type: TypedOrderDict(Plane,)

        eg. Observation.planes.add(Plane("SCIENCE"))
        """
        return self._planes

    @planes.setter
    def planes(self, value):
        if value is None:
            self.planes = caom_util.TypedOrderedDict(Plane, )
        else:
            caom_util.type_check(value, caom_util.TypedOrderedDict, 'planes')
            self._planes = value

    @property
    def algorithm(self):
        """The process that was used to select this observation.

        normally 'exposure' for an observation or the name of a group
        process for a composite observation"""
        return self._algorithm

    @algorithm.setter
    def algorithm(self, value):
        if isinstance(value, str):
            value = Algorithm(value)
        caom_util.type_check(value, Algorithm, 'algorithm')
        self._algorithm = value

    @property
    def intent(self):
        """The original intent of having this data.

        type: ObservationIntentType

        see ObservationIntentType for allowed values

        """
        return self._intent

    @intent.setter
    def intent(self, value):
        if isinstance(value, str):
            value = ObservationIntentType(value)
        caom_util.type_check(value, ObservationIntentType, 'intent')
        self._intent = value

    @property
    def sequence_number(self):
        """An integer counter that reprsents this observation within a
        collection.  eg.  EXPNUM type: int
        """
        return self._sequence_number

    @sequence_number.setter
    def sequence_number(self, value):
        caom_util.type_check(value, int_32, 'sequence_number')
        self._sequence_number = int_32(value) if value is not None else None

    @property
    def type(self):
        """The OBSTYPE of the observation being recorded.

        eg. OBJECT, FLAT, BIAS
        type: unicode str
        """
        return self._type

    @type.setter
    def type(self, value):
        caom_util.type_check(value, str, 'type')
        self._type = value

    @property
    def proposal(self):
        """Refence to a Proposal object that describe the science proposal
        that lead to this observation.

        can be None
        see caom2.Proposal for help on building a Proposal object
        type: caom2.Proposal
        """
        return self._proposal

    @proposal.setter
    def proposal(self, value):
        caom_util.type_check(value, Proposal, "proposal")
        self._proposal = value

    @property
    def telescope(self):
        """Reference to a Telescope object associated with this observation.

        can be None
        type: caom2.Telescope
        """
        return self._telescope

    @telescope.setter
    def telescope(self, value):
        caom_util.type_check(value, Telescope, 'telescope')
        self._telescope = value

    @property
    def instrument(self):
        """Reference to an Instrument object associated with this observation.

        can be None
        type: caom2.Instrument
        """
        return self._instrument

    @instrument.setter
    def instrument(self, value):
        if isinstance(value, str):
            value = Instrument(value)
        caom_util.type_check(value, Instrument, "instrument")
        self._instrument = value

    @property
    def target(self):
        """Reference to a Target object associated with this observation.

        can be None
        type: caom2.Target
        """
        return self._target

    @target.setter
    def target(self, value):
        if isinstance(value, str):
            value = Target(value)
        caom_util.type_check(value, Target, 'target')
        self._target = value

    @property
    def environment(self):
        """Reference to an Environment object associated with this
        observation.

        can be None
        type: caom2.Environment
        """
        return self._environment

    @environment.setter
    def environment(self, value):
        caom_util.type_check(value, Environment, 'environment')
        self._environment = value

    @property
    def target_position(self):
        """Reference to a TargetPosition object associated
        with this observation.

        can be None
        type: caom2.TargetPosition
        """
        return self._target_position

    @target_position.setter
    def target_position(self, value):
        caom_util.type_check(value, TargetPosition, 'target_position')
        self._target_position = value

    @property
    def requirements(self):
        """Reference to a Requirements object associated
        with this observation.

        can be None
        type: caom2.Requirements
        """
        return self._requirements

    @requirements.setter
    def requirements(self, value):
        caom_util.type_check(value, Requirements, 'requirements')
        self._requirements = value

    @property
    def meta_release(self):
        """A datetime value indicating when the meta-data of this observation
        is publicly accessible.

        This only controls access to the information about the
        observation. Access to the observational data is controlled
        via the Plane.data_release attribute (see the planes
        attribute).

        eg. '2012/11/28 12:00:00'
        type: datatime
        """
        return self._meta_release

    @meta_release.setter
    def meta_release(self, value):
        caom_util.type_check(value, datetime, 'meta_release')
        self._meta_release = value

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


class Algorithm(CaomObject):
    """
    The concept of Algorithm is to provide a way for users to find all
    composite observation sets that have been built using a particular
    grouping algorithm (eg. the MegaPipe stacks).  For simple
    observations the algorithm is 'exposure'.
    """

    def __init__(self, name):
        """
        Initializes an Algorithm instance

        Arguments:
        name : name of the algorithm.  Should be 'exposure' if this is a
               simple observation, otherwise name of algorithm that selected
               composite members or just 'composite' works too.
        """
        caom_util.type_check(name, str, 'name', override=False)
        self._name = str(name)

    def _key(self):
        return self._name

    def __ne__(self, y):
        return not self.__eq__(y)

    # Properties

    @property
    def name(self):
        """
        algorithm that defines the composite grouping; 'exposure' for
        simple observations
        """
        return self._name


class SimpleObservation(Observation):
    """A convenience class for building Observations where
    algorithm='exposure'.

    Simple Observation : an object to store the metadata associated
    with a telescopic observation.

    Simple Observations are the basic package for the CAOM2
    model. Each simple observation refers to an observation by a
    telescope.

    The metadata content is stored in a variety of objects, see the
    caom2_observtion class for details.

    The algorithm value for a simple observation is forced set to
    "exposure"

    """
    _DEFAULT_ALGORITHM_NAME = 'exposure'

    def __init__(self,
                 collection,
                 observation_id,
                 algorithm=_DEFAULT_ALGORITHM_NAME,
                 sequence_number=None,
                 intent=None,
                 type=None,
                 proposal=None,
                 telescope=None,
                 instrument=None,
                 target=None,
                 meta_release=None,
                 meta_read_groups=None,
                 planes=None,
                 environment=None,
                 target_position=None
                 ):
        """
        collection - A name that describes a collection of data,
        nominally the name of a telescope

        observation_id - A UNIQUE identifier with in that collection
        """
        super(SimpleObservation, self).__init__(collection,
                                                observation_id,
                                                algorithm,
                                                sequence_number,
                                                intent,
                                                type,
                                                proposal,
                                                telescope,
                                                instrument,
                                                target,
                                                meta_release,
                                                meta_read_groups,
                                                planes,
                                                environment,
                                                target_position
                                                )

    @property
    def algorithm(self):
        """The algorithm that built the observation, for SimpleObservation
        this is always 'exposure'"""
        return super(SimpleObservation, self).algorithm

    @algorithm.setter
    def algorithm(self, value):
        # build an Algorithm type if passed a string...
        if value is None:
            raise ValueError('Algorithm name required')
        if isinstance(value, str):
            value = Algorithm(value)
        caom_util.type_check(value, Algorithm, 'algorithm', override=False)
        self._algorithm = value


class DerivedObservation(Observation):
    """
    Derived Observation

     A DerivedObservation is created by collecting data from multiple
    SimpleObservations together. They could represent:
        - stacked observations
        - observations that are extracted subsets of other observations
        - virtual observations that are define groups but don't have products

    """

    def __init__(self,
                 collection,
                 observation_id,
                 algorithm,
                 sequence_number=None,
                 intent=None,
                 type=None,
                 proposal=None,
                 telescope=None,
                 instrument=None,
                 target=None,
                 meta_release=None,
                 meta_read_groups=None,
                 planes=None,
                 environment=None,
                 target_position=None):
        super(DerivedObservation, self).__init__(
              collection=collection,
              observation_id=observation_id,
              algorithm=algorithm,
              sequence_number=sequence_number,
              intent=intent,
              type=type,
              proposal=proposal,
              telescope=telescope,
              instrument=instrument,
              target=target,
              meta_release=meta_release,
              meta_read_groups=meta_read_groups,
              planes=planes,
              environment=environment,
              target_position=target_position)
        self._members = caom_util.TypedSet(ObservationURI, )

    @property
    def algorithm(self):
        return super(DerivedObservation, self).algorithm

    @algorithm.setter
    def algorithm(self, value):
        if value is None:
            raise ValueError('Algorithm name required')
        if isinstance(value, str):
            value = Algorithm(value)
        if value.name == SimpleObservation._DEFAULT_ALGORITHM_NAME:
            raise ValueError("cannot set DerivedObservation.algorithm to {0}"
                             " (reserved for SimpleObservation)".format(value))
        self._algorithm = value

    @property
    def members(self):
        return self._members


@deprecated(version='CAOM2.4', reason='Replaced by DerivedObservation')
class CompositeObservation(DerivedObservation):
    """
    deprecated class
    """

    def __init__(self,
                 collection,
                 observation_id,
                 algorithm,
                 sequence_number=None,
                 intent=None,
                 type=None,
                 proposal=None,
                 telescope=None,
                 instrument=None,
                 target=None,
                 meta_release=None,
                 planes=None,
                 environment=None,
                 target_position=None):
        super(CompositeObservation, self).__init__(
            collection=collection,
            observation_id=observation_id,
            algorithm=algorithm,
            sequence_number=sequence_number,
            intent=intent,
            type=type,
            proposal=proposal,
            telescope=telescope,
            instrument=instrument,
            target=target,
            meta_release=meta_release,
            planes=planes,
            environment=environment,
            target_position=target_position)
        warnings.warn("CompositeObservation has been deprecated. Please use "
                      "DerivedObservation instead.")


class Environment(CaomObject):
    """A CAOM2 Environment Object.

    This object contains the various values that can be set in the environment
    table entry. Normally each Observation object will have an associate
    Environment Object."""

    def __init__(self):
        """
        Initializes an Environment instance

        Arguments:
        """
        self._seeing = None
        self._humidity = None
        self._elevation = None
        self._tau = None
        self._wavelength_tau = None
        self._ambient_temp = None
        self._photometric = None

    # Properties
    @property
    def seeing(self):
        """atmospheric seeing (FWHM) in arcsec.

        units: arcseconds
        """
        return self._seeing

    @seeing.setter
    def seeing(self, value):
        caom_util.type_check(value, float, 'seeing')
        caom_util.value_check(value, 0, 360 * 3600.0, 'seeing')
        self._seeing = value

    @property
    def humidity(self):
        """Relative humidity, expressed as a fraction between 0 and 200.

        units: fraction
        """
        return self._humidity

    @humidity.setter
    def humidity(self, value):
        # set the humidity to value, which  must be +ve fraction less than 200
        caom_util.type_check(value, float, 'humidity')
        caom_util.value_check(value, 0, 200, 'humidity')
        self._humidity = value

    @property
    def elevation(self):
        """ Elevation above horizon (0 to 90 deg) at which tau was measured.

        units: deg
        """
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        caom_util.type_check(value, float, 'elevation')
        caom_util.value_check(value, 0, 90, 'elevation')
        self._elevation = value

    @property
    def tau(self):
        """The tau at zennith at the time of observation.

        units:  fraction

        This tau can be used, in combination with the elevation,
        to determine the actual tau."""
        return self._tau

    @tau.setter
    def tau(self, value):
        caom_util.type_check(value, float, 'tau')
        # Value must be >= 0, but has no upper limit
        self._tau = value

    @property
    def wavelength_tau(self):
        """Wavelength at which tau was measured
        (normally 225GHz converted to wavelength).

        units: meters

        """
        return self._wavelength_tau

    @wavelength_tau.setter
    def wavelength_tau(self, value):
        caom_util.type_check(value, float, 'wavelength_tau')
        caom_util.value_check(value, 0, 1E3, 'wavelength_tau')
        self._wavelength_tau = value

    @property
    def ambient_temp(self):
        """The ambient air temperature at time of observation.

        unit: Celsius degrees
        """
        return self._ambient_temp

    @ambient_temp.setter
    def ambient_temp(self, value):
        caom_util.type_check(value, float, 'ambient_temp')
        self._ambient_temp = value

    @property
    def photometric(self):
        """A boolean flag (True/False) indicating if
        the observational conditions were photometric."""
        return self._photometric

    @photometric.setter
    def photometric(self, value):
        caom_util.type_check(value, bool, 'photometric')
        self._photometric = value


class Instrument(CaomObject):
    """The telescopic instrument that recorded a given observation.

    Each observation should have an associated instrument.

    eg:
    inst=caom2.Observation('CFHT','123456p','exposure').Instrument("MEGAPRIME")

    Each instrument can have list of keywords.
    eg.
    inst.keywords.append("shutter=closed")
    """

    def __init__(self, name):
        """
        Initializes a Instrument instance

        Arguments:
        name - name of the instrument
        """
        caom_util.type_check(name, str, 'name', override='none')
        self._name = name
        self._keywords = set()

    # Properties
    @property
    def name(self):
        """The name of the instrument.

        type: unicode string
        """
        return self._name

    @property
    def keywords(self):
        """A set of strings that are keywords associated with the instrument.

        eg.  keywords.append(("ccd=off","hot","shutter broken"))

        Keywords are stored as text in the database and are searched as text,
        so in the above example you can search for ccd=off but that is just
        matching the text 'ccd=off' and you can not do 'ccd!=on' to find that
        key/value pair.

        Also, the list is concated together and inserted into the
        model as a single text field, so keywords like 'shutter
        broken' will match 'shutter' and 'broken'.  If these keywords
        appear in a pick list 'shutter' and 'broken' will be seperate
        entries.  So shutter_broken is a better choice.
        """
        return self._keywords


class Proposal(CaomObject):
    """ Proposal """

    def __init__(self,
                 id,
                 pi_name=None,
                 project=None,
                 title=None):
        """
        Initializes a Proposal instance

        Arguments:
        myId : id of the proposal
        """

        self.id = id
        self.pi_name = pi_name
        self.project = project
        self.title = title

        self.keywords = set()

    # Properties

    @property
    def id(self):
        """The proposal ID.  Sometimes also called a RUNID.

        type: unicode string
        """
        return self._id

    @id.setter
    def id(self, value):
        caom_util.type_check(value, str, 'id')
        self._id = value

    @property
    def keywords(self):
        """A Set of keywords connected to this proposal.

        keywords are stored as a string of words and do not need to be
        key/value pairs.

        eg. Proposal.keywords.add('galaxies')

        type: set
        """
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        caom_util.type_check(value, set, 'keywords', override=False)
        self._keywords = value

    @property
    def pi_name(self):
        """The name (First Last) of the Principle Investigator of the
        Proposal.

        type: unicode string
        """
        return self._pi_name

    @pi_name.setter
    def pi_name(self, value):
        caom_util.type_check(value, str, 'pi_name')
        self._pi_name = value

    @property
    def project(self):
        """The name of a project associated with this proposal.

        type: unicode string
        """
        return self._project

    @project.setter
    def project(self, value):
        caom_util.type_check(value, str, 'project')
        self._project = value

    @property
    def title(self):
        """The title of the proposal.

        type: unicode string
        """
        return self._title

    @title.setter
    def title(self, value):
        caom_util.type_check(value, str, 'title')
        self._title = value


class Requirements(CaomObject):
    """ Requirements """

    def __init__(self, flag):
        """
        Construct an Requirements instance

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
        caom_util.type_check(value, Status, "flag")
        self._flag = value


class Target(CaomObject):
    """ Target """

    def __init__(self, name,
                 target_type=None,
                 standard=None,
                 redshift=None,
                 keywords=None,
                 moving=None,
                 target_id=None):
        """
        Initializes a Target instance

        Arguments:
        name : name of the target
        type : type of the target
        """

        self.name = name
        self.target_type = target_type
        self.standard = standard
        self.redshift = redshift
        if keywords is None:
            keywords = set()
        self.keywords = keywords
        self.moving = moving
        self.target_id = target_id

    # Properties

    @property
    def name(self):
        """A name for the target

        eg 'NGC 3115'
        type: unicode string

        """
        return self._name

    @name.setter
    def name(self, value):
        caom_util.type_check(value, str, "name", override=False)
        self._name = value

    @property
    def target_type(self):
        """A keyword describing the type of target.
        must be from the list
        """ + str(list(TargetType)) + """
        type: TargetType

        """
        return self._type

    @target_type.setter
    def target_type(self, value):
        if isinstance(value, str):
            value = TargetType(value)
        caom_util.type_check(value, TargetType, "target_type")
        self._type = value

    @property
    def keywords(self):
        """A set of keywords associated with this target.

        eg. keywords.add('galaxy')
        type: set

        """
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        caom_util.type_check(value, set, 'keywords', override=False)
        self._keywords = value

    @property
    def standard(self):
        """Is this a standard field?

        eg True
        type: bool

        """
        return self._standard

    @standard.setter
    def standard(self, value):
        caom_util.type_check(value, bool, 'standard')
        self._standard = value

    @property
    def redshift(self):
        """The redshift of the observed target.

        eg 1.2  (can be None)
        type: float

        """
        return self._redshift

    @redshift.setter
    def redshift(self, value):
        caom_util.type_check(value, float, 'redshift')
        caom_util.value_check(value, -0.5, 1200, 'redshift')
        self._redshift = value

    @property
    def moving(self):
        """Is this a moving target?

        eg True
        type: bool

        """
        return self._moving

    @moving.setter
    def moving(self, value):
        caom_util.type_check(value, bool, 'moving')
        self._moving = value

    @property
    def target_id(self):
        return self._target_id

    @target_id.setter
    def target_id(self, value):
        caom_util.type_check(value, str, 'target_id')
        if value is not None:
            tmp = urlsplit(value)
            if tmp.geturl() != value:
                raise ValueError("Invalid URI: " + value)
        self._target_id = value


class TargetPosition(CaomObject):
    """ TargetPosition """

    def __init__(self, coordinates, coordsys, equinox=None):
        """
        Initialize a TargetPosition instance.

        Arguments:
        coordinates : target position as a Point.
        coordsys: target coordsys
        equinox: target equinox
        """
        self.coordinates = coordinates
        self.coordsys = coordsys
        self.equinox = equinox

    # Properties

    @property
    def coordinates(self):
        """ Coordinates """
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        caom_util.type_check(value, Point, "coordinates")
        self._coordinates = value

    @property
    def coordsys(self):
        """ Coordsys """
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        caom_util.type_check(value, str, "coordsys")
        self._coordsys = value

    @property
    def equinox(self):
        """ Equinox """
        return self._equinox

    @equinox.setter
    def equinox(self, value):
        caom_util.type_check(value, float, "equinox")
        self._equinox = value


class Telescope(CaomObject):
    """ Telescope """

    def __init__(self, name,
                 geo_location_x=None,
                 geo_location_y=None,
                 geo_location_z=None,
                 keywords=None
                 ):
        """
        Initializes a Telescope instance

        Arguments:
        name : name of the telescope
        """
        if not name:
            raise AttributeError("Telescope name required")
        caom_util.type_check(name, str, 'name', override=False)
        self._name = name
        self.geo_location_x = geo_location_x
        self.geo_location_y = geo_location_y
        self.geo_location_z = geo_location_z
        if keywords is None:
            keywords = set()
        self.keywords = keywords

    # Properties

    @property
    def name(self):
        """a name for this facility.

        eg. CFHT
        type: unicode string

        """
        return self._name

    @property
    def keywords(self):
        """A set that contains keywords associated with this telescope

        eg.  keywords.add('big')
        type: set

        """
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        caom_util.type_check(value, set, 'keywords', override=False)
        self._keywords = value

    @property
    def geo_location_x(self):
        """The x geocentric location of the telescope.

        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_x

    @geo_location_x.setter
    def geo_location_x(self, value):
        caom_util.type_check(value, float, 'geo_location_x')
        self._geo_location_x = value

    @property
    def geo_location_y(self):
        """the y geocentric (ECR) location of the telescope.

        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_y

    @geo_location_y.setter
    def geo_location_y(self, value):
        caom_util.type_check(value, float, 'geo_location_y')
        self._geo_location_y = value

    @property
    def geo_location_z(self):
        """the z geocentric (ECR) location of the telescope.
        This should be valid at time of MJD-OBS.

        These coordinates should be in the ITRS reference,
        basically whatever a GPS device is saying.

        The directions of the x/y/z follow the
        Earth Centred Rotation frame.

        units: m
        type: float

        """
        return self._geo_location_z

    @geo_location_z.setter
    def geo_location_z(self, value):
        caom_util.type_check(value, float, 'geo_location_z')
        self._geo_location_z = value
