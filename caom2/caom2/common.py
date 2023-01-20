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

import inspect
import uuid
from datetime import datetime

from builtins import int, str
from urllib.parse import SplitResult, urlparse, urlsplit
import logging

from . import caom_util
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aenum import Enum


__all__ = ['CaomObject', 'AbstractCaomEntity', 'ObservationURI', 'ChecksumURI',
           'VocabularyTerm']

_OBSCORE_VOCAB_NS = "http://www.ivoa.net/std/ObsCore"
_CAOM_VOCAB_NS = "http://www.opencadc.org/caom2/DataProductType"

logger = logging.getLogger('caom2')


def get_current_ivoa_time():
    """Generate a datetime with 3 digit microsecond precision.

    return: datatime
        IVOA date format to millisecond precision.
    """
    now = datetime.now()
    return datetime(now.year, now.month, now.day, now.hour, now.minute,
                    now.second, int(str(now.microsecond)[:-3] + '000'))


class OrderedEnum(Enum):
    """
    Enums are in the order of their definition.

    TODO: not sure this is required in Python 3
    enum.Enum is supposed to support this.
    """

    def __init__(self, *args):
        super(Enum, self).__init__()
        self._order = len(self.__class__.__members__) + 1

    def __ge__(self, other):
        if self.__class__ is other.__class__:
            return self._order >= other._order
        return NotImplemented

    def __gt__(self, other):
        if self.__class__ is other.__class__:
            return self._order > other._order
        return NotImplemented

    def __le__(self, other):
        if self.__class__ is other.__class__:
            return self._order <= other._order
        return NotImplemented

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self._order < other._order
        return NotImplemented


class CaomObject(object):
    """
    setup all objects with the same generic equality, str and repr methods
    """

    def __init__(self):
        pass

    def __str__(self):
        args = inspect.getfullargspec(self.__init__).args[1:]
        class_name = self.__class__.__name__
        return "\n".join(["{}.{} : {}".
                         format(class_name, arg, getattr(self, arg, None))
                          for arg in args])

    def __eq__(self, other):
        if type(other) == type(self):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __repr__(self):
        args = inspect.getfullargspec(self.__init__).args[1:]
        class_name = ""
        if self.__class__.__module__ != '__main__':
            class_name += self.__class__.__module__ + "."
        class_name += self.__class__.__name__
        pading = " " * (len(class_name) + 1)
        return class_name + "(" + (
            ",\n" + pading).join(
            ["%s=%r" % (arg, getattr(self, arg, None)
                        ) for arg in args]) + ")"


class AbstractCaomEntity(CaomObject):
    """Class that defines the persistence unique ID and last mod date """

    def __init__(self, fulluuid=True, meta_producer=None):
        super(CaomObject, self).__init__()
        self._id = AbstractCaomEntity._gen_id(fulluuid)
        self.last_modified = None
        self.max_last_modified = None
        self.meta_checksum = None
        self.acc_meta_checksum = None
        self._meta_producer = meta_producer

    @classmethod
    def _gen_id(cls, fulluuid=True):
        """Generate a 128 but UUID by default. For backwards compatibility
        allow creation of a 64 bit UUID using a rand number for the
        lower 64 bits. First two bytes of the random number are generated
        with the random and the last 6 bytes from the current time
        in microseconds.

        return: UUID
        """

        gen_id = uuid.uuid4()
        if fulluuid:
            return gen_id
        else:
            return uuid.UUID(fields=(0x00000000, 0x0000, 0x0000,
                                     gen_id.clock_seq_hi_variant,
                                     gen_id.clock_seq_low, gen_id.node))

    def compute_meta_checksum(self):
        raise NotImplementedError(
            "meta checksum calculation not yet implemented.")

    @property
    def last_modified(self):
        return self._last_modified

    @last_modified.setter
    def last_modified(self, value):
        if value is None:
            self._last_modified = None
        else:
            caom_util.type_check(value, datetime, "last_modified", False)
            self._last_modified = value

    @property
    def max_last_modified(self):
        return self._max_last_modified

    @max_last_modified.setter
    def max_last_modified(self, value):
        if value is None:
            self._max_last_modified = None
        else:
            caom_util.type_check(value, datetime, "max_last_modified", False)
            self._max_last_modified = value

    @property
    def meta_checksum(self):
        """the meta checksum value

        type: ChecksumURI

        """
        return self._meta_checksum

    @meta_checksum.setter
    def meta_checksum(self, value):
        if value is None:
            self._meta_checksum = None
        else:
            caom_util.type_check(value, ChecksumURI, "meta_checksum", False)
            self._meta_checksum = value

    @property
    def acc_meta_checksum(self):
        """the accumulated meta checksum value

        type: ChecksumURI

        """
        return self._acc_meta_checksum

    @acc_meta_checksum.setter
    def acc_meta_checksum(self, value):
        if value is None:
            self._acc_meta_checksum = None
        else:
            caom_util.type_check(value, ChecksumURI, "acc_meta_checksum",
                                 False)
            self._acc_meta_checksum = value

    @property
    def meta_producer(self):
        """
        Returns meta producer
        type: URI

        """
        return self._meta_producer

    @meta_producer.setter
    def meta_producer(self, value):
        if value is None:
            self._meta_producer = None
        else:
            try:
                urlparse(value)
            except ValueError:
                raise TypeError('Expected any IVOA URI for meta_producer, '
                                'received {}'.format(value))
            self._meta_producer = value


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
        if not tmp.geturl() == value:
            raise AttributeError("Invalid URI: " + value)
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


class ObservationURI(CaomObject):
    """ Observation URI """

    _SCHEME = str("caom")

    def __init__(self, uri):
        """
        Initializes an Observation instance

        Arguments:
        uri : URI corresponding to observation
        """
        super(CaomObject, self).__init__()
        tmp = urlparse(uri)

        if tmp.scheme != ObservationURI._SCHEME:
            raise ValueError(
                "uri must be have scheme of {}. received: {}"
                .format(ObservationURI._SCHEME, uri))
        if tmp.geturl() != uri:
            raise ValueError(
                "uri parsing failure.  received: {}".format(uri))

        self._uri = tmp.geturl()
        (collection, observation_id) = tmp.path.split("/")
        if collection is None:
            raise ValueError(
                "uri did not contain a collection part. received: {}"
                .format(uri))
        caom_util.validate_path_component(self, "collection", collection)
        if observation_id is None:
            raise ValueError(
                "uri did not contain an observation_id part. received: {}"
                .format(uri))
        caom_util.validate_path_component(self, "observation_id",
                                          observation_id)
        (self._collection, self._observation_id) = (collection, observation_id)
        self._print_attributes = ['uri', 'collection', 'observation_id']

    def _key(self):
        return self.uri

    def __hash__(self):
        return hash(self._key())

    def __lt__(self, other):
        if not isinstance(other, ObservationURI):
            raise ValueError(
                'Cannot compare ObservationURI with {}'.format(type(other)))
        return self.uri < other.uri

    def __eq__(self, other):
        if not isinstance(other, ObservationURI):
            raise ValueError(
                'Cannot compare ObservationURI with {}'.format(type(other)))
        return self.uri == other.uri

    @classmethod
    def get_observation_uri(cls, collection, observation_id):
        """
        Initializes an Observation URI instance

        Arguments:
        collection : collection
        observation_id : ID of the observation
        """

        caom_util.type_check(collection, str, "collection", override=False)
        caom_util.type_check(observation_id, str, "observation_id",
                             override=False)

        caom_util.validate_path_component(cls, "collection", collection)
        caom_util.validate_path_component(cls, "observation_id",
                                          observation_id)

        uri = SplitResult(ObservationURI._SCHEME, "",
                          collection + "/" + observation_id,
                          "", "").geturl()
        return cls(uri)

    # Properties

    @property
    @classmethod
    def scheme(cls):
        """The scheme defines where this Observation can be looked up.

        Only 'caom' is currently supported."""
        return cls._SCHEME

    @property
    def uri(self):
        """The uri that the caom service can use to find the observation"""
        return self._uri

    @property
    def collection(self):
        """The collection part of this Observations uri"""
        return self._collection

    @property
    def observation_id(self):
        """The observation_id of this Observations uri"""
        return self._observation_id


class ChecksumURI(CaomObject):
    """ Checksum URI """

    def __init__(self, uri):
        """
        Initializes an Checksum URI instance

        Arguments:
        uri : Checksum URI in the format Algorithm:ChecksumValue
        """
        super(CaomObject, self).__init__()
        # note: urlparse does not recognize scheme in uri of form scheme:val
        tmp = uri.split(':', 1)

        # TODO change this raise a ValueError when the rule is being enforced
        if len(tmp) < 2:
            logger.warning(("A checksum scheme noting the algorithm is "
                            "required.. received: {}").format(uri))
            algorithm = None
            checksum = tmp[0]
        else:
            algorithm = tmp[0]
            checksum = tmp[1]

        if checksum is None:
            raise ValueError(
                "checksum uri did not contain an checksum part. received: {}"
                .format(uri))
        caom_util.validate_path_component(self, "checksum", checksum)

        (self._uri, self._algorithm, self._checksum) = (
            uri, algorithm, checksum)
        self._print_attributes = ['uri', 'algorithm', 'checksum']

    def _key(self):
        return self.uri

    def __eq__(self, y):
        if isinstance(y, ChecksumURI):
            return self._key() == y._key()
        return False

    def __hash__(self):
        return hash(self._key())

    def get_bytes(self):
        return bytearray.fromhex(self._checksum)

    # Properties
    @property
    def uri(self):
        """The uri that the caom service can use to find the observation"""
        return self._uri

    @property
    def algorithm(self):
        """The checksum algorithm"""
        return self._algorithm

    @property
    def checksum(self):
        """The checksum value"""
        return self._checksum
