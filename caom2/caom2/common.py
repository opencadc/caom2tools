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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import uuid
from datetime import datetime

from builtins import int, str
from six.moves.urllib.parse import SplitResult, urlsplit

from . import caom_util

__all__ = ['CaomObject', 'AbstractCaomEntity', 'ObservationURI', 'ChecksumURI']


def get_current_ivoa_time():
    """Generate a datetime with 3 digit microsecond precision.

    return: datatime
        IVOA date format to millisecond precision.
    """
    now = datetime.now()
    return datetime(now.year, now.month, now.day, now.hour, now.minute,
                    now.second, int(str(now.microsecond)[:-3] + '000'))


class CaomObject(object):
    """
    setup all objects with the same generic equality, str and repr methods
    """

    def __init__(self):
        pass

    def __str__(self):
        args = inspect.getargspec(self.__init__).args[1:]
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
        args = inspect.getargspec(self.__init__).args[1:]
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

    def __init__(self, fulluuid=False):
        super(CaomObject, self).__init__()
        self._id = AbstractCaomEntity._gen_id(fulluuid)
        self.last_modified = None
        self.max_last_modified = None
        self.meta_checksum = None
        self.acc_meta_checksum = None

    @classmethod
    def _gen_id(cls, fulluuid=False):
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
        tmp = urlsplit(uri)

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
                'Canot compare ObservationURI with {}'.format(type(other)))
        return self.uri < other.uri

    def __eq__(self, other):
        if not isinstance(other, ObservationURI):
            raise ValueError(
                'Canot compare ObservationURI with {}'.format(type(other)))
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
        tmp = urlsplit(uri)

        algorithm = tmp.scheme
        checksum = tmp.path

        if algorithm is None:
            raise ValueError(
                ("A checksum scheme noting the algorithm is "
                 "required.. received: {}")
                .format(uri))

        if checksum is None:
            raise ValueError(
                "checksum uri did not contain an checksum part. received: {}"
                .format(uri))
        caom_util.validate_path_component(self, "checksum", checksum)

        (self._uri, self._algorithm, self._checksum) = (
            tmp.geturl(), algorithm, checksum)
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
