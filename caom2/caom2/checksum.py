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

import hashlib
import logging
import struct
import uuid
from datetime import datetime
import argparse
import sys
from . import obs_reader_writer
import warnings
from builtins import bytes, int, str

from caom2.caom_util import TypedSet, TypedList, TypedOrderedDict, int_32
from caom2.common import CaomObject, AbstractCaomEntity, ObservationURI
from caom2.common import ChecksumURI
from caom2.observation import Observation
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aenum import Enum

__all__ = ['get_meta_checksum', 'get_acc_meta_checksum',
           'update_meta_checksum']

"""
This modules contains the functionality for calculating the checksums
corresponding to entites in CAOM2. Two types of checksums are available:
    get_meta_checksum - returns the checksum of the entity attributes except
                        its CAOM2 children
    get_acc_meta_checksum - returns the accumulated checksum of the entity
                        attributes including those of the CAOM2 children
    update_meta_checksum - updates the meta_checksum and acc_meta_checksum of
                        all the elements of an observation

The checksums can be used to get the state of the entity: the metachecksum
represents the state the entity itself while the accumulated metachecksum
represents the state of the entity and all the children below. These states
can be saved and later compared to quickly detect state changes.

IMPORTANT NODE: The checksum algorithms use introspection to automatically
find the attributes that are part of the model. It is therefore very important
that attributes that are not part of the CAOM model
(http://www.opencadc.org/caom2) be prefixed with an '_' so that the checksum
algorithms ignore them in their computations. Equally important
is to use the names of the attributes from the model (with the Python syntax)
as the algorithms parses them in alphabetical order.

Gotchas to look for when calculating the checksum:
- There are two types of integers in the model: int (32 bit) and long (64 bit).
In this implementation the corresponding types are int_32 and int.
- Sets are ordered alphabetical. Therefore their members have to implement,
at the minimum, __eq__ and __lt__ that will result in proper sorting

"""
logger = logging.getLogger('checksum')
logging.basicConfig()


def get_meta_checksum(entity):
    """
    Uses md5sum algorithm to calculate the checksum of the caom2_object
    excluding its CAOM2 children.
    calculation order:
    1. id for entities
    2.state fields in alphabetic order; depth - first recursion so foo.abc.x
    comes before foo.def

    value handling:
    - Date: truncate time to whole number of seconds and treat as a long
    - String: UTF - 8 encoded bytes
    - float: IEEE754 double(8 bytes)
    - boolean: convert to single byte, false = 0, true = 1(1bytes)
    - byte: as- is (1 byte)
    - integer: 9(8 bytes, network byte order == big endian)
    - unrecognized classes: encode their string representation

    :param entity: CAOM2 entity
    :return: md5 checksum corresponding to the entity metadata
    """
    if not isinstance(entity, AbstractCaomEntity):
        raise AttributeError('AbstractCaomEntity expected')
    md5 = hashlib.md5()
    update_caom_checksum(md5, entity)
    return ChecksumURI('md5:{}'.format(md5.hexdigest()))


def get_acc_meta_checksum(entity, no_logging=False):
    """
    Similar to get_meta_checksum except that the accumulated checksum of
    the CAOM2 children are also included in alphabetical order of their ids

    :param entity: CAOM2 entity
    :param no_logging: if True turns off any logging while running this method
    :return: md5 checksum corresponding to the entity metadata
    """
    if not isinstance(entity, AbstractCaomEntity):
        raise AttributeError("AbstractCaomEntity class expected")
    if no_logging:
        log_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
    md5 = hashlib.md5()
    update_acc_checksum(md5, entity)
    if no_logging:
        logger.setLevel(log_level)
    return ChecksumURI('md5:{}'.format(md5.hexdigest()))


def update_meta_checksum(obs):
    """
    Convenience function that updates the meta_checksum and acc_meta_checksum
    of an observation and all it's entities
    :param obs: observation to be updated
    """
    if not isinstance(obs, Observation):
        raise AttributeError('Observation required')
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    logger.debug('*** START CHUNK ***')
                    chunk.meta_checksum = get_meta_checksum(chunk)
                    chunk.acc_meta_checksum = \
                        get_acc_meta_checksum(chunk, no_logging=True)
                    logger.debug('*** END CHUNK ***')
                logger.debug('*** START PART ***')
                part.meta_checksum = get_meta_checksum(part)
                part.acc_meta_checksum = get_acc_meta_checksum(part,
                                                               no_logging=True)
                logger.debug('*** END PART ***')
            logger.debug('*** START ARTIFACT ***')
            artifact.meta_checksum = get_meta_checksum(artifact)
            artifact.acc_meta_checksum = get_acc_meta_checksum(artifact,
                                                               no_logging=True)
            logger.debug('*** END ARTIFACT ***')
        logger.debug('*** START PLANE ***')
        plane.meta_checksum = get_meta_checksum(plane)
        plane.acc_meta_checksum = get_acc_meta_checksum(plane, no_logging=True)
        logger.debug('*** END PLANE ***')
    logger.debug('*** START OBSERVATION ***')
    obs.meta_checksum = get_meta_checksum(obs)
    obs.acc_meta_checksum = get_acc_meta_checksum(obs, no_logging=True)
    logger.debug('*** END OBSERVATION ***')


def update_acc_checksum(checksum, entity):
    """
    Updates the checksum alogrithm with the bytes corresponding to the entity.
    It first generates the corresponding bytes for the meta checksum of the
    entity. To that, it then adds the bytes corresponding to the accumulated
    meta checksum of each child listed in alphabeticals order of their id.

    :param checksum: checksum algorithm that consumes the bytes (md5)
    :param entity: entity to generate the bytes for and consume them
    """
    if not isinstance(entity, AbstractCaomEntity):
        raise AttributeError('AbstractCaomEntity expected')
    md5 = hashlib.md5()
    update_caom_checksum(md5, entity)
    checksum.update(md5.digest())
    # go through the children and calculate the acc checksums
    children = {}
    for i in dir(entity):
        if not callable(getattr(entity, i)) and not i.startswith('_'):
            attrib = getattr(entity, i)
            if (isinstance(attrib, TypedOrderedDict) or
                    isinstance(attrib, TypedList) or
                    isinstance(attrib, TypedSet)) and (
                    issubclass(attrib.key_type, AbstractCaomEntity)):
                if (isinstance(attrib, TypedOrderedDict)):
                    values = attrib.values()
                else:
                    values = attrib
                for j in values:
                    children[j._id] = j
    for i in sorted(children):
        md5 = hashlib.md5()
        update_acc_checksum(md5, children[i])
        checksum.update(md5.digest())


def update_checksum(checksum, value, attribute=''):
    """
    Updates the checksum algorithm (md5) of (mostly) native types with
    corresponding bytes (network byte order == big endian)
    :param checksum: the checksum algorithm to update (md5)
    :param value: value to translate into bytes in order to update the
    checksum algorithm
    :param attribute: name of the attribute this value belongs to
    (used for debugging only)
    """

    if type(value) is None:
        logger.debug('Empty attribute {}'.format(attribute))
        return

    b = None

    if isinstance(value, ObservationURI) or isinstance(value, ChecksumURI):
        b = value.uri.encode('utf-8')
    elif isinstance(value, CaomObject):
        logger.debug('Process object {}'.format(attribute))
        update_caom_checksum(checksum, value, attribute)
    elif isinstance(value, bytes):
        b = value
    elif isinstance(value, bool):
        b = struct.pack('!?', value)
    elif isinstance(value, float):
        b = struct.pack('!d', value)
    elif isinstance(value, int_32):
        # must be before int
        b = struct.pack('!l', value)
    elif isinstance(value, int):
        b = struct.pack('!q', value)
    elif isinstance(value, str):
        b = value.strip().encode('utf-8')
    elif isinstance(value, datetime):
        b = struct.pack('!q', int(
            (value - datetime(1970, 1, 1)).total_seconds()))
    elif isinstance(value, set) or \
            (isinstance(value, TypedSet) and not
                isinstance(value.key_type, AbstractCaomEntity)):
        for i in sorted(value):
            update_checksum(checksum, i, attribute)
    elif isinstance(value, list) or isinstance(value, TypedList):
        for i in value:
            if not isinstance(i, AbstractCaomEntity):
                update_checksum(checksum, i, attribute)
    elif isinstance(value, Enum):
        update_checksum(checksum, value.value, attribute)
    elif isinstance(value, uuid.UUID):
        b = value.bytes
    elif isinstance(value, TypedOrderedDict):
        # calculate the checksum of each component and add them in
        # alphabetical order of their ids
        # Note: ignore dictionaries of AbstractCaomEntity types
        checksums = []
        for i in value:
            if not isinstance(value[i], AbstractCaomEntity):
                checksums.append(value[i]._id)
        for i in sorted(checksums):
            update_checksum(checksum, checksum[i], attribute)
    else:
        raise ValueError(
            'Cannot transform in bytes: {}({})'.format(value, type(value)))

    if b is not None:
        checksum.update(b)
        if logger.isEnabledFor(logging.DEBUG):
            md5 = hashlib.md5()
            md5.update(b)
            logger.debug('Encoded attribute ({}) {} = {} -- {}'.
                         format(type(value), attribute,
                                value, md5.hexdigest()))


def update_caom_checksum(checksum, entity, parent=None):
    """
    Method to go through the attributes of a CAOM class and update the
    checksum of each of them.
    Uses introspection and goes through the attributes in alphabetical order.
    :param checksum: checksum algorithm to update (md5)
    :param entity: entity to go through
    :param parent: parent of the entity (used for debugging only)
    """
    if not isinstance(entity, CaomObject):
        raise AttributeError('CaomObject expected')
    # get the id first
    if isinstance(entity, AbstractCaomEntity):
        update_checksum(checksum, entity._id)
        if entity._meta_producer:
            update_checksum(checksum, entity._meta_producer)

    # determine the excluded fields if necessary
    checksum_excluded_fields = []
    if not hasattr(update_caom_checksum, "checksum_excluded_fields"):
        # this is a way to do it statically
        checksum_excluded_fields = [i for i in dir(AbstractCaomEntity)
                                    if not callable(
                getattr(AbstractCaomEntity, i))
                                    and not i.startswith('_')]
        # get the fields in alphabetical order, depth first recursion
        # but remove the
    for i in sorted(dir(entity)):
        if not callable(getattr(entity, i)) and not i.startswith('_') and \
                        i not in checksum_excluded_fields:
            if getattr(entity, i) is not None:
                atrib = '{}.{}'.format(parent, i) if parent is not None else i
                update_checksum(checksum, getattr(entity, i), atrib)


def checksum_diff():
    """
    Calculates the checksum of elements in an observation and compares it
    with the ones specified in that observation
    """

    parser = argparse.ArgumentParser(
        description='Compare observation checksum')
    parser.add_argument('file', help='Observation file')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Display details')
    parser.add_argument('-o', '--output', required=False,
                        help='save checked file in this file')

    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("caom2_checksum: error: too few arguments")
        sys.exit(-1)

    if args.debug:
        logger.setLevel(logging.DEBUG)

    reader = obs_reader_writer.ObservationReader(True)
    orig = reader.read(args.file)

    # read again for the observation that would have the checksums updated
    actual = reader.read(args.file)

    update_meta_checksum(actual)

    print('** metaChecksum **\n')
    mistmatches = 0
    for plane in zip(orig.planes.values(), actual.planes.values()):
        for artifact in zip(plane[0].artifacts.values(),
                            plane[1].artifacts.values()):
            for part in zip(artifact[0].parts.values(),
                            artifact[1].parts.values()):
                for chunk in zip(part[0].chunks, part[1].chunks):
                    mistmatches += _print_diff(chunk[0], chunk[1])
                mistmatches += _print_diff(part[0], part[1])
            mistmatches += _print_diff(artifact[0], artifact[1])
        mistmatches += _print_diff(plane[0], plane[1])
    mistmatches += _print_diff(orig, actual)

    if args.output:
        writer = obs_reader_writer.ObservationWriter(validate=True)
        writer.write(actual, args.output)

    print("Total: {} mistmatches".format(mistmatches))
    if mistmatches > 0:
        sys.exit(-1)


def _print_diff(orig, actual):
    elem_type = str(type(orig)).split('.')[1]
    mistmatches = 0
    if orig.meta_checksum == actual.meta_checksum:
        print('{}: {} {} == {}'.format(elem_type, orig._id,
                                       orig.meta_checksum.checksum,
                                       actual.meta_checksum.checksum))
    else:
        print('{}: {} {} != {} [MISMATCH]'.
              format(elem_type, orig._id, orig.meta_checksum.checksum,
                     actual.meta_checksum.checksum))
        mistmatches += 1

    if elem_type != 'chunk':
        # do the accummulated checksums
        if orig.acc_meta_checksum == actual.acc_meta_checksum:
            print('{}: {} {} == {}'.
                  format(elem_type, orig._id, orig.acc_meta_checksum.checksum,
                         actual.acc_meta_checksum.checksum))
        else:
            print('{}: {} {} != {} [MISMATCH]'.
                  format(elem_type, orig._id, orig.acc_meta_checksum.checksum,
                         actual.acc_meta_checksum.checksum))
            mistmatches += 1
    return mistmatches
