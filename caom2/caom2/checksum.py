# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2017.                            (c) 2017.
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

import hashlib
import logging
import struct
import uuid
from datetime import datetime

from aenum import Enum
from builtins import bytes, int, str

from caom2.caom_util import TypedSet, TypedList, TypedOrderedDict, int_32
from caom2.common import CaomObject, AbstractCaomEntity, ObservationURI, \
    ChecksumURI

__all__ = ['get_meta_checksum', 'get_acc_meta_checksum']

"""
This modules contains the functionality for calculating the checksums
corresponding to entites in CAOM2. Two types of checksums are available:
    get_meta_checksum - returns the checksum of the entity attributes except
                        its CAOM2 children
    get_acc_meta_checksum - returns the accumulated checksum of the entity
                        attributes including those of the CAOM2 children

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
    assert (isinstance(entity, AbstractCaomEntity))
    md5 = hashlib.md5()
    update_caom_checksum(md5, entity)
    return ChecksumURI('md5:{}'.format(md5.hexdigest()))


def get_acc_meta_checksum(entity):
    """
    Similar to get_meta_checksum except that the accumulated checksum of
    the CAOM2 children are also included in alphabetical order of their ids

    :param entity: CAOM2 entity
    :return: md5 checksum corresponding to the entity metadata
    """
    assert (isinstance(entity, AbstractCaomEntity))
    md5 = hashlib.md5()
    update_acc_checksum(md5, entity)
    return ChecksumURI('md5:{}'.format(md5.hexdigest()))


def update_acc_checksum(checksum, entity):
    """
    Updates the checksum alogrithm with the bytes corresponding to the entity.
    It first generates the corresponding bytes for the meta checksum of the
    entity. To that, it then adds the bytes corresponding to the accumulated
    meta checksum of each child listed in alphabeticals order of their id.

    :param checksum: checksum algorithm that consumes the bytes (md5)
    :param entity: entity to generate the bytes for and consume them
    """
    assert (isinstance(entity, AbstractCaomEntity))
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
        logger.debug('Encoded empty attribute {}'.format(attribute))
        return

    if isinstance(value, ObservationURI) or isinstance(value, ChecksumURI):
        logger.debug('Encoded attribute uri {} = {}'.format(attribute, value))
        checksum.update(value.uri.encode('utf-8'))
    elif isinstance(value, CaomObject):
        logger.debug('Encoded attribute {}'.format(attribute))
        update_caom_checksum(checksum, value, attribute)
    elif isinstance(value, bytes):
        logger.debug(
            'Encoded attribute bytes {} = {}'.format(attribute, value))
        checksum.update(value)
    elif isinstance(value, bool):
        logger.debug('Encoded attribute bool {} = {}'.format(attribute, value))
        checksum.update(struct.pack('!?', value))
        # elif isinstance(value, float_32):
        # must be before float
        # checksum.update(struct.pack('!f', value))
    elif isinstance(value, float):
        logger.debug(
            'Encoded attribute float {} = {}'.format(attribute, value))
        checksum.update(struct.pack('!d', value))
    elif isinstance(value, int_32):
        # must be before int
        logger.debug(
            'Encoded attribute int_32 {} = {}'.format(attribute, value))
        checksum.update(struct.pack('!l', value))
    elif isinstance(value, int):
        logger.debug('Encoded attribute int {} = {}'.format(attribute, value))
        checksum.update(struct.pack('!q', value))
    elif isinstance(value, str):
        logger.debug('Encoded attribute str {} = {}'.format(attribute, value))
        checksum.update(value.encode('utf-8'))
    elif isinstance(value, datetime):
        logger.debug(
            'Encoded attribute datetime {} = {}'.format(attribute, value))
        checksum.update(struct.pack('!q', int(
            (value - datetime(1970, 1, 1)).total_seconds())))
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
        logger.debug('Encoded attribute uuid {} = {}'.format(attribute, value))
        checksum.update(value.bytes)
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


def update_caom_checksum(checksum, entity, parent=None):
    """
    Method to go through the attributes of a CAOM class and update the
    checksum of each of them.
    Uses introspection and goes through the attributes in alphabetical order.
    :param checksum: checksum algorithm to update (md5)
    :param entity: entity to go through
    :param parent: parent of the entity (used for debugging only)
    """
    assert isinstance(entity, CaomObject)
    # get the id first
    if isinstance(entity, AbstractCaomEntity):
        update_checksum(checksum, entity._id)

    # determine the excluded fields if necessary
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
