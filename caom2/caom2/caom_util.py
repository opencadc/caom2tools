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

"""
A collection of utility methods for CAOM2.

This module contains a number of 'TYPED' objects for use in CAOM2.
These Typed versions were implemented so that content checking near
the first point of use could be implemented.  This helps the data
engineer get the correct meta data more quickly.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import math
import struct
import sys
import uuid
from datetime import datetime

import six
from builtins import bytes, int

from caom2.common import AbstractCaomEntity


__all__ = ['TypedList', 'TypedSet', 'TypedOrderedDict', 'ClassProperty',
           'get_differences']

# TODO both these are very bad, implement more sensibly
IVOA_DATE_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"


class int_32(int):
    """
    The checksum algorithm must distinguished between 32 bit integers and 64
    bit integers. This subtype of int is used to tell the algorithm to use
    only 4 bytes in the checksum of an attribute of this type.
    """

    def __new__(cls, *args, **kwargs):
        return int.__new__(cls, *args, **kwargs)


def validate_path_component(caller, name, test):
    """
    Function to validate a URI path component. Component is invalid
    if it contains space ( ), slash (/), escape (\\) or percent (%) characters.

    Arguments:
    caller : caller object
    name : name of the component
    test : component to be tested

    An assertionError is thrown when the the provided test argument
    is invalid
    """

    assert (' ' not in test and
            '/' not in test and
            '||' not in test and
            '%' not in test), (
        caller.__class__.__name__ + ": invalid " + name +
        ": may not contain space ( ), slash (/), escape (\\), or percent (%)")


def date2ivoa(d):
    """
    Takes a datetime and returns a string formatted
    to the IVOA date format yyyy-MM-dd'T'HH:mm:ss.SSS
    """

    if d is None:
        return None
    return d.strftime(IVOA_DATE_FORMAT)[:23]


def str2ivoa(s):
    """Takes a IVOA date formatted string and returns a datetime"""

    if s is None:
        return None
    return datetime.strptime(s, IVOA_DATE_FORMAT)


def attr2str(s):
    pass


def repr2str(s):
    pass


def uuid2long(uid):
    """
    UUID is 128 bits (32 bytes). Unpack the 32 bytes into two
    16 byte longs. For CAOM-2.0 compatibility only the least significant
    16 bytes in the UUID should have a value.

    return the UUID least significant bytes as a long.
    """
    longs = struct.unpack(str('>qq'), bytes(uid.bytes))
    if longs[0] != 0:
        raise ValueError("lossy conversion from UUID to long: {}".format(uid))
    return longs[1]


def long2uuid(lng):
    """
    Takes a long and creates a UUID using the 16 byte long
    as the least significant bytes in the 32 byte UUID.
    """
    if lng.bit_length() > 63:
        raise ValueError("expected 64 bit long {}".format(lng))
    if lng < 0:
        lng = (1 << 64) + lng
    return uuid.UUID(int=lng)


def type_check(value, value_type, variable, override=None):
    """Check value is of type value_type, or is override"""

    sys.tracebacklimit = None
    # int_32 is an internal type so for the purpose of this external checking
    # it's OK to use the parent type (int).
    vtype = value_type
    if value_type == int_32:
        vtype = int
    if not isinstance(value, vtype) and value is not override:
        if override is not False:
            raise TypeError(
                "Expected {} or {} for {}, received {}".format(vtype,
                                                               override,
                                                               variable,
                                                               type(value)))
        else:
            raise TypeError(
                "Expected {} for {}, received {}".format(vtype,
                                                         variable,
                                                         type(value)))
    return True


def value_check(value, min_value, max_value, variable, override=None):
    """Check if value is inside allowed range, or override"""

    sys.tracebacklimit = None
    if value != override and not (
                (min_value is not None) and (min_value <= value) and
            (max_value is not None) and (value <= max_value)):
        if override is not False:
            raise ValueError(
                "Expected {} <= {} <= {} or {}, received {}".format(
                    min_value, variable, max_value, override, value))
        else:
            raise ValueError(
                "Expected {} <= {} <= {}, received {}".format(
                    min_value, variable, max_value, value))

    return True


def get_differences(expected, actual, parent=None):
    """
    Compare two entities. Provide a report if differences exist between the two
    entities.

    :param expected: What is expected to exist. May be AbstractCaomEntity,
    TypedOrderedDict, TypedList, or TypedSet.
    :param actual: What exists. May be AbstractCaomEntity,
    TypedOrderedDict, TypedList, or TypedSet.
    :return: None if the entities are the same, or a text report of the
    individual differences.
    """
    report = []

    if expected == actual:
        return None

    if type(expected) != type(actual):
        report.append(
            'Types:: expected \'{}\' actual \'{}\''.format(type(expected),
                                                           type(actual)))
        return report

    if (isinstance(expected, TypedOrderedDict) or
        isinstance(expected, TypedList) or
            isinstance(expected, TypedSet)):
        temp_report = _get_collection_differences(expected, actual, parent)
    else:
        assert isinstance(expected, AbstractCaomEntity)
        assert isinstance(actual, AbstractCaomEntity)
        temp_report = _get_object_differences(expected, actual, parent)

    if temp_report:
        report.extend(temp_report)

    return report if len(report) > 0 else None


def _get_object_differences(expected, actual, parent=None):
    """Reports on the differences between both attributes and
    their values for object differences."""
    report = []
    expected_dict, expected_decompose = _get_dict(expected)
    actual_dict, actual_decompose = _get_dict(actual)

    if expected_dict != actual_dict:
        new_parent = '{}.{}'.format(
            parent, expected.__class__) if parent else expected.__class__
        temp_report = _get_dict_differences(expected_dict, actual_dict, new_parent)
        if temp_report:
            report.extend(temp_report)

    for expected_key, expected_value in expected_decompose.items():
        if expected_key in actual_decompose:
            actual_value = actual_decompose[expected_key]
            label = '{}.{}'.format(parent, expected_key)
            temp_report = get_differences(expected_value, actual_value, label)
            actual_decompose.pop(expected_key)
            if temp_report:
                report.extend(temp_report)
        else:
            report.append('Member:: {} missing from {}'.format(
                expected_key, actual.__class__))

    for actual_key in actual_decompose.items():
        report.append(
            'Member::{} missing from {}'.format(actual_key, expected.__class__))

    return report if len(report) > 0 else None


def _get_collection_differences(expected, actual, parent=None):
    """Reports on the differences between two collections. Ignores collection
    ordering."""
    report = []
    if len(expected) != len(actual):
        report.append(
            'Collection:: length of expected {} != actual {}'.format(
                len(expected), len(actual)))

    for expected_key, expected_value in expected.items():
        label = '{}[\'{}\']'.format(parent, expected_key)
        if expected_key in actual.keys():
            actual_value = actual[expected_key]
            temp_report = get_differences(expected_value,
                                          actual_value, label)
            actual.pop(expected_key)
            if temp_report:
                report.extend(temp_report)
            break
        else:
            report.append('Collection:: {} not in actual.'.format(label))

    for key in actual.keys():
        label = '{}[\'{}\']'.format(parent, key)
        report.append('Collection:: actual {} not in expected.'.format(
            label))

    return report if len(report) > 0 else None


def _get_dict_differences(expected, actual, parent):
    """Reports on how two dictionaries are different."""
    report = []
    for expected_key, expected_value in expected.items():
        if expected_key in actual:
            actual_value = actual[expected_key]
            if _not_equal(expected_value, actual_value):
                report.append(
                    'Member value:: {}.{}: expected {} actual {}'.format(
                        parent,
                        expected_key,
                        expected_value,
                        actual_value))
            actual.pop(expected_key)
        else:
            report.append(
                'Member:: {}.{}: expected missing from actual.'.format(
                    parent, expected_key))

    for key in actual.items():
        report.append(
            'Member:: {}.{}: actual missing from expected.'.format(parent,
                                                                   key))

    return report if len(report) > 0 else None


def _not_equal(rhs, lhs):
    """Handling for not-quite-equals float comparisons."""
    if isinstance(rhs, float) and isinstance(lhs, float):
        result = math.isclose(rhs, lhs, rel_tol=1e-10)
    else:
        result = rhs == lhs
    return not result


def _get_dict(entity):
    """This removes all the entity attributes that are not considered part of
    an AbstractCaomEntity comparison, and tracks the entities that are not
    straight-to-dictionary conversions for later handling."""

    attributes = {}
    caom_collections = {}
    for i in dir(entity):
        attribute = getattr(entity, i)
        if (i.startswith('_') or
                callable(attribute) or
                i.find('checksum') != -1):
            continue
        if (isinstance(attribute, TypedOrderedDict) or
                isinstance(i, TypedList) or
                isinstance(i, TypedSet) or
                isinstance(i, AbstractCaomEntity)):
            caom_collections[i] = attribute
        else:
            attributes[i] = getattr(entity, i)
    return attributes, caom_collections


class TypedList(collections.MutableSequence):
    """
    Class that implements a typed list in Python. Supported types
    are specified when instance is created. Example:

       obstype = TypedList((str), "calibration", "science")
       emptylist = TypedList((str), )
       multipletypes = TypedList((str, int, int), 1, 12L, "ABC")

    An AssertionError is thrown when the caller attempts to insert
    objects with the wrong type.
    """

    def __init__(self, oktypes, *args):
        """
        Initializes a TypedList.

        Arguments:
        oktypes : a list of types accepted in this list
        *args : values to be stored in the list

        AssertionError will be raised if the types of the arguments
        are not found in the list of oktypes.
        """
        self._oktypes = oktypes
        self.list = list()
        self.extend(list(args))

    def __str__(self):
        return "\n".join(["{}".format(v) for v in self])

    def __repr__(self):
        return "TypedList((%r))," % self._oktypes + (
            "(".join(["(%r)" % v for v in self]) + ")")

    def check(self, v):
        assert isinstance(v, self._oktypes), (
            "Wrong type in list. OK Types: {0}".format(self._oktypes))

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

    def __setitem__(self, i, v):
        self.check(v)
        self.list[i] = v

    def insert(self, i, v):
        self.check(v)
        self.list.insert(i, v)

    @property
    def key_type(self):
        """ Returns the type of the elements of this list"""
        return self._oktypes


class TypedSet(collections.MutableSet):
    """
    Class that implements a typed set in Python. Supported types
    are specified when instance is created. Example:

       obstype = TypedSet((str), "calibration", "science")
       emptyset = TypedSet((str), )
       multipletypes = TypedSet((str, int, int), 1, 12L, "ABC")

    An AssertionError is thrown when the caller attempts to insert
    objects with the wrong type.
    """

    def __init__(self, oktypes, *args):
        """
        Initializes a TypedSet.

        Arguments:
        oktypes : a list of types accepted in this set
        *args : values to be stored in the set

        AssertionError will be raised if the types of the arguments
        are not found in the list of oktypes.
        """
        self._oktypes = oktypes
        self._set = set()
        for arg in args:
            self.add(arg)

    def check(self, v):
        assert isinstance(v, self._oktypes), (
            "Wrong type in list. OK Types: {0}".format(self._oktypes))

    def add(self, v):
        """Add an element."""
        self.check(v)
        self._set.add(v)

    def update(self, other):
        """Add all elements in other."""
        for v in other:
            self.check(v)
        self._set.update(other)

    def discard(self, v):
        """Remove an element. Do not raise an exception if absent."""
        self._set.discard(v)

    @property
    def key_type(self):
        """ Returns the type of the elements of this set"""
        return self._oktypes

    def __iter__(self):
        return iter(self._set)

    def __len__(self):
        return len(self._set)

    def __contains__(self, item):
        try:
            return item in self._set
        except AttributeError:
            return False


class TypedOrderedDict(collections.OrderedDict):
    """
    Class that implements a typed ordered dictionary in Python. Supported
    type is specified when an instance is created. A supported type
    must have a key property that is used as the key in the dictionary.

    Example:

       obsType = TypedOrderedDict((Plane), ("productID1", Plane('productID1')))
       emptyDict = TypedSet((Plane), )

    A TypeError is thrown when the caller attempts to insert objects
    with the wrong type.
    """

    def __init__(self, key_type=None, *args):
        """
        Initializes a TypedOrderedDict.

        Arguments:
        keyType : the type values that will be stored in this dictionary.
        NOTE: the default values has been added only because the
        parent class (OrderDict) requires a default constructor in
        copydeep operations.
        *args : (key, value) tuples to be stored in the dictionary
        (value must match type)

        AssertionError will be raised if the type of the arguments
        is not the same as keyType.
        """
        super(TypedOrderedDict, self).__init__(self)
        self._oktypes = key_type
        for arg in args:
            self.__setitem__(arg[0], arg[1])

    def __str__(self):
        return "\n".join(["{} => {}".format(k, v)
                          for k, v in six.iteritems(self)])

    def __repr__(self):
        return "TypedOrderedDict((%r))," % self._oktypes + (
            "(".join(
                ["(%r,%r)" % (k, v) for k, v in six.iteritems(self)]) + ")")

    def check(self, key, value):
        """
        Check that the value is of the correct type for this typed
        dictionary and that the key attribute of value matches the 'key'
        passed to check
        """

        if not isinstance(value, self._oktypes):
            raise TypeError("Wrong type in dictionary. Key Type: {0}"
                            .format(self._oktypes))
        if not value._key() == key:
            raise ValueError("Key mismatch, key in object: {0}"
                             .format(value._key()))

    @property
    def key_type(self):
        """ Returns the type of the elements of this dictionary"""
        return self._oktypes

    def add(self, value):
        """Adds a new entry, value must be of type {0}"""
        self.__setitem__(value._key(), value)

    def pop(self, key):
        """Removes an entry"""
        return super(TypedOrderedDict, self).pop(key)

    def __setitem__(self, key, value):
        """Add an item."""
        self.check(key, value)
        super(TypedOrderedDict, self).__setitem__(key, value)


class ClassProperty(property):
    """ """

    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()
