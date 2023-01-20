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

"""
A difference method for CAOM2 entities.

"""

import math

from caom2.common import CaomObject
from caom2.caom_util import TypedSet, TypedOrderedDict, TypedList
from caom2 import Chunk
from . import caom_util


__all__ = ['get_differences']


# a list of dict keys that are not checked as a result of a get_differences
# call:
# 'meta_producer' is ignored because it can cause a lot of noise in github
#                 commits for not a lot of value
#
ignore_keys = ['meta_producer']


def get_differences(expected, actual, parent=None):
    """
    Compare two entities. Provide a report if differences exist between the two
    entities.

    :param expected: What is expected to exist. May be AbstractCaomEntity,
    TypedOrderedDict, TypedList, or TypedSet.
    :param actual: What exists. May be AbstractCaomEntity,
    TypedOrderedDict, TypedList, or TypedSet.
    :param parent: str for message content.
    :return: None if the entities are the same, or a text report of the
    individual differences.
    """
    report = []

    if type(expected) != type(actual):
        report.append(
            'Types:: expected \'{}\' actual \'{}\''.format(type(expected),
                                                           type(actual)))
        return report

    if parent:
        parent = '{}.{}'.format(parent, expected.__class__.__name__)
    else:
        parent = expected.__class__.__name__
    temp_report = None
    if (isinstance(expected, TypedOrderedDict) or
        isinstance(expected, TypedList) or
            isinstance(expected, TypedSet) or
            isinstance(expected, set) or
            isinstance(expected, list)):
        temp_report = _get_collection_differences(expected, actual, parent)
    elif isinstance(expected, CaomObject):
        caom_util.type_check(actual, CaomObject, "CaomObject")
        temp_report = _get_object_differences(expected, actual, parent)
    else:
        if expected != actual:
            temp_report = ['Value:: {} expected is {}, actual is {}'.format(
                expected, actual, parent)]

    if temp_report:
        report.extend(temp_report)

    return report if len(report) > 0 else None


def _get_object_differences(expected, actual, parent):
    """Reports on the differences between both attributes and
    their values for object differences."""
    report = []
    expected_dict, expected_decompose = _get_dict(expected)
    actual_dict, actual_decompose = _get_dict(actual)

    if expected_dict != actual_dict:
        temp_report = _get_dict_differences(expected_dict, actual_dict, parent)
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
            report.append('Member:: {}.{}: missing from {}'.format(
                parent, expected_key, actual.__class__))

    for actual_key in actual_decompose.items():
        report.append('Member:: {}.{}: missing from {}'.format(
            parent, actual_key, expected.__class__))

    return report if len(report) > 0 else None


def _get_collection_differences(expected, actual, parent):
    """Reports on the differences between two collections. Ignores collection
    ordering."""
    report = []
    if not expected and actual:
        report.append('Collection:: {} expected not found'.format(parent))
        return report

    if expected and not actual:
        report.append('Collection:: {} actual not found'.format(parent))
        return report

    if not expected and not actual:
        return report

    if len(expected) != len(actual):
        report.append(
            'Collection:: {}: length of expected {} != length of actual {}'.
            format(parent, len(expected), len(actual)))

    if (isinstance(actual, TypedList) or isinstance(actual, TypedSet)
            or isinstance(actual, set)):
        temp_report = _get_sequence_differences(expected, actual, parent)
    elif isinstance(actual, list):
        temp_report = _get_list_differences(expected, actual, parent)
    else:
        temp_report = _get_mapping_differences(expected, actual, parent)

    if temp_report:
        report.extend(temp_report)

    return report if len(report) > 0 else None


def _get_mapping_differences(expected, actual, parent):
    """Reports on how two maps are different."""
    report = []
    actual_keys_copy = list(actual.keys())

    for expected_key, expected_value in expected.items():
        label = '{}[\'{}\']'.format(parent, expected_key)
        if expected_key in actual.keys():
            actual_value = actual[expected_key]
            temp_report = get_differences(expected_value,
                                          actual_value, label)
            actual_keys_copy.remove(expected_key)
            if temp_report:
                report.extend(temp_report)
        else:
            report.append('Map:: {} not in actual.'.format(label))

    for key in actual_keys_copy:
        label = '{}[\'{}\']'.format(parent, key)
        report.append('Map:: actual {} not in expected.'.format(
            label))

    return report if len(report) > 0 else None


def _get_list_differences(expected, actual, parent):
    """Reports on how two lists are different. A lot like
    _get_sequence_differences, except that the second for loop will start
    at an advanced index, and it does not expect Chunk-type entities."""

    report = []

    # this method removes list elements, so make copies for modification,
    # leaving the originals unchanged - kind of immutable inputs
    actual_copy = list(actual)
    expected_copy = list(expected)

    start_index = 0
    for ex_index, e in enumerate(expected):
        label = '{}[\'{}\']'.format(parent, ex_index)
        match_found = False
        tracking_report = None
        tracking_actual = None
        tracking_expected = None

        for act_index, a in enumerate(actual):
            # index comparison skips already-matched list entries
            if act_index >= start_index:
                temp_report = get_differences(e, a, label)
                if temp_report is None:
                    match_found = True
                    start_index = act_index + 1
                    if a in actual_copy:
                        actual_index = actual_copy.index(a)
                        actual_copy.pop(actual_index)
                    if e in expected_copy:
                        expected_index = expected_copy.index(e)
                        expected_copy.pop(expected_index)
                    break
                else:
                    # every non-matching comparison will have failures,
                    # so pick the report with the least number of errors
                    if tracking_report:
                        if len(temp_report) < len(tracking_report):
                            tracking_report = temp_report
                            tracking_actual = a
                            tracking_expected = e
                    else:
                        tracking_report = temp_report
                        tracking_actual = a
                        tracking_expected = e

        if not match_found:
            report.extend(tracking_report)
            if tracking_actual in actual_copy:
                actual_index = actual_copy.index(tracking_actual)
                actual_copy.pop(actual_index)
            if tracking_expected in expected_copy:
                expected_index = expected_copy.index(tracking_expected)
                expected_copy.pop(expected_index)

    for e in enumerate(expected_copy):
        label = '{}[\'{}\']'.format(parent, e)
        report.append(
            'List:: {} expected not found in actual'.format(label))

    for a in enumerate(actual_copy):
        label = '{}[\'{}\']'.format(parent, a)
        report.append(
            'List:: {} actual not found in expected'.format(label))

    return report if len(report) > 0 else None


def _get_sequence_differences(expected, actual, parent):
    """Reports on how two sequences are different."""
    report = []

    # this method removes sequence elements, so make copies of the collections
    # under comparison and modify the copies, leaving the originals
    # unchanged - kind of immutable inputs
    actual_copy = list(actual)
    expected_copy = list(expected)

    for ex_index, e in enumerate(expected):
        label = '{}[\'{}\']'.format(parent, ex_index)
        if isinstance(e, Chunk):
            temp_report = get_differences(e, actual[0])
            if temp_report is not None:
                report.extend(temp_report)
            return report if len(report) > 0 else None
        else:
            match_found = False
            tracking_report = None
            tracking_actual = None
            tracking_expected = None
            for act_index, a in enumerate(actual):
                temp_report = get_differences(e, a, label)
                if temp_report is None:
                    match_found = True
                    actual_index = actual_copy.index(a)
                    actual_copy.pop(actual_index)
                    expected_index = expected_copy.index(e)
                    expected_copy.pop(expected_index)
                    break
                else:
                    # pick the report with the least number of errors
                    if tracking_report:
                        if len(temp_report) < len(tracking_report):
                            tracking_report = temp_report
                            tracking_actual = a
                            tracking_expected = e
                    else:
                        tracking_report = temp_report
                        tracking_actual = a
                        tracking_expected = e

            if not match_found:
                report.extend(tracking_report)
                if tracking_actual in actual_copy:
                    actual_index = actual_copy.index(tracking_actual)
                    actual_copy.pop(actual_index)
                if tracking_expected in expected_copy:
                    expected_index = expected_copy.index(tracking_expected)
                    expected_copy.pop(expected_index)

    for e in enumerate(expected_copy):
        label = '{}[\'{}\']'.format(parent, e)
        report.append(
            'Sequence:: {} expected not found in actual'.format(label))

    for a in enumerate(actual_copy):
        label = '{}[\'{}\']'.format(parent, a)
        report.append(
            'Sequence:: {} actual not found in expected'.format(label))

    return report if len(report) > 0 else None


def _get_dict_differences(expected, actual, parent):
    """Reports on how two dictionaries are different."""
    report = []
    for expected_key, expected_value in expected.items():
        found = False
        for ignore_key in ignore_keys:
            if expected_key == ignore_key:
                actual.pop(expected_key)
                found = True
        if found:
            continue
        if expected_key in actual:
            actual_value = actual[expected_key]
            if _is_composite_instance_type(actual_value):
                temp_report = get_differences(expected_value, actual_value,
                                              expected_key)
                if temp_report:
                    report.append(temp_report)
            elif _not_equal(expected_value, actual_value):
                report.append(
                    'Value:: {}.{}: expected {} actual {}'.format(
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
            'Member:: {}.{}: unexpected.'.format(parent, key))

    return report if len(report) > 0 else None


def _not_equal(rhs, lhs):
    """Handling for not-quite-equals float comparisons."""
    if isinstance(rhs, float) and isinstance(lhs, float):
        if math.isnan(rhs) or math.isnan(lhs):
            if math.isnan(rhs) and math.isnan(lhs):
                result = True  # the opposite of what python will say
            else:
                result = rhs == lhs
        else:
            # if only using python 3.5+, use math.isclose, instead of this
            # description of math.isclose from the python documentation
            result = abs(rhs-lhs) <= max(1e-12 * max(abs(rhs), abs(lhs)), 1e-11)
    else:
        result = rhs == lhs
    return not result


def _get_dict(entity):
    """This removes all the entity attributes that are not considered part of
    an CaomObject comparison, and tracks the entities that are not
    straight-to-dictionary conversions for later handling."""

    attributes = {}
    caom_collections = {}
    for i in dir(entity):
        try:
            attribute = getattr(entity, i)
        except TypeError:
            pass
        if (i.startswith('_') or
                callable(attribute) or
                i.find('checksum') != -1 or
                i.find('last_modified') != -1):
            # ignore the built-ins,
            # checksum and last modified will not be the same
            continue
        if _is_composite_instance_type(attribute):
            caom_collections[i] = attribute
        else:
            attributes[i] = attribute
    return attributes, caom_collections


def _is_composite_instance_type(entity):
    """Common location to test for any of the varying types that might make
    up a CAOM2 entity. """
    return (isinstance(entity, TypedOrderedDict) or
            isinstance(entity, TypedList) or
            isinstance(entity, TypedSet) or
            isinstance(entity, CaomObject) or
            isinstance(entity, set) or
            isinstance(entity, list))
