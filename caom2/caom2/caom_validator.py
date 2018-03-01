# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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
This module contains the functionlity for validation of entities in CAOM2.
Two types of validation are available:
    validate - returns the validation of the entity, except its CAOM2
               children
    validate_acc - returns the accumulated validation of the entity attributes
                   including those of the CAOM2 children.

The validation represents the correctness of the state of the entity itself,
while the accumulated validation represents the state of the entity, and
all the children below.

IMPORTANT NOTE: The validate algorithms use introspection to automatically
find the attributes that are part of the CAOM2 model. It is therefore very
important that attributes that are not part of the model
(http://www.opencadc.org/caom2) be prefixed with an '_' so that the validate
algorithms ignore them in their computations.

Perform validation of the content of an Observation."""

from __future__ import (absolute_import, print_function, unicode_literals)

from builtins import str, int


import logging
import uuid

from aenum import Enum
from datetime import datetime

from caom2.caom_util import TypedList, TypedOrderedDict, TypedSet, int_32
from caom2.common import CaomObject, ObservationURI
from caom2.common import ChecksumURI
from caom2.observation import Observation


__all__ = ['validate_obs', 'validate', 'validate_acc']


logger = logging.getLogger('caom_validate')


def validate(entity, acc=False, name=None):
    """
    Validate the content of any CAOM2 entity.

    Throws an AssertionError if the validate call fails.

    :param entity: CAOM2 entity, descended from CaomObject.
    :param acc: Boolean to indicate whether to validate the entities members.
    :param name: Better logging for a better future.
    """

    if type(entity) is None:
        logger.debug('Validated empty entity {}'.format(entity.__name__))
        return

    if (isinstance(entity, ObservationURI) or isinstance(entity, ChecksumURI)
            or isinstance(entity, bytes) or isinstance(entity, bool)
            or isinstance(entity, float) or isinstance(entity, int_32)
            or isinstance(entity, int) or isinstance(entity, str)
            or isinstance(entity, datetime) or isinstance(entity, Enum)
            or isinstance(entity, uuid.UUID)):
        return
    elif isinstance(entity, CaomObject):
        if isinstance(entity, Observation):
            validate_obs(entity)
        if hasattr(entity, 'validate') and callable(
                getattr(entity, 'validate')):
            logger.debug('Direct invocation of validate for {}'.format(
                entity.__class__.__name__))
            entity.validate()
        if acc:
            validate_acc(entity)
    elif isinstance(entity, set) or \
            isinstance(entity, list) or isinstance(entity, TypedList) or \
            (isinstance(entity, TypedSet) and not
                isinstance(entity.key_type, CaomObject)):
        for i in entity:
            validate(i, acc)
    elif isinstance(entity, TypedOrderedDict):
        for i in entity:
            validate(entity[i], acc)
    else:
        raise AssertionError(
            'Cannot find a validate behaviour for {} of type {}.'.format(
                name, type(entity)))


def validate_acc(entity):
    """
    Validate the content of any CAOM2 entity, and it's children.

    Throws an AssertionError if any of the validate calls fail.

    :param entity: CAOM2 entity.
    """
    assert isinstance(entity, CaomObject)

    for i in dir(entity):
        if not callable(getattr(entity, i)) and not i.startswith('_'):
            attribute = getattr(entity, i)
            if attribute is not None:
                logger.debug('Calling validate for {} of type {}'.format(
                    i, type(attribute)))
                validate(attribute, True, i)


def validate_obs(observation):
    """
    Perform validation of the content of an Observation.

    Throws AssertionError if the keywords that are members of various
    observation and plane metadata do not meet structure criteria.

    :param observation: The Observation (SimpleObservation or
        CompositeObservation) to validate.
    """
    assert isinstance(observation, Observation)
    if observation is not None:
        _validate_keywords(observation)


def _validate_keywords(observation):
    """
    Perform validation for all the places where keywords are members of
    Observations and Planes.

    Throws AssertionError if the keywords that are members of various
    observation and plane metadata do not meet structure criteria.

    :param observation: The parent Observation (SimpleObservation or
        CompositeObservation) to validate.
    """
    if observation.proposal is not None:
        _validate_keyword('proposal.keywords', observation.proposal.keywords)
    if observation.target is not None:
        _validate_keyword('target.keywords', observation.target.keywords)
    if observation.telescope is not None:
        _validate_keyword('telescope.keywords', observation.telescope.keywords)
    if observation.instrument is not None:
        _validate_keyword('instrument.keywords',
                           observation.instrument.keywords)
    for plane_key in observation.planes.keys():
        if observation.planes[plane_key].provenance is not None:
            _validate_keyword('provenance.keywords',
                              observation.planes[plane_key].provenance.keywords)


def _validate_keyword(name, keywords):
    """
    Perform validation for a Set of keywords.

    :param name: Logging parameter to identify where the keywords originate.
    :param keywords: Set of keyword values to validate.
    """
    for keyword in keywords:
        _assert_validate_keyword('', name, keyword)


def _assert_validate_keyword(caller, name, keyword):
    """
    Keywords can contain any valid UTF-8 character except the pipe (|). The
    pipe character is reserved for use as a separator in persistence
    implementations so the list of keywords can be serialized in a single string
    to support querying.

    :param caller: Logging parameter - identifies which class called this
        method.
    :param name: Logging parameter - narrows the calling identification for
        this method.
    :param keyword: The parameter to check.
    """
    if keyword.find('|') != -1:
        raise AssertionError(
            '{}: invalid {}: may not contain pipe (|)'.format(caller, name))
