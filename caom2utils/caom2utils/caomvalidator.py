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

import logging
import six
import numpy as np
from spherical_geometry import polygon
from aenum import Enum
import numpy as np

from caom2 import Observation, Plane, Artifact, Part, Chunk, Provenance
from caom2 import Point, Polygon, MultiPolygon, Polarization, SegmentType
from caom2utils import validate_wcs


__all__ = ['validate']


logger = logging.getLogger('caom_validate')


def validate(caom2_entity, deep=True):
    if caom2_entity is not None:
        if isinstance(caom2_entity, Observation):
            _validate_observation(caom2_entity, deep)
        elif isinstance(caom2_entity, Plane):
            _validate_plane(caom2_entity, deep)
        elif isinstance(caom2_entity, Artifact):
            _validate_artifact(caom2_entity)
        elif isinstance(caom2_entity, Part):
            _validate_part(caom2_entity)
        elif isinstance(caom2_entity, Chunk):
            _validate_chunk(caom2_entity)
        else:
            raise AssertionError("Not a CAOM2 entity")


def _validate_observation(caom2_entity, deep=True):
    """
    Perform validation of the content of an Observation.

    Throws AssertionError if the keywords that are members of various
    observation and plane metadata do not meet structure criteria.

    :param caom2_entity: The Observation (SimpleObservation or
        CompositeObservation) to validate.
    :param deep if True, also validate the 'has-a' members of an Observation.
    """
    assert isinstance(caom2_entity, Observation), 'Must be an Observation'
    if caom2_entity is not None:
        _validate_keywords(caom2_entity)

    if deep:
        pass  # TODO


def _validate_plane(caom2_entity, deep=True):
    """
    Perform validation of the content of a Plane.

    Throws AssertionError if the members of the plane metadata do not meet
    structure criteria.

    :param caom2_entity: The Plane to validate.
    :param deep if True, also validate the 'has-a' members of a Plane.
    """
    assert caom2_entity, 'Must provide a plane for validation'
    assert isinstance(caom2_entity, Plane), 'Must be a Plane'

    if deep:
        if caom2_entity.provenance is not None:
            _validate_provenance(caom2_entity.provenance)
        # Quality member is an enumerated type
        # all Metrics members are already validated by constructor
        if (caom2_entity._position is not None and
                caom2_entity._position.bounds is not None):
            _validate_bounds(caom2_entity._position)
        if (caom2_entity._energy is not None and
                caom2_entity._energy.bounds is not None):
            _validate_bounds(caom2_entity._energy)
        if (caom2_entity._time is not None and
                caom2_entity._time.bounds is not None):
            _validate_bounds(caom2_entity._time)
        if caom2_entity._polarization is not None:
            _validate_polarization(caom2_entity._polarization)


def _validate_artifact(caom2_entity):
    """
    Perform validation of the content of an Artifact.

    Throws AssertionError if the members of the artifact metadata do not meet
    structure criteria.

    :param caom2_entity: The Artifact to validate.
    """
    assert caom2_entity, 'Must provide an artifact for validation'
    assert isinstance(caom2_entity, Artifact), 'Must be a Artifact'
    if caom2_entity.parts is not None:
        for value in caom2_entity.parts.items():
            _validate_part(value)


def _validate_part(caom2_entity):
    """
    Perform validation of the content of a Part.

    Throws AssertionError if the members of the part metadata do not meet
    structure criteria.

    :param caom2_entity: The Part to validate.
    """
    assert caom2_entity, 'Must provide a part for validation'
    assert isinstance(caom2_entity, Part), 'Must be a Part'
    if caom2_entity.chunks is not None:
        for k, value in enumerate(caom2_entity.chunks):
            _validate_chunk(value)


def _validate_chunk(caom2_entity):
    """
    Perform validation of the content of a Chunk.

    Throws AssertionError if the members of the chunk metadata do not meet
    structure criteria.

    :param caom2_entity: The Chunk to validate.
    """
    assert caom2_entity, 'Must provide a chunk for validation'
    assert isinstance(caom2_entity, Chunk), 'Must be a Chunk'
    validate_wcs(caom2_entity)


def _validate_provenance(caom2_entity):
    """
    Perform validation of the content of a Provenance.

    Throws AssertionError if the members of the plane metadata do not meet
    structure criteria.

    :param caom2_entity: The Provenance to validate.
    """
    assert caom2_entity, 'Must provide a provenance for validation'
    assert isinstance(caom2_entity, Provenance), 'Must be a Provenance'
    if caom2_entity.keywords is not None:
        _validate_keyword('plane.provenance', caom2_entity.keywords)


def _validate_bounds(caom2_entity):
    """
    Perform validation of the content of a bounds entry.

    Throws AssertionError if the members of the bounds metadata do not meet
    structure criteria.

    :param caom2_entity: The bounds to validate.
    """
    assert caom2_entity, 'Must provide a bounds for validation'
    caom2_entity.validate()


def _validate_polarization(caom2_entity):
    """
    Perform validation of the content of a Polarization.

    Throws AssertionError if the members of the polarization metadata do not
    meet structure criteria.

    :param caom2_entity: The Polarization to validate.
    """
    assert caom2_entity, 'Must provide a polarization for validation'
    assert isinstance(caom2_entity, Polarization), 'Must be a Polarization'


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
        _validate_keyword(
            'instrument.keywords', observation.instrument.keywords)
    for plane_key in observation.planes.keys():
        if observation.planes[plane_key].provenance is not None:
            _validate_keyword(
                'provenance.keywords',
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
    implementations so the list of keywords can be serialized in a single
    string to support querying.

    :param caller: Logging parameter - identifies which class called this
        method.
    :param name: Logging parameter - narrows the calling identification for
        this method.
    :param keyword: The parameter to check.
    """
    if keyword.find('|') != -1:
        raise AssertionError(
            '{}: invalid {}: may not contain pipe (|)'.format(caller, name))


def _validate_multipolygon(mp):
    """
    Performs a basic validation of the current object.

    An AssertionError is thrown if the object does not represent a
    multi polygon
    """
    # perform a quick validation of this multipolygon object to fail early

    assert not isinstance(mp, MultiPolygon), \
        'MultiPoligon expected in validation received {}'.format(type(mp))

    _validate_size_and_end_vertices()

    # perform a more detailed validation of this multipolygon object
    mp_validator = MultiPolygonValidator()
    for i in range(len(mp.vertices)):
        mp_validator.validate(mp.vertices[i])


def _validate_size_and_end_vertices(self):
    if len(self._vertices) < 4:
        # triangle
        raise AssertionError('invalid polygon: {} vertices (min 4)'.format(
            len(self._vertices)))

    if self._vertices[0].type != SegmentType.MOVE:
        raise AssertionError(
            'invalid polygon: first vertex is not a MOVE vertex')

    if self._vertices[-1].type != SegmentType.CLOSE:
        raise AssertionError(
            'invalid polygon: last vertex is not a CLOSE vertex')


class MultiPolygonValidator():
    """
    A class to validate the sequencing of vertices in a polygon,
    as well as constructing and validating the polygon.

    An AssertionError is thrown if an incorrect polygon is detected.
    """

    def __init__(self):
        self._lines = 0
        self._open_loop = False
        self._polygon = Polygon()

    def validate(self, vertex):
        if vertex.type == SegmentType.MOVE:
            self.validate_move(vertex)
        elif vertex.type == SegmentType.CLOSE:
            self.validate_close(vertex)
        else:
            self.validate_line(vertex)

    def validate_move(self, vertex):
        if self._open_loop:
            raise AssertionError(
                'invalid polygon: MOVE vertex when loop open')
        self._lines = 0
        self._open_loop = True
        self._polygon.points.append(Point(vertex.cval1, vertex.cval2))

    def validate_close(self, vertex):
        # close the polygon
        if not self._open_loop:
            raise AssertionError(
                'invalid polygon: CLOSE vertex when loop close')
        if self._lines < 2:
            raise AssertionError(
                'invalid polygon: minimum 2 lines required')
        self._open_loop = False
        # SphericalPolygon requires point[0] == point[-1]
        point = self._polygon.points[0]
        self._polygon.points.append(Point(point.cval1, point.cval2))
        # validate the polygons in the multipolygon
        self._polygon.validate()
        # instantiate a new Polygon for the next iteration
        self._polygon = Polygon()

    def validate_line(self, vertex):
        if not self._open_loop:
            raise AssertionError(
                'invalid polygon: LINE vertex when loop close')
        self._lines += 1
        self._polygon.points.append(Point(vertex.cval1, vertex.cval2))