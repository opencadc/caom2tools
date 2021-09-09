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
validate performs a validation on a CAOM2 element. By default deep=True
triggering a validation of all the sub-elements in the CAOM2 tree.
"""


import logging

from caom2 import Observation, Plane, Artifact, Part, Chunk
from caom2utils.polygonvalidator import validate_polygon
from caom2utils.wcsvalidator import validate_wcs


__all__ = ['validate']


logger = logging.getLogger('caomvalidator')


def validate(caom2_entity, deep=True):
    """
    Perform validation of the content of a CAOM element.

    Throws AssertionError if validation fails.

    :param caom2_entity: CAOM element to perform the validation on. It
    could be Observation, Plane, Artifact, Part or Chunk
    :param deep if True, also validate the 'has-a' members of an element.
    """
    if caom2_entity is not None:
        if isinstance(caom2_entity, Observation):
            _validate_observation(caom2_entity, deep)
        elif isinstance(caom2_entity, Plane):
            _validate_plane(caom2_entity, deep)
        elif isinstance(caom2_entity, Artifact):
            _validate_artifact(caom2_entity, deep)
        elif isinstance(caom2_entity, Part):
            _validate_part(caom2_entity, deep)
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
    _check_param(caom2_entity, Observation)
    if caom2_entity.proposal:
        _validate_keyword('proposal.keywords', caom2_entity.proposal.keywords)
    if caom2_entity.target:
        _validate_keyword('target.keywords', caom2_entity.target.keywords)
    if caom2_entity.telescope:
        _validate_keyword('telescope.keywords',
                          caom2_entity.telescope.keywords)
    if caom2_entity.instrument:
        _validate_keyword('telescope.instrument',
                          caom2_entity.instrument.keywords)
    if deep:
        for plane in caom2_entity.planes.values():
            _validate_plane(plane)


def _validate_plane(caom2_entity, deep=True):
    """
    Perform validation of the content of a Plane.

    Throws AssertionError if the members of the plane metadata do not meet
    structure criteria.

    :param caom2_entity: The Plane to validate.
    :param deep if True, also validate the 'has-a' members of a Plane.
    """
    _check_param(caom2_entity, Plane)
    if caom2_entity.provenance:
        _validate_keyword('provenance.keywords',
                          caom2_entity.provenance.keywords)
    if caom2_entity.position:
        validate_polygon(caom2_entity.position.bounds)

    if deep:
        for artifact in caom2_entity.artifacts.values():
            _validate_artifact(artifact)


def _validate_artifact(caom2_entity, deep=True):
    """
    Perform validation of the content of an Artifact.

    Throws AssertionError if the members of the artifact metadata do not meet
    structure criteria.

    :param caom2_entity: The Artifact to validate.
    """
    _check_param(caom2_entity, Artifact)
    if deep and caom2_entity.parts is not None:
        for value in caom2_entity.parts.values():
            _validate_part(value)


def _validate_part(caom2_entity, deep=True):
    """
    Perform validation of the content of a Part.

    Throws AssertionError if the members of the part metadata do not meet
    structure criteria.

    :param caom2_entity: The Part to validate.
    """
    _check_param(caom2_entity, Part)
    if deep and (caom2_entity.chunks is not None):
        for chunk in caom2_entity.chunks:
            _validate_chunk(chunk)


def _validate_chunk(caom2_entity):
    """
    Perform validation of the content of a Chunk.

    Throws AssertionError if the members of the chunk metadata do not meet
    structure criteria.

    :param caom2_entity: The Chunk to validate.
    """
    _check_param(caom2_entity, Chunk)
    validate_wcs(caom2_entity)


def _validate_keyword(name, keywords):
    """
    Perform validation for a Set of keywords.

    :param name: Logging parameter to identify where the keywords originate.
    :param keywords: Set of keyword values to validate.
    """
    if not keywords:
        return
    for keyword in keywords:
        if keyword is not None and keyword.find('|') != -1:
            raise AssertionError(
                f'invalid {name}: may not contain pipe (|)')


def _check_param(param, param_type):
    if param is None or not isinstance(param, param_type):
        raise ValueError(
            f'{param} must be a valid {param_type.__name__}.')
