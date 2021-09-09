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


import os

from caom2utils import validate
from caom2utils.caomvalidator import _validate_keyword, _check_param
from caom2 import ObservationReader
from caom2 import SimpleObservation, DerivedObservation, Proposal
from caom2 import Algorithm, Telescope, Instrument, Target
from caom2 import Plane, Provenance

import pytest


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA = 'data'


def test_assert_validate_keyword():
    _validate_keyword('test', 'foo')
    _validate_keyword('test', 'foo=42')
    _validate_keyword('test', 'foo:42')
    _validate_keyword('test', "tick'marks")
    _validate_keyword('test', 'has multiple spaces')
    exception_raised = False
    try:
        _validate_keyword('test', 'pipe|denied')
    except AssertionError:
        # successful test case
        exception_raised = True
    assert exception_raised


def test_validate_observation():
    obs = SimpleObservation('test_collection', 'test_obs_id',
                            Algorithm('test_name'))
    validate(obs)
    obs = DerivedObservation('test_collection', 'test_obs_id',
                             Algorithm('test_name'),
                             proposal=Proposal('test_proposal'),
                             telescope=Telescope('test_telescope'),
                             instrument=Instrument('test_instrument'),
                             target=Target('test_targets'))
    obs.algorithm.keywords = 'foo'
    obs.proposal.keywords = set('foo=42')
    obs.telescope.keywords = set('foo:42')
    obs.instrument.keywords.add("tick'marks")
    obs.target.keywords = set('has multiple spaces')
    test_plane = Plane('test_plane')
    test_plane.provenance = Provenance('test_provenance')
    test_plane.provenance.keywords.add('pipe|denied')
    obs.planes['test_plane'] = test_plane
    with pytest.raises(AssertionError):
        validate(obs)


def test_compatibility():
    # tests a previously generated observation and validates the
    # entities, and the entities with children

    source_file_path = os.path.join(THIS_DIR, TEST_DATA,
                                    'SampleComposite-CAOM-2.3.xml')
    reader = ObservationReader(True)
    with open(source_file_path):
        obs = reader.read(source_file_path)

    # shallow validates first
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    validate(chunk, False)
                validate(part, False)
            validate(artifact, False)
        validate(plane, False)

    validate(obs, False)

    # deep validate
    validate(obs, True)


def test_failures():
    test_object = type('', (), {})()
    with pytest.raises(AssertionError):
        validate(test_object)

    with pytest.raises(ValueError):
        _check_param(test_object, SimpleObservation)
