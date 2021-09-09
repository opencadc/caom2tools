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

from caom2utils import ObsBlueprint
from caom2utils.legacy import ConvertFromJava, load_config, apply_java_config
from caom2utils.legacy import _JAVA_CAOM2_CONFIG

import os
import pytest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')
cfhtwircam_override = os.path.join(TESTDATA_DIR, 'test.override')


@pytest.mark.parametrize('override_file', [cfhtwircam_override])
def test_class_apply_defaults(override_file):
    ob = ObsBlueprint(position_axes=(1, 2), energy_axis=3,
                      polarization_axis=4, time_axis=5)
    usc = {'Plane.dataProductType': 'plane.dataProductType',
           'Plane.provenance.producer': 'provenance.producer',
           'Plane.provenance.project': 'provenance.project',
           'Plane.metaRelease': 'plane.metaRelease',
           'Plane.dataRelease': 'plane.dataRelease',
           'Plane.calibrationLevel': 'plane.calibrationLevel',
           'Observation.metaRelease': 'obs.metaRelease',
           'Observation.intent': 'obs.intent',
           'Observation.type': 'obs.type',
           'Observation.proposal.pi': 'proposal.pi',
           'Observation.proposal.project': 'proposal.project',
           'Observation.proposal.title': 'proposal.title',
           'Observation.sequenceNumber': 'obs.sequenceNumber',
           'Observation.target.standard': 'target.standard',
           'Artifact.productType': 'artifact.productType',
           'Chunk.time.resolution': 'time.resolution',
           'Chunk.time.exposure': 'time.exposure',
           'Chunk.energy.resolvingPower': 'resolvingPower',
           'Chunk.energy.bandpassName': 'filtername',
           'Artifact.contentChecksum': 'artifact.contentChecksum'
           }

    convert = ConvertFromJava(ob, usc)
    test_overrides = load_config(override_file)

    for key, value in test_overrides.items():
        try:
            # artifacts is a substructure to be dealt with separately,
            # WCSAXES should work .... ;)
            if key == 'artifacts' or key == 'WCSAXES':
                continue

            result = convert.get_caom2_elements(key)
            for r in result:
                ob._get(r)
        except ValueError:
            assert False, f'Could not find key {key} in ObsBlueprint'


def test_apply_java_config():
    nonexistent_file = ''
    result = apply_java_config(nonexistent_file, use_only_defaults=True)
    assert result is not None, 'expect a result'
    assert result == _JAVA_CAOM2_CONFIG, 'should return the default'
