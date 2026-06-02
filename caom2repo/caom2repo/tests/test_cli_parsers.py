# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2026.                            (c) 2026.
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
#  General Public License for           Générale Publique GNU Affero pour
#  more details.                        plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est pas le cas,
#  <http://www.gnu.org/licenses/>.      consultez :
#                                       <http://www.gnu.org/licenses/>.
#
# ***********************************************************************

"""Contract tests for caom2-repo CLI parsers."""

import pytest

from caom2repo.core import DEFAULT_RESOURCE_ID, build_parser
from cadcutils.util.tests.parser_helpers import (
    assert_epilog_contains, assert_has_base_dests, assert_has_dests,
    assert_help_contains, get_subparser, subparser_names,
)

_SUBCOMMANDS = ('create', 'read', 'update', 'delete', 'visit')


@pytest.fixture
def root_parser():
    return build_parser()


def test_root_parser_contract(root_parser):
    assert set(subparser_names(root_parser)) == set(_SUBCOMMANDS)
    assert_help_contains(
        root_parser,
        'Client for a CAOM2 repo',
        'CRUD',
    )


@pytest.mark.parametrize('subcmd,extra_dests', [
    ('create', ('observation',)),
    ('read', ('output', 'collection', 'observationID')),
    ('update', ('observation',)),
    ('delete', ('collection', 'observationID')),
    ('visit', ('plugin', 'start', 'end', 'obs_file', 'threads',
               'halt_on_error', 'collection')),
])
def test_subcommand_parser_contract(root_parser, subcmd, extra_dests):
    parser = get_subparser(root_parser, subcmd)
    assert_has_base_dests(parser)
    assert_has_dests(parser, *extra_dests)
    assert_help_contains(parser, DEFAULT_RESOURCE_ID)


def test_visit_epilog(root_parser):
    parser = get_subparser(root_parser, 'visit')
    assert_epilog_contains(
        parser,
        'ObservationUpdater',
        'def update(self, observation, **kwargs)',
    )
