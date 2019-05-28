# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

import logging

from caom2 import TypedSet, ObservationURI, PlaneURI

from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc

__all__ = ['update_plane_provenance', 'update_observation_members']


def update_plane_provenance(plane, headers, lookup, collection,
                            repair, obs_id):
    """Add inputs to Planes, based on a particular keyword prefix.

    :param plane Plane instance to add inputs to
    :param headers FITS keyword headers that have lookup values.
    :param lookup The keyword pattern to find in the FITS header keywords for
        input files.
    :param collection The collection name for URI construction
    :param repair The function to fix input values, to ensure they match
        input observation ID values.
    :param obs_id String value for logging only.
    """
    plane_inputs = TypedSet(PlaneURI,)

    for header in headers:
        for keyword in header:
            if keyword.startswith(lookup):
                value = header.get(keyword)
                prov_obs_id, prov_prod_id = repair(value, obs_id)
                if prov_obs_id is not None and prov_prod_id is not None:
                    obs_member_uri_str = \
                        ec.CaomName.make_obs_uri_from_obs_id(
                            collection, prov_obs_id)
                    obs_member_uri = ObservationURI(obs_member_uri_str)
                    plane_uri = PlaneURI.get_plane_uri(
                        obs_member_uri, prov_prod_id)
                    plane_inputs.add(plane_uri)
                    logging.debug('Adding PlaneURI {}'.format(plane_uri))

    mc.update_typed_set(plane.provenance.inputs, plane_inputs)


def update_observation_members(observation):
    """Add members to Observation from all its Planes.

    :param observation Observation instance to add members to
    """
    members_inputs = TypedSet(ObservationURI,)
    for plane in observation.planes.values():
        if (plane.provenance is not None and
                plane.provenance.inputs is not None):
            for inpt in plane.provenance.inputs:
                members_inputs.add(inpt.get_observation_uri())
                logging.debug('Adding Observation URI {}'.format(
                    inpt.get_observation_uri()))
    mc.update_typed_set(observation.members, members_inputs)
