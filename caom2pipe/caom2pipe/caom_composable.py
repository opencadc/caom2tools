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
import os

from caom2 import TypedSet, ObservationURI, PlaneURI, Chunk, CoordPolygon2D
from caom2 import ValueCoord2D

from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc

__all__ = ['exec_footprintfinder', 'update_plane_provenance',
           'update_observation_members', 'reset_energy', 'reset_position',
           'reset_observable']


def exec_footprintfinder(chunk, science_fqn, log_file_directory, obs_id,
                         params='-f'):
    """Execute the footprintfinder on a file. All preconditions for successful
    execution should be in place i.e. the file exists, and is unzipped (because
    that is faster).

    :param chunk The CAOM Chunk that will have Position Bounds information
        added
    :param science_fqn A string of the fully-qualified file name for
        footprintfinder to run on
    :param log_file_directory A string of the fully-qualified name for the log
        directory, where footprintfinder output files will be moved to, after
        execution
    :param obs_id specifies location where footprintfinder log files end up
    :param params specific footprintfinder parameters by collection - default
        forces full-chip, regardless of illumination
    """
    logging.debug('Begin _update_position')
    mc.check_param(chunk, Chunk)

    # local import because footprintfinder depends on matplotlib being
    # installed, which is not declared as a caom2pipe dependency
    import footprintfinder

    if (chunk.position is not None
            and chunk.position.axis is not None):
        logging.debug('position exists, calculate footprints for {}.'.format(
            science_fqn))
        full_area, footprint_xc, footprint_yc, ra_bary, dec_bary, \
            footprintstring, stc = footprintfinder.main(
                '-r {} {}'.format(params, science_fqn))
        logging.debug('footprintfinder result: full area {} '
                      'footprint xc {} footprint yc {} ra bary {} '
                      'dec_bary {} footprintstring {} stc {}'.format(
                          full_area, footprint_xc, footprint_yc, ra_bary,
                          dec_bary, footprintstring, stc))
        bounds = CoordPolygon2D()
        coords = None
        fp_results = stc.split('Polygon FK5')
        if len(fp_results) > 1:
            coords = fp_results[1].split()
        else:
            fp_results = stc.split('Polygon ICRS')
            if len(fp_results) > 1:
                coords = fp_results[1].split()

        if coords is None:
            raise mc.CadcException(
                'Do not recognize footprint {}'.format(stc))

        index = 0
        while index < len(coords):
            vertex = ValueCoord2D(mc.to_float(coords[index]),
                                  mc.to_float(coords[index + 1]))
            bounds.vertices.append(vertex)
            index += 2
            logging.debug('Adding vertex\n{}'.format(vertex))
        chunk.position.axis.bounds = bounds

        return_file = '{}_footprint.txt'.format(obs_id)
        return_string_file = '{}_footprint_returnstring.txt'.format(obs_id)
        _handle_footprint_logs(log_file_directory, return_file)
        _handle_footprint_logs(log_file_directory, return_string_file)

    else:
        logging.info('No position information for footprint generation.')
    logging.debug('Done _update_position.')


def _handle_footprint_logs(log_file_directory, log_file):
    """Move footprintfinder logs to specific log directory, if there
    is one."""
    orig_log_fqn = os.path.join(os.getcwd(), log_file)
    if log_file_directory is not None and os.path.exists(log_file_directory):
        if os.path.exists(orig_log_fqn):
            log_fqn = os.path.join(log_file_directory, log_file)
            os.rename(orig_log_fqn, log_fqn)
            logging.debug('Moving footprint log file from {} to {}'.format(
                orig_log_fqn, log_fqn))
    else:
        logging.debug('Removing footprint log file {}'.format(orig_log_fqn))
        os.unlink(orig_log_fqn)


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


def reset_energy(chunk):
    """
    :param chunk: Set the energy component of a chunk to None as a side-effect.
    """
    chunk.energy = None
    chunk.energy_axis = None


def reset_position(chunk):
    """
    :param chunk: Set the position component of a chunk to None as a
    side-effect.
    """
    chunk.position = None
    chunk.position_axis_1 = None
    chunk.position_axis_2 = None


def reset_observable(chunk):
    """
    :param chunk: Set the observable component of a chunk to None as a
    side-effect.
    """
    chunk.observable = None
    chunk.observable_axis = None
