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

import logging
import sys

from . import caom2blueprint
import traceback

APP_NAME = 'fits2caom2'


class ConvertFromJava:
    """
    Do the work that makes the input from a Java fits2caom2 run usable by the
    ObsBlueprint class in this python implementation.
    """

    def __init__(self, blueprint, user_supplied_config):
        # invert the dict for a quick lookup of keywords referenced by the
        # plan, because the blueprint relies on the values, not the keys
        self._inverse_plan = {}
        for key, value in blueprint._plan.items():
            if isinstance(value, tuple):
                for ii in value[0]:
                    if ii in self._inverse_plan:
                        self._inverse_plan[ii].append(key)
                    else:
                        self._inverse_plan[ii] = [key]

        # invert the dict for a quick lookup of config reference values,
        # because the blueprint relies on the values, not the keys
        self._inverse_user_supplied_config = {}
        if user_supplied_config:
            for k, v in user_supplied_config.items():
                if v in self._inverse_user_supplied_config:
                    self._inverse_user_supplied_config[v].append(k)
                else:
                    self._inverse_user_supplied_config[v] = [k]

    def get_caom2_elements(self, lookup):
        if lookup in caom2blueprint.ObsBlueprint._CAOM2_ELEMENTS:
            return [lookup]
        elif lookup in self._inverse_user_supplied_config.keys():
            return self._inverse_user_supplied_config[lookup]
        elif lookup in self._inverse_plan.keys():
            return self._inverse_plan[lookup]
        else:
            raise ValueError(
                '{} caom2 element not found in the plan (spelling?).'.
                format(lookup))


# Mimic the default java fits2caom2.config file content, to support the
# indirection from named config values to named defaults and overrides.
#
# Drop-in use seems to expect that the default config file exists, and
# therefore certain indirections also exist.
#
# This is for drop-in functionality support only, and should not be relied on
# going forward.
_JAVA_CAOM2_CONFIG = {
    'DerivedObservation.members': 'members',

    'Observation.type': 'OBSTYPE',
    'Observation.intent': 'obs.intent',
    'Observation.sequenceNumber': 'obs.sequenceNumber',
    'Observation.metaRelease': 'obs.metaRelease',

    'Observation.algorithm.name': 'algorithm.name',

    'Observation.instrument.name': 'instrument.name',
    'Observation.instrument.keywords': 'instrument.keywords',

    'Observation.proposal.id': 'proposal.id',
    'Observation.proposal.pi': 'proposal.pi',
    'Observation.proposal.project': 'proposal.project',
    'Observation.proposal.title': 'proposal.title',
    'Observation.proposal.keywords': 'proposal.keywords',

    'Observation.target.name': 'target.name',
    'Observation.target.type': 'target.type',
    'Observation.target.standard': 'target.standard',
    'Observation.target.redshift': 'target.redshift',
    'Observation.target.keywords': 'target.keywords',
    'Observation.target.moving': 'target.moving',

    'Observation.telescope.name': 'telescope.name',
    'Observation.telescope.geoLocationX': 'telescope.geoLocationX',
    'Observation.telescope.geoLocationY': 'telescope.geoLocationY',
    'Observation.telescope.geoLocationZ': 'telescope.geoLocationZ',
    'Observation.telescope.keywords': 'telescope.keywords',

    'Observation.environment.seeing': 'environment.seeing',
    'Observation.environment.humidity': 'environment.humidity',
    'Observation.environment.elevation': 'environment.elevation',
    'Observation.environment.tau': 'environment.tau',
    'Observation.environment.wavelengthTau': 'environment.wavelengthTau',
    'Observation.environment.ambientTemp': 'environment.ambientTemp',
    'Observation.environment.photometric': 'environment.photometric',

    'Plane.metaRelease': 'plane.metaRelease',
    'Plane.dataRelease': 'plane.dataRelease',
    'Plane.dataProductType': 'plane.dataProductType',
    'Plane.calibrationLevel': 'plane.calibrationLevel',

    'Plane.provenance.name': 'provenance.name',
    'Plane.provenance.version': 'provenance.version',
    'Plane.provenance.project': 'provenance.project',
    'Plane.provenance.producer': 'provenance.producer',
    'Plane.provenance.runID': 'provenance.runID',
    'Plane.provenance.reference': 'provenance.reference',
    'Plane.provenance.lastExecuted': 'provenance.lastExecuted',
    'Plane.provenance.keywords': 'provenance.keywords',
    'Plane.provenance.inputs': 'provenance.inputs',

    'Plane.metrics.sourceNumberDensity': 'metrics.sourceNumberDensity',
    'Plane.metrics.background': 'metrics.background',
    'Plane.metrics.backgroundStddev': 'metrics.backgroundStddev',
    'Plane.metrics.fluxDensityLimit': 'metrics.fluxDensityLimit',
    'Plane.metrics.magLimit': 'metrics.magLimit',

    'Artifact.productType': 'artifact.productType',
    'Artifact.releaseType': 'artifact.releaseType',

    'Part.name': 'part.name',
    'Part.productType': 'part.productType',

    'Chunk.naxis': 'ZNAXIS,NAXIS',
    'Chunk.observableAxis': 'chunk.observableAxis',
    'Chunk.positionAxis1': 'getPositionAxis()',
    'Chunk.positionAxis2': 'getPositionAxis()',
    'Chunk.energyAxis': 'getEnergyAxis()',
    'Chunk.timeAxis': 'getTimeAxis()',
    'Chunk.polarizationAxis': 'getPolarizationAxis()',

    'Chunk.observable.dependent.bin': 'observable.dependent.bin',
    'Chunk.observable.dependent.axis.ctype': 'observable.dependent.ctype',
    'Chunk.observable.dependent.axis.cunit': 'observable.dependent.cunit',
    'Chunk.observable.independent.bin': 'observable.independent.bin',
    'Chunk.observable.independent.axis.ctype': 'observable.independent.ctype',
    'Chunk.observable.independent.axis.cunit': 'observable.independent.cunit',

    'Chunk.position.coordsys': 'RADECSYS,RADESYS',
    'Chunk.position.equinox': 'EQUINOX,EPOCH',
    'Chunk.position.resolution': 'position.resolution',
    'Chunk.position.axis.axis1.ctype': 'CTYPE{positionAxis1}',
    'Chunk.position.axis.axis1.cunit': 'CUNIT{positionAxis1}',
    'Chunk.position.axis.axis2.ctype': 'CTYPE{positionAxis2}',
    'Chunk.position.axis.axis2.cunit': 'CUNIT{positionAxis2}',
    'Chunk.position.axis.error1.syser': 'CSYER{positionAxis1}',
    'Chunk.position.axis.error1.rnder': 'CRDER{positionAxis1}',
    'Chunk.position.axis.error2.syser': 'CSYER{positionAxis2}',
    'Chunk.position.axis.error2.rnder': 'CRDER{positionAxis2}',
    'Chunk.position.axis.function.cd11': 'CD{positionAxis1}_{positionAxis1}',
    'Chunk.position.axis.function.cd12': 'CD{positionAxis1}_{positionAxis2}',
    'Chunk.position.axis.function.cd21': 'CD{positionAxis2}_{positionAxis1}',
    'Chunk.position.axis.function.cd22': 'CD{positionAxis2}_{positionAxis2}',
    'Chunk.position.axis.function.dimension.naxis1':
        'ZNAXIS{positionAxis1},NAXIS{positionAxis1}',
    'Chunk.position.axis.function.dimension.naxis2':
        'ZNAXIS{positionAxis2},NAXIS{positionAxis2}',
    'Chunk.position.axis.function.refCoord.coord1.pix': 'CRPIX{positionAxis1}',
    'Chunk.position.axis.function.refCoord.coord1.val': 'CRVAL{positionAxis1}',
    'Chunk.position.axis.function.refCoord.coord2.pix': 'CRPIX{positionAxis2}',
    'Chunk.position.axis.function.refCoord.coord2.val': 'CRVAL{positionAxis2}',
    'Chunk.position.axis.range.start.coord1.pix':
        'position.range.start.coord1.pix',
    'Chunk.position.axis.range.start.coord1.val':
        'position.range.start.coord1.val',
    'Chunk.position.axis.range.start.coord2.pix':
        'position.range.start.coord2.pix',
    'Chunk.position.axis.range.start.coord2.val':
        'position.range.start.coord2.val',
    'Chunk.position.axis.range.end.coord1.pix':
        'position.range.end.coord1.pix',
    'Chunk.position.axis.range.end.coord1.val':
        'position.range.end.coord1.val',
    'Chunk.position.axis.range.end.coord2.pix':
        'position.range.end.coord2.pix',
    'Chunk.position.axis.range.end.coord2.val':
        'position.range.end.coord2.val',

    'Chunk.energy.specsys': 'SPECSYS',
    'Chunk.energy.ssysobs': 'SSYSOBS',
    'Chunk.energy.restfrq': 'RESTFRQ',
    'Chunk.energy.restwav': 'RESTWAV',
    'Chunk.energy.velosys': 'VELOSYS',
    'Chunk.energy.zsource': 'ZSOURCE',
    'Chunk.energy.ssyssrc': 'SSYSSRC',
    'Chunk.energy.velang': 'VELANG',
    'Chunk.energy.bandpassName': 'bandpassName',
    'Chunk.energy.resolvingPower': 'resolvingPower',
    'Chunk.energy.transition.species': 'energy.transition.species',
    'Chunk.energy.transition.transition': 'energy.transition.transition',
    'Chunk.energy.axis.axis.ctype': 'CTYPE{energyAxis}',
    'Chunk.energy.axis.axis.cunit': 'CUNIT{energyAxis}',
    'Chunk.energy.axis.bounds.samples': 'energy.samples',
    'Chunk.energy.axis.error.syser': 'CSYER{energyAxis}',
    'Chunk.energy.axis.error.rnder': 'CRDER{energyAxis}',
    'Chunk.energy.axis.function.naxis': 'NAXIS{energyAxis}',
    'Chunk.energy.axis.function.delta': 'CDELT{energyAxis}',
    'Chunk.energy.axis.function.refCoord.pix': 'CRPIX{energyAxis}',
    'Chunk.energy.axis.function.refCoord.val': 'CRVAL{energyAxis}',
    'Chunk.energy.axis.range.start.pix': 'energy.range.start.pix',
    'Chunk.energy.axis.range.start.val': 'energy.range.start.val',
    'Chunk.energy.axis.range.end.pix': 'energy.range.end.pix',
    'Chunk.energy.axis.range.end.val': 'energy.range.end.val',

    'Chunk.polarization.axis.axis.ctype': 'CTYPE{polarizationAxis}',
    'Chunk.polarization.axis.axis.cunit': 'CUNIT{polarizationAxis}',
    'Chunk.polarization.axis.bounds.samples': 'polarization.samples',
    'Chunk.polarization.axis.error.syser': 'polarization.error.syser',
    'Chunk.polarization.axis.error.rnder': 'polarization.error.reder',
    'Chunk.polarization.axis.function.naxis': 'NAXIS{polarizationAxis}',
    'Chunk.polarization.axis.function.delta': 'CDELT{polarizationAxis}',
    'Chunk.polarization.axis.function.refCoord.pix': 'CRPIX{polarizationAxis}',
    'Chunk.polarization.axis.function.refCoord.val': 'CRVAL{polarizationAxis}',
    'Chunk.polarization.axis.range.start.pix': 'polarization.range.start.pix',
    'Chunk.polarization.axis.range.start.val': 'polarization.range.start.val',
    'Chunk.polarization.axis.range.end.pix': 'polarization.range.end.pix',
    'Chunk.polarization.axis.range.end.val': 'polarization.range.end.val',

    'Chunk.time.exposure': 'time.exposure',
    'Chunk.time.resolution': 'time.resolution',
}


def apply_java_config(file_name, use_only_defaults=False):
    """
    Override CONFIG with externally-supplied values.

    The override file can contain information for more than one input file,
    as well as providing information for different HDUs.

    :param file_name Name of the configuration file to load.
    :param use_only_defaults if True, rely on _JAVA_CAOM2_CONFIG content for
        config information.
    :return: dict representation of file content.
    """
    if use_only_defaults:
        d = _JAVA_CAOM2_CONFIG
    else:
        d = load_config(file_name)
        copy = _JAVA_CAOM2_CONFIG.copy()
        copy.update(d)
        d = copy
    return d


def load_config(file_name):
    """
    Override CONFIG with externally-supplied values.

    The override file can contain information for more than one input file,
    as well as providing information for different HDUs.

    :param file_name Name of the configuration file to load.
    :return: dict representation of file content.
    """
    d = {}
    extension = 0
    artifact = ''
    ptr = d
    with open(file_name) as file:
        for line in file:
            if line.startswith('?'):
                # file, extension-specific content
                if line.find('#[') == -1:
                    artifact = line.split('?')[1].strip()
                    extension = 0
                    if 'artifacts' not in d.keys():
                        d['artifacts'] = {}
                else:
                    artifact = line.split('?')[1].split('#[')[0].strip()
                    extension = \
                        int(line.split('#[')[1].split(']')[0].strip())
                    logging.debug(
                        'Adding overrides for artifact {} in extension {}'.
                        format(artifact, extension))
                if artifact not in d['artifacts'].keys():
                    d['artifacts'][artifact] = {}
                if extension not in d['artifacts'][artifact].keys():
                    d['artifacts'][artifact][extension] = {}
                ptr = d['artifacts'][artifact][extension]
            elif line.find('=') == -1 or line.startswith('#'):
                continue
            else:
                key, value = line.split('=', 1)
                ptr[key.strip()] = value.strip()
    return d


def _update_axis_info(parser, defaults, overrides, config):
    # look for info regarding axis types in the default and override file
    if config is None:
        raise ValueError('Empty config when updating axis info.')
    energy_axis = None
    polarization_axis = None
    time_axis = None
    ra_axis = None
    dec_axis = None
    obs_axis = None
    for i in defaults, overrides:
        for key, value in i.items():
            if (key.startswith('CTYPE')) and key[-1].isdigit():
                value = value.split('-')[0]
                if value in caom2blueprint.ENERGY_CTYPES:
                    energy_axis = key[-1]
                elif value in caom2blueprint.POLARIZATION_CTYPES:
                    polarization_axis = key[-1]
                elif value in caom2blueprint.TIME_KEYWORDS:
                    time_axis = key[-1]
                elif value in caom2blueprint.POSITION_CTYPES[0]:
                    ra_axis = key[-1]
                elif value in caom2blueprint.POSITION_CTYPES[1]:
                    dec_axis = key[-1]
                elif value in caom2blueprint.OBSERVABLE_CTYPES:
                    obs_axis = key[-1]
                else:
                    raise ValueError(f'Unrecognized CTYPE: {value}')

    ignore = '{ignore}'
    if ('Chunk.position' not in config) or \
            (config['Chunk.position'] != ignore):
        if ra_axis and dec_axis:
            parser.configure_position_axes((ra_axis, dec_axis))
        elif ra_axis or dec_axis:
            raise ValueError('Only one positional axis found (ra/dec): {}/{}'.
                             format(ra_axis, dec_axis))
        else:
            # assume that positional axis are 1 and 2 by default
            if time_axis in ['1', '2'] or energy_axis in ['1', '2'] or \
               polarization_axis in ['1', '2'] or obs_axis in ['1', '2']:
                raise ValueError('Cannot determine the positional axis')
            else:
                parser.configure_position_axes(('1', '2'))

    if time_axis and (('Chunk.time' not in config) or
       (config['Chunk.time'] != ignore)):
        parser.configure_time_axis(time_axis)

    if energy_axis and (('Chunk.energy' not in config) or
       (config['Chunk.energy'] != ignore)):
        parser.configure_energy_axis(energy_axis)

    if polarization_axis and (('Chunk.polarization' not in config) or
       (config['Chunk.polarization'] != ignore)):
        parser.configure_polarization_axis(polarization_axis)

    if obs_axis and (('Chunk.observable' not in config) or
       (config['Chunk.observable'] != ignore)):
        parser.configure_observable_axis(obs_axis)


def update_blueprint(obs_blueprint, artifact_uri=None, config=None,
                     defaults=None, overrides=None):
    """
    Update an observation blueprint according to defaults and/or overrides as
    configured by the user.
    :param obs_blueprint: ObsBlueprint to update
    :param artifact_uri: Where the overrides come from, and where
    to apply them.
    :param config: Input configuration in a dict.
    :param defaults: FITS header and configuration default values in a dict.
    :param overrides: FITS header keyword and configuration default overrides
    in a dict.
    :return: String containing error messages. Result is None if no errors
    encountered.
    """

    _update_axis_info(obs_blueprint, defaults, overrides, config)

    convert = ConvertFromJava(obs_blueprint, config)
    errors = []
    if config:
        logging.debug(
            f'Setting user-supplied configuration for {artifact_uri}.')
        for key, value in config.items():
            try:
                if value.isupper() and value.find('.') == -1:
                    # assume FITS keywords, in the 0th extension,
                    # and add them to the blueprint
                    for caom2_key in convert.get_caom2_elements(key):
                        obs_blueprint.add_attribute(caom2_key, value)
            except ValueError:
                errors.append(f'{key}: {sys.exc_info()[1]}')
        logging.debug(
            f'User-supplied configuration applied for {artifact_uri}.')

    if defaults:
        logging.debug(f'Setting defaults for {artifact_uri}')
        for key, value in defaults.items():
            try:
                for caom2_key in convert.get_caom2_elements(key):
                    obs_blueprint.set_default(caom2_key, value)
                    logging.debug(
                        '{} setting default value to {}'.format(
                            caom2_key, value))
            except ValueError:
                errors.append(f'{key}: {sys.exc_info()[1]}')
        logging.debug(f'Defaults set for {artifact_uri}.')

    if overrides:
        logging.debug(f'Setting overrides for {artifact_uri}.')
        for key, value in overrides.items():
            if key == 'BITPIX':
                logging.debug(f'01/11/18 Chris said ignore {key!r}.')
                continue
            if key == 'artifacts' and artifact_uri in overrides['artifacts']:
                logging.debug('Found extension overrides for URI {}.'.format(
                    artifact_uri))
                for extension in overrides['artifacts'][artifact_uri].keys():
                    for ext_key, ext_value in \
                      overrides['artifacts'][artifact_uri][extension].items():
                        if ext_key == 'BITPIX':
                            logging.debug(
                                f'01/11/18 Chris said ignore {key!r}.')
                            continue
                        try:
                            for caom2_key in \
                                    convert.get_caom2_elements(ext_key):
                                obs_blueprint.set(caom2_key, ext_value,
                                                  extension)
                                logging.debug(('{} set override value to {} '
                                               'in extension {}.').format(
                                    caom2_key, ext_value, extension))
                        except ValueError:
                            errors.append('{}: ext {} {}'.format(
                                key, extension, sys.exc_info()[1]))
            else:
                try:
                    for caom2_key in convert.get_caom2_elements(key):
                        obs_blueprint.set(caom2_key, value)
                except ValueError:
                    errors.append(f'{key}: {sys.exc_info()[1]}')
        logging.debug(f'Overrides set for {artifact_uri}.')

        if errors:
            return '\n'.join(errors)
        else:
            return None


def main_app():
    parser = caom2blueprint.get_arg_parser()

    # add legacy fits2caom2 arguments
    parser.add_argument('--config', required=False,
                        help=('optional CAOM2 utype to keyword config file to '
                              'merge with the internal configuration'))

    parser.add_argument('--default',
                        help='file with default values for keywords')
    parser.add_argument('--override',
                        help='file with override values for keywords')

    if len(sys.argv) < 2:
        # correct error message when running python3
        parser.print_usage(file=sys.stderr)
        sys.stderr.write(f"{APP_NAME}: error: too few arguments\n")
        sys.exit(-1)

    args = parser.parse_args()

    config = None
    if args.config:
        config = apply_java_config(args.config)
        logging.debug(f'Apply configuration from {args.config}.')

    defaults = {}
    if args.default:
        defaults = load_config(args.default)
        logging.debug(f'Apply defaults from {args.default}.')

    overrides = {}
    if args.override:
        overrides = load_config(args.override)
        logging.debug(f'Apply overrides from {args.override}.')

    obs_blueprint = {}
    for i, uri in enumerate(args.fileURI):
        if '.h5' in uri:
            obs_blueprint[uri] = caom2blueprint.Hdf5ObsBlueprint()
        else:
            obs_blueprint[uri] = caom2blueprint.ObsBlueprint()
        if config:
            result = update_blueprint(obs_blueprint[uri], uri,
                                      config, defaults, overrides)
            if result:
                logging.debug(
                    f'Errors parsing the config files: {result}')

    try:
        caom2blueprint.proc(args, obs_blueprint)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)

    logging.info("DONE")
