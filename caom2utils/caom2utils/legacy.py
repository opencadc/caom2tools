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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from builtins import str

import logging
import sys

from . import fits2caom2

APP_NAME = 'fits2caom2'

__all__ = ['main_app', 'update_blueprint']


class ConvertFromJava(object):
    """
    Do the work that makes the input from a Java fits2caom2 run usable by the
    ObsBlueprint class in this python implementation.
    """

    def __init__(self, blueprint, user_supplied_config):
        # for a quick lookup of keywords referenced by the plan
        self._inverse_plan = {}
        for key, value in blueprint._plan.items():
            if isinstance(value, tuple):
                for ii in value[0]:
                    if ii in self._inverse_plan:
                        self._inverse_plan[ii].append(key)
                    else:
                        self._inverse_plan[ii] = [key]

        # for a quick lookup of config reference values
        self._inverse_user_supplied_config = {}
        if user_supplied_config:
            for k, v in user_supplied_config.items():
                if v in self._inverse_user_supplied_config:
                    self._inverse_user_supplied_config[v].append(k)
                else:
                    self._inverse_user_supplied_config[v] = [k]

    def get_caom2_elements(self, lookup):
        if lookup in fits2caom2.ObsBlueprint._CAOM2_ELEMENTS:
            return [lookup]
        elif lookup in self._inverse_user_supplied_config.keys():
            return self._inverse_user_supplied_config[lookup]
        elif lookup in self._inverse_plan.keys():
            return self._inverse_plan[lookup]
        else:
            raise ValueError(
                '{} caom2 element not found in the plan (spelling?).'.
                format(lookup))


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
    assert config is not None
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
                if value in fits2caom2.ENERGY_CTYPES:
                    energy_axis = key[-1]
                elif value in fits2caom2.POLARIZATION_CTYPES:
                    polarization_axis = key[-1]
                elif value in fits2caom2.TIME_KEYWORDS:
                    time_axis = key[-1]
                elif value in fits2caom2.POSITION_CTYPES[0]:
                    ra_axis = key[-1]
                elif value in fits2caom2.POSITION_CTYPES[1]:
                    dec_axis = key[-1]
                elif value in fits2caom2.OBSERVABLE_CTYPES:
                    obs_axis = key[-1]
                else:
                    raise ValueError('Unrecognized CTYPE: {}'.format(value))

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
            'Setting user-supplied configuration for {}.'.format(artifact_uri))
        for key, value in config.items():
            try:

                if value.isupper() and value.find('.') == -1:
                    # assume FITS keywords, in the 0th extension,
                    # and add them to the blueprint
                    for caom2_key in convert.get_caom2_elements(key):
                        obs_blueprint.set_fits_attribute(caom2_key, [value])
            except ValueError:
                errors.append(('{}: {}'.format(key, sys.exc_info()[1])))
        logging.debug(
            'User-supplied configuration applied for {}.'.format(artifact_uri))

    if defaults:
        logging.debug('Setting defaults for {}'.format(artifact_uri))
        for key, value in defaults.items():
            try:
                for caom2_key in convert.get_caom2_elements(key):
                    obs_blueprint.set_default(caom2_key, value)
                    logging.debug(
                        '{} setting default value to {}'.format(
                            caom2_key, value))
            except ValueError:
                errors.append('{}: {}'.format(key, sys.exc_info()[1]))
        logging.debug('Defaults set for {}.'.format(artifact_uri))

    if overrides:
        logging.debug('Setting overrides for {}.'.format(artifact_uri))
        for key, value in overrides.items():
            if key == 'BITPIX':
                logging.debug('01/11/18 Chris said ignore {!r}.'.format(key))
                continue
            if key == 'artifacts' and artifact_uri in overrides['artifacts']:
                logging.debug('Found extension overrides for URI {}.'.format(
                    artifact_uri))
                for extension in overrides['artifacts'][artifact_uri].keys():
                    for ext_key, ext_value in \
                      overrides['artifacts'][artifact_uri][extension].items():
                        if ext_key == 'BITPIX':
                            logging.debug(
                                '01/11/18 Chris said ignore {!r}.'.format(key))
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
                    errors.append('{}: {}'.format(key, sys.exc_info()[1]))
        logging.debug('Overrides set for {}.'.format(artifact_uri))

        if errors:
            return '\n'.join(errors)
        else:
            return None


def _dump_config(parser, uri):
    f = None
    try:
        temp = uri.split('/')
        mod_uri = temp[len(temp) - 1]
        fname = './{}.mod.fits'.format(mod_uri)
        logging.debug('Writing modified fits file to {}.'.format(fname))
        f = open(fname, 'w')
        for index, extension in enumerate(parser._headers):
            f.write('\nHeader {}\n'.format(index))
            f.write(extension.tostring('\n'))
        f.close()
        fname = './{}.blueprint.out'.format(mod_uri)
        logging.debug('Writing blueprint to {}.'.format(fname))
        f = open(fname, 'w')
        f.write(str(parser.blueprint))
        f.close()
        fname = './{}.errors.out'.format(mod_uri)
        logging.debug('Writing errors to {}.'.format(fname))
        f = open(fname, 'w')
        for ii in parser._errors:
            f.write(ii)
            f.write('\n')
        f.close()
    except EnvironmentError:
        logging.warning('Failed to dump config. {}'.format(sys.exc_info()[1]))
    finally:
        if f:
            f.close()


def main_app():
    parser = fits2caom2.get_arg_parser()

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
        sys.stderr.write("{}: error: too few arguments\n".format(APP_NAME))
        sys.exit(-1)

    args = parser.parse_args()

    config = None
    if args.config:
        config = load_config(args.config)
        logging.debug('Apply configuration from {}.'.format(args.config))

    defaults = {}
    if args.default:
        defaults = load_config(args.default)
        logging.debug('Apply defaults from {}.'.format(args.default))

    overrides = {}
    if args.override:
        overrides = load_config(args.override)
        logging.debug('Apply overrides from {}.'.format(args.override))

    obs_blueprint = {}
    for i, uri in enumerate(args.fileURI):
        obs_blueprint[uri] = fits2caom2.ObsBlueprint()
        if config:
            result = update_blueprint(obs_blueprint[uri], uri,
                                      config, defaults, overrides)
            if result:
                logging.warning(
                    'Errors parsing the config files: {}'.format(result))

    try:
        fits2caom2.proc(args, obs_blueprint)
    except Exception as e:
        logging.error(e)
        sys.exit(-1)

    logging.info("DONE")
