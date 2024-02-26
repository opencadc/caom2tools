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

"""
There is a dual inheritance hierarchy in this module:


         BlueprintParser
               ^
               |
          ContentParser <1:n>--------------------------- WcsParser
               ^                                            ^
               |                                            |
        ---------------------                    ------------------------
        |                   |                    |                      |
        |               Hdf5Parser <1:n>----- Hdf5WcsParser             |
        |                                                               |
     FitsParser <1:n>--------------------------------------------- FitsWcsParser

The *WcsParser hierarchy uses astropy.wcs for WCS construction and correctness when building CAOM records.

When building CAOM records, use:
- the BlueprintParser for records with no WCS information,
- the ContentParser for populating records with WCS information from database queries, etc,
- the FitsParser for populating records with WCS information from FITS files, and
- the Hdf5Parser for populating with WCS information from HDF5 files.

"""
import argparse
from logging.handlers import TimedRotatingFileHandler

from cadcutils import version
from caom2 import (
    Artifact,
    Algorithm,
    ChecksumURI,
    CompositeObservation,
    DerivedObservation,
    ObservationReader,
    ObservationWriter,
    Plane,
    ProductType,
    ReleaseType,
    SimpleObservation,
)
from caom2utils import data_util
from caom2utils.caomvalidator import validate
from caom2utils.wcsvalidator import InvalidWCSError
from caom2utils.blueprints import Hdf5ObsBlueprint, ObsBlueprint, _to_int, _to_str
from caom2utils.parsers import BlueprintParser, ContentParser, FitsParser, Hdf5Parser
import importlib
import logging
import os
import requests
import sys
import tempfile
import traceback
from urllib.parse import urlparse
from cadcutils import net, util
from cadcdata import FileInfo
from vos import Client
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

APP_NAME = 'caom2gen'

__all__ = [
    'augment',
    'DispatchingFormatter',
    'gen_proc',
    'get_arg_parser',
    'get_external_headers',
    'get_gen_proc_arg_parser',
    'get_vos_headers',
    'proc',
    'update_artifact_meta',
]


GLOBAL_STORAGE_RESOURCE_ID = "ivo://cadc.nrc.ca/global/raven"


class DispatchingFormatter:
    """Dispatch formatter for logger and it's sub-logger, so there can be multiple formatters."""

    def __init__(self, formatters, default_formatter):
        self._formatters = formatters
        self._default_formatter = default_formatter

    def format(self, record):
        logger = logging.getLogger(record.name)
        while logger:
            # check if suitable formatter for current logger exists
            if logger.name in self._formatters:
                formatter = self._formatters[logger.name]
                break
            else:
                logger = logger.parent
        else:
            # if no formatter found, just use the default
            formatter = self._default_formatter
        return formatter.format(record)


def get_external_headers(external_url):
    try:
        session = requests.Session()
        retries = 10
        retry = Retry(total=retries, read=retries, connect=retries, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        r = session.get(external_url, timeout=20)
        if r.status_code == requests.codes.ok:
            headers = data_util.make_headers_from_string(r.text)
        else:
            headers = None
            logging.warning('Error {} when retrieving {} headers.'.format(r.status_code, external_url))
        r.close()
        return headers
    except Exception as e:
        logging.error(f'Connection failed to {external_url}.\n{e}')
        raise RuntimeError(e)


def get_vos_headers(uri, subject=None):
    """
    Creates the FITS headers object from a vospace file.
    :param uri: vos URI
    :param subject: user credentials. Anonymous if subject is None
    :return: List of headers corresponding to each extension. Each header is of astropy.wcs.Header type - essentially
      a dictionary of FITS keywords.
    """
    if uri.startswith('vos'):
        if subject is not None and subject.certificate is not None:
            client = Client(subject.certificate)
        else:
            client = Client()

        temp_filename = tempfile.NamedTemporaryFile()
        client.copy(uri, temp_filename.name, head=True)
        return data_util.get_local_file_headers(f'file://{temp_filename.name}')
    else:
        # this should be a programming error by now
        raise NotImplementedError('Only vos type URIs supported')


def _get_and_update_artifact_meta(uri, artifact, subject=None, connected=True, client=None):
    """
    Updates contentType, contentLength and contentChecksum of an artifact
    :param artifact:
    :param subject: User credentials
    :param connected: True if there's a network connection
    :param client: connection to CADC storage
    :return:
    """
    logging.debug(f'Begin _get_and_update_artifact_meta for {uri}')
    file_url = urlparse(uri)
    if file_url.scheme == 'gemini' and '.jpg' not in file_url.path:
        # will get file metadata from Gemini JSON summary for fits, because the metadata is available long before
        # the data will be stored at CADC
        return
    elif file_url.scheme == 'vos':
        metadata = _get_vos_meta(subject, uri)
    elif file_url.scheme == 'file':
        if file_url.path.endswith('.header') and subject is not None and connected:
            if artifact.uri.startswith('vos'):
                metadata = _get_vos_meta(subject, artifact.uri)
            else:
                # if header is on disk, get the content_* from CADC
                metadata = client.info(artifact.uri)
                if metadata is None:
                    logging.info('Could not find {} at CADC. No Artifact ' 'metadata.'.format(artifact.uri))
                    return
        else:
            metadata = data_util.get_local_file_info(file_url.path)
    else:
        metadata = client.info(uri)
        if metadata is None:
            logging.info('Could not find {} at CADC. No Artifact ' 'metadata.'.format(artifact.uri))
            return

    update_artifact_meta(artifact, metadata)


def update_artifact_meta(artifact, file_info):
    """
    Updates contentType, contentLength and contentChecksum of an artifact
    :param artifact:
    :param file_info
    :return:
    """
    logging.debug(
        'old artifact metadata - '
        'uri({}), encoding({}), size({}), type({})'.format(
            artifact.uri, artifact.content_checksum, artifact.content_length, artifact.content_type
        )
    )
    if file_info.md5sum is not None:
        if file_info.md5sum.startswith('md5:'):
            checksum = ChecksumURI(file_info.md5sum)
        else:
            checksum = ChecksumURI(f'md5:{file_info.md5sum}')
        artifact.content_checksum = checksum
    artifact.content_length = _to_int(file_info.size)
    artifact.content_type = _to_str(file_info.file_type)
    logging.debug(
        'updated artifact metadata - '
        'uri({}), encoding({}), size({}), type({})'.format(
            artifact.uri, artifact.content_checksum, artifact.content_length, artifact.content_type
        )
    )


def _get_vos_meta(subject, uri):
    """
    Gets contentType, contentLength and contentChecksum of a VOS artifact
    :param subject: user credentials
    :param uri:
    :return:
    """
    if subject is not None and subject.certificate is not None:
        client = Client(subject.certificate)
    else:
        client = Client()
    node = client.get_node(uri, limit=None, force=False)
    return FileInfo(
        id=uri, size=node.props['length'], md5sum=node.props['MD5'], file_type=data_util.get_file_type(uri)
    )


def _lookup_blueprint(blueprints, uri):
    """
    Blueprint handling may be one-per-observation, or one-per-URI. Find the correct one here.
    :param blueprints: The collection of blueprints provided by the user.
    :param uri: Which blueprint to look for
    :return: the blueprint to apply to Observation creation.
    """
    if len(blueprints) == 1:
        key, value = blueprints.popitem()
        return value
    else:
        return blueprints[uri]


def _extract_ids(cardinality):
    """
    Localize cardinality structure knowledge.

    :param cardinality:
    :return: product_id, artifact URI
    """
    return cardinality.split('/', 1)


def _augment(
    obs,
    product_id,
    uri,
    blueprint,
    subject,
    dumpconfig=False,
    validate_wcs=True,
    plugin=None,
    local=None,
    external_url=None,
    connected=True,
    use_blueprint_parser=False,
    client=None,
    **kwargs,
):
    """
    Find or construct a plane and an artifact to go with the observation under augmentation.

    :param obs: Observation - target of CAOM2 model augmentation
    :param product_id: Unique identifier for a plane in an Observation
    :param uri: Unique identifier for an artifact in a plane
    :param blueprint: Which blueprint to use when mapping from a telescope data model to CAOM2
    :param subject: authorization for any metdata access
    :param dumpconfig: print the blueprint to stdout
    :param validate_wcs: if true, call the validate method on the constructed observation, which checks that the WCS
        in the CAOM model is valid,
    :param plugin: what code to use for modifying a CAOM instance
    :param local: the input is the name of a file on disk
    :param external_url: if header information should be retrieved externally, this is where to find it
    :param client: StorageClientWrapper
    :return: an updated Observation
    """
    if dumpconfig:
        print(f'Blueprint for {uri}: {blueprint}')

    if product_id not in obs.planes.keys():
        obs.planes.add(Plane(product_id=str(product_id)))

    plane = obs.planes[product_id]

    if uri not in plane.artifacts.keys():
        plane.artifacts.add(Artifact(uri=str(uri), product_type=ProductType.SCIENCE, release_type=ReleaseType.DATA))

    meta_uri = uri
    visit_local = None
    if use_blueprint_parser:
        logging.debug(f'Using a BlueprintParser as requested for {uri}')
        parser = BlueprintParser(blueprint, uri=uri)
    elif local:
        if uri.startswith('vos'):
            if '.fits' in local or '.fits.gz' in local:
                meta_uri = f'file://{local}'
                logging.debug(f'Using a FitsParser for vos local {local}')
                headers = data_util.get_local_file_headers(local)
                parser = FitsParser(headers, blueprint, uri=uri)
            elif '.csv' in local:
                logging.debug(f'Using a BlueprintParser for vos local {local}')
                parser = BlueprintParser(blueprint, uri=uri)
            else:
                raise ValueError(f'Unexpected file type {local}')
        else:
            meta_uri = f'file://{local}'
            visit_local = local
            if '.header' in local or data_util.get_file_type(local) == 'application/fits':
                if uri.startswith('cadc'):
                    logging.debug(f'Using a FitsParser for local file {local}')
                    parser = FitsParser(local, blueprint, uri=uri)
                else:
                    logging.debug(f'Using a ContentParser for local file {local}')
                    parser = ContentParser(blueprint, uri)
            elif '.h5' in local:
                logging.debug(f'Using an Hdf5Parser for local file {local}')
                # h5py is an extra in this package since most collections do not require it
                import h5py

                temp = h5py.File(local)
                parser = Hdf5Parser(blueprint, uri, temp)
            else:
                # explicitly ignore headers for txt and image files
                logging.debug(f'Using a BlueprintParser for {local}')
                parser = BlueprintParser(blueprint, uri=uri)
    elif external_url:
        headers = get_external_headers(external_url)
        if headers is None:
            logging.debug('Using a BlueprintParser for un-retrievable remote headers ' '{}'.format(uri))
            parser = BlueprintParser(blueprint, uri=uri)
        else:
            logging.debug(f'Using a FitsParser for remote headers {uri}')
            parser = FitsParser(headers, blueprint, uri=uri)
    else:
        if '.fits' in uri:
            if uri.startswith('vos'):
                headers = get_vos_headers(uri, subject)
            elif uri.startswith('file'):
                headers = data_util.get_local_headers_from_fits(uri)
            else:
                headers = client.get_head(uri)
            logging.debug(f'Using a FitsParser for remote file {uri}')
            parser = FitsParser(headers, blueprint, uri=uri)
        else:
            # explicitly ignore headers for txt and image files
            logging.debug(f'Using a BlueprintParser for remote file {uri}')
            parser = BlueprintParser(blueprint, uri=uri)

    if parser is None:
        result = None
    else:
        _get_and_update_artifact_meta(meta_uri, plane.artifacts[uri], subject, connected, client)
        parser.augment_observation(observation=obs, artifact_uri=uri, product_id=plane.product_id)

        result = _visit(plugin, parser, obs, visit_local, product_id, uri, subject, **kwargs)

        if result is not None:
            if validate_wcs:
                try:
                    validate(obs)
                except InvalidWCSError as e:
                    logging.error(e)
                    tb = traceback.format_exc()
                    logging.debug(tb)
                    raise e

        if len(parser._errors) > 0:
            logging.debug('{} errors encountered while processing {!r}.'.format(len(parser._errors), uri))
            logging.debug(f'{parser._errors}')

    return result


def _load_module(module):
    """If a user provides code for execution during blueprint configuration, add that code to the execution
    environment of the interpreter here.

    :param module the fully-qualified path name to the source code from a user.
    """
    mname = os.path.basename(module)
    if '.' in mname:
        # remove extension from the provided name
        mname = mname.split('.')[0]
    pname = os.path.dirname(module)
    sys.path.append(pname)
    try:
        return importlib.import_module(mname)
    except ImportError as e:
        logging.debug(f'Looking for {mname!r} in {pname!r}')
        raise e


def caom2gen():
    parser = get_gen_proc_arg_parser()
    parser.add_argument(
        '--blueprint',
        nargs='+',
        required=True,
        help=(
            'list of files with blueprints for CAOM2 construction, in serialized format. If the list is of length 1, '
            'the same blueprint will be applied to all lineage entries. Otherwise, there must be a blueprint file '
            'per lineage entry.'
        ),
    )

    if len(sys.argv) < 2:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write(f'{APP_NAME}: error: too few arguments\n')
        sys.exit(-1)

    args = parser.parse_args()
    _set_logging(args.verbose, args.debug, args.quiet)

    module = None
    if args.module:
        module = _load_module(args.module)

    blueprints = {}
    if len(args.blueprint) == 1:
        # one blueprint to rule them all
        temp = ' '.join(ii for ii in args.lineage)
        if '.h5' in temp:
            blueprint = Hdf5ObsBlueprint(module=module)
        else:
            blueprint = ObsBlueprint(module=module)
        blueprint.load_from_file(args.blueprint[0])
        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            blueprints[uri] = blueprint
    else:
        # there needs to be the same number of blueprints as plane/artifact identifiers
        if len(args.lineage) != len(args.blueprint):
            logging.debug(f'Lineage: {args.lineage}')
            logging.debug(f'Blueprints: {args.blueprint}')
            sys.stderr.write(
                '{}: error: different number of blueprints '
                '{}  and files {}.'.format(APP_NAME, len(args.blueprint), len(args.lineage))
            )
            sys.exit(-1)

        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            logging.debug('Loading blueprint for {} from {}'.format(uri, args.blueprint[i]))
            if '.h5' in uri:
                blueprint = Hdf5ObsBlueprint(module=module)
            else:
                blueprint = ObsBlueprint(module=module)
            blueprint.load_from_file(args.blueprint[i])
            blueprints[uri] = blueprint

    try:
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed caom2gen execution.')
        logging.error(e)
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug(f'Done {APP_NAME} processing.')


def _gen_obs(obs_blueprints, in_obs_xml, collection=None, obs_id=None):
    """
    Determine whether to create a Simple or Derived Observation, or to read an existing Observation from an input
    file.

    :param obs_blueprints: Collection of blueprints provided to application.
    :param in_obs_xml: Existing observation information, contains the collection and obs_id values.
    :param collection: This plus the obs_id is a unique key for an observation.
    :param obs_id: This plus the collection is a unique key for an observation.
    :return: Initially constructed Observation.
    """
    obs = None
    if in_obs_xml:
        # append to existing observation
        reader = ObservationReader(validate=True)
        obs = reader.read(in_obs_xml)
    else:
        # determine the type of observation to create by looking for the the DerivedObservation.members in the
        # blueprints. If present in any of it assume derived
        for bp in obs_blueprints.values():
            if bp._get('DerivedObservation.members') is not None:
                logging.debug('Build a DerivedObservation')
                obs = DerivedObservation(
                    collection=collection, observation_id=obs_id, algorithm=Algorithm('composite')
                )
                break
            elif bp._get('CompositeObservation.members') is not None:
                logging.debug('Build a CompositeObservation with obs_id {}'.format(obs_id))
                obs = CompositeObservation(
                    collection=collection, observation_id=obs_id, algorithm=Algorithm('composite')
                )
                break
    if not obs:
        # build a simple observation
        logging.debug(f'Build a SimpleObservation with obs_id {obs_id}')
        obs = SimpleObservation(collection=collection, observation_id=obs_id, algorithm=Algorithm('exposure'))
    return obs


def _set_logging(verbose, debug, quiet):
    logger = logging.getLogger()
    # replace the StreamHandler with one that has custom formatters
    if logger.handlers:
        for handler in logger.handlers:
            if not isinstance(handler, TimedRotatingFileHandler):
                logger.removeHandler(handler)

    handler = logging.StreamHandler()
    handler.setFormatter(
        DispatchingFormatter(
            {
                'caom2utils.fits2caom2.FitsWcsParser': logging.Formatter(
                    '%(asctime)s:%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:' '%(lineno)d:%(message)s'
                ),
                'astropy': logging.Formatter(
                    '%(asctime)s:%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:' '%(lineno)d:%(message)s'
                ),
            },
            logging.Formatter('%(asctime)s:%(levelname)s:%(name)-12s:' '%(lineno)d:%(message)s'),
        )
    )
    logger.addHandler(handler)
    if verbose:
        logger.setLevel(logging.INFO)
        handler.setLevel(logging.INFO)
    elif debug:
        logger.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.ERROR)
        handler.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.WARN)
        handler.setLevel(logging.WARN)


def _get_common_arg_parser():
    """
    Returns the arg parser with common arguments between fits2caom2 and caom2gen
    :return: args parser
    """
    parser = util.get_base_parser(
        subparsers=False, version=version.version, default_resource_id=GLOBAL_STORAGE_RESOURCE_ID
    )

    parser.description = 'Augments an observation with information in one or more fits files.'

    parser.add_argument(
        '--dumpconfig', action='store_true', help=('output the utype to keyword mapping to ' 'the console')
    )

    parser.add_argument(
        '--not_connected',
        action='store_true',
        help=('if set, there is no internet connection, so ' 'skip service invocations.'),
    )

    parser.add_argument(
        '--no_validate',
        action='store_true',
        help=(
            'by default, the application will validate the WCS information for an observation. Specifying this flag '
            'skips that step.'
        ),
    )

    parser.add_argument(
        '-o', '--out', dest='out_obs_xml', help='output of augmented observation in XML', required=False
    )

    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument(
        '-i',
        '--in',
        dest='in_obs_xml',
        type=argparse.FileType('r'),
        help='input of observation to be augmented in XML',
    )
    in_group.add_argument(
        '--observation', nargs=2, help='observation in a collection', metavar=('collection', 'observationID')
    )
    parser.add_argument('--local', nargs='+', help=('list of files in local filesystem (same order ' 'as uri)'))
    return parser


def get_arg_parser():
    """
    Returns the arg parser with minimum arguments required to run fits2caom2
    :return: args parser
    """
    parser = _get_common_arg_parser()
    parser.add_argument('--productID', help='product ID of the plane in the observation', required=False)
    parser.add_argument('fileURI', help='URI of a fits file', nargs='+')
    return parser


def proc(args, obs_blueprints):
    """
    Function to process an observation according to command line arguments and a dictionary of blueprints.

    This implementation mirrors the Java implementation of fits2caom2, and the command line arguments it
    handles are productID and fileURI or local.

    There is no support for plugin execution to modify the blueprint with this access point.

    :param args: argparse args object containing the user supplied arguments. Arguments correspond to the parser
        returned by the get_arg_parser function
    :param obs_blueprints: dictionary of blueprints reguired to process the observation. The fileURIs represent the
        keys in this dictionary. Every fileURI in args.fileURI should have a corresponding blueprint.
    :return:
    """

    _set_logging(args.verbose, args.debug, args.quiet)

    if args.local and (len(args.local) != len(args.fileURI)):
        msg = ('number of local arguments not the same with file ' 'URIs ({} vs {})').format(
            len(args.local), args.fileURI
        )
        raise RuntimeError(msg)

    if args.in_obs_xml:
        obs = _gen_obs(obs_blueprints, args.in_obs_xml)
    else:
        obs = _gen_obs(obs_blueprints, None, args.observation[0], args.observation[1])

    if args.in_obs_xml and len(obs.planes) != 1:
        if not args.productID:
            msg = (
                'A productID parameter is required if there are zero or more than one planes in the input '
                'observation.',
            )
            raise RuntimeError(msg)

    subject = net.Subject.from_cmd_line_args(args)
    client = data_util.StorageClientWrapper(subject, resource_id=args.resource_id)
    validate_wcs = True
    if args.no_validate:
        validate_wcs = False

    for i, uri in enumerate(args.fileURI):
        blueprint = obs_blueprints[uri]
        # override the command-line argument for the plane product ID value
        product_id = blueprint._get('Plane.productID')
        if product_id is None or isinstance(product_id, tuple):
            if args.productID:
                product_id = args.productID
            else:
                msg = '{}{}'.format(
                    'A productID parameter is required if one is not ', 'identified in the blueprint.'
                )
                raise RuntimeError(msg)

        file_name = None
        if args.local:
            file_name = args.local[i]

        obs = _augment(
            obs,
            product_id,
            uri,
            blueprint,
            subject,
            args.dumpconfig,
            validate_wcs,
            plugin=None,
            local=file_name,
            client=client,
        )

    _write_observation(obs, args)


def _load_plugin(plugin_name):
    plgin = _load_module(plugin_name)
    if hasattr(plgin, 'update'):
        pass
    elif hasattr(plgin, 'ObservationUpdater'):
        # for backwards compatibility with caom2repo
        plgin = getattr(plgin, 'ObservationUpdater')()

    if not hasattr(plgin, 'update'):
        msg = (
            f'The plugin {plugin_name} is not correct.  It must provide one of:\n'
            '1 - a function named update, or\n'
            '2 - a class ObservationUpdater with a function named update.\n '
            'In either case, the update signature needs to be (Observation, **kwargs).'
        )
        raise ImportError(msg)
    return plgin


def _visit(plugin_name, parser, obs, visit_local, product_id=None, uri=None, subject=None, **kwargs):
    result = obs
    if plugin_name is not None and len(plugin_name) > 0:
        # TODO make a check that's necessary under both calling conditions here
        logging.debug(
            'Begin plugin execution {!r} update method on ' 'observation {!r}'.format(plugin_name, obs.observation_id)
        )
        plgin = _load_plugin(plugin_name)
        if isinstance(parser, FitsParser):
            kwargs['headers'] = parser.headers
        if visit_local is not None:
            kwargs['fqn'] = visit_local
        if product_id is not None:
            kwargs['product_id'] = product_id
        if uri is not None:
            kwargs['uri'] = uri
        if subject is not None:
            kwargs['subject'] = subject
        try:
            result = plgin.update(observation=obs, **kwargs)
            if result is not None:
                logging.debug(
                    'Finished executing plugin {!r} update '
                    'method on observation {!r}'.format(plugin_name, obs.observation_id)
                )
        except Exception as e:
            logging.error(e)
            tb = traceback.format_exc()
            logging.debug(tb)
            raise e
    return result


def _write_observation(obs, args):
    writer = ObservationWriter()
    if args.out_obs_xml:
        writer.write(obs, args.out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)


def gen_proc(args, blueprints, **kwargs):
    """The implementation that expects a product ID to be provided as part of the lineage parameter, and blueprints
    as input parameters, and a plugin parameter, that supports external programmatic blueprint modification."""
    _set_logging(args.verbose, args.debug, args.quiet)
    result = 0

    if args.in_obs_xml:
        obs = _gen_obs(blueprints, args.in_obs_xml)
    else:
        obs = _gen_obs(blueprints, None, args.observation[0], args.observation[1])

    validate_wcs = True
    if args.no_validate:
        validate_wcs = False
    connected = True
    client = None
    subject = None
    if args.not_connected:
        connected = False
    else:
        subject = net.Subject.from_cmd_line_args(args)
        if args.resource_id is None:
            # if the resource_id is Undefined, using CadcDataClient
            client = data_util.StorageClientWrapper(subject, using_storage_inventory=False)
        else:
            # if the resource_id is defined, assume that the caller intends to use the Storage Inventory system, as
            # it's the CADC storage client that depends on a resource_id
            client = data_util.StorageClientWrapper(subject, resource_id=args.resource_id)

    for ii, cardinality in enumerate(args.lineage):
        product_id, uri = _extract_ids(cardinality)
        blueprint = _lookup_blueprint(blueprints, uri)
        logging.debug('Begin augmentation for product_id {}, uri {}'.format(product_id, uri))

        file_name = None
        if args.local:
            file_name = args.local[ii]

        external_url = None
        if args.external_url:
            external_url = args.external_url[ii]

        use_blueprint_parser = False
        if args.use_blueprint_parser:
            use_blueprint_parser = uri in args.use_blueprint_parser

        obs = _augment(
            obs,
            product_id,
            uri,
            blueprint,
            subject,
            args.dumpconfig,
            validate_wcs,
            args.plugin,
            file_name,
            external_url,
            connected,
            use_blueprint_parser,
            client,
            **kwargs,
        )

        if obs is None:
            logging.warning('No observation. Stop processing.')
            break

    if obs is None:
        if args.in_obs_xml:
            log_id = args.lineage
        else:
            log_id = args.observation
        logging.warning(f'No Observation generated for {log_id}')
        result = -1
    else:
        _write_observation(obs, args)
    return result


def get_gen_proc_arg_parser():
    """
    Returns the arg parser with minimum arguments required to run caom2gen
    :return: args parser
    """
    parser = _get_common_arg_parser()
    parser.add_argument(
        '--external_url',
        nargs='+',
        help=(
            'service endpoint(s) that return(s) a string that can be made into FITS headers. Cardinality should '
            'be consistent with lineage.'
        ),
    )
    parser.add_argument(
        '--module',
        help=(
            'if the blueprint contains function calls, call importlib.import_module for the named module. Provide a '
            'fully qualified name. Parameter choices are the artifact URI (uri) or a list of astropy Header '
            'instances (header). This will allow the update of a single blueprint entry with a single call.'
        ),
    )
    parser.add_argument(
        '--plugin',
        help=(
            'if this parameter is specified, call importlib.import_module for the named module. Then execute the '
            'method "update", with the signature (Observation, **kwargs). This will allow for the update of '
            'multiple observation data members with one call.'
        ),
    )
    parser.add_argument(
        '--lineage',
        nargs='+',
        help=(
            'productID/artifactURI. List of plane/artifact identifiers that will be created for the identified '
            'observation.'
        ),
    )
    parser.add_argument(
        '--use_blueprint_parser',
        nargs='+',
        help=(
            'productID/artifactURI. List of lineage entries that will be processed with a BlueprintParser. '
            'Good for files with no metadata in the content.'
        ),
    )
    return parser


def augment(
    blueprints,
    no_validate=False,
    dump_config=False,
    plugin=None,
    out_obs_xml=None,
    in_obs_xml=None,
    collection=None,
    observation=None,
    product_id=None,
    uri=None,
    netrc=False,
    file_name=None,
    verbose=False,
    debug=False,
    quiet=False,
    **kwargs,
):
    _set_logging(verbose, debug, quiet)
    logging.debug('Begin augmentation for product_id {}, uri {}'.format(product_id, uri))

    # The 'visit_args' are a dictionary within the 'params' dictionary. They are set by the collection-specific
    # implementation, as they are dependent on that collection-specific implementation. The args to the visit
    # function are not set in fits2caom2.

    params = kwargs.get('params')
    kwargs = {}
    if params is not None:
        kwargs = params.get('visit_args')

    obs = _gen_obs(blueprints, in_obs_xml, collection, observation)
    subject = net.Subject(username=None, certificate=None, netrc=netrc)
    validate_wcs = True
    if no_validate is not None:
        validate_wcs = not no_validate

    for ii in blueprints:
        obs = _augment(
            obs, product_id, uri, blueprints[ii], subject, dump_config, validate_wcs, plugin, file_name, **kwargs
        )

    writer = ObservationWriter()
    if out_obs_xml:
        writer.write(obs, out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)
    logging.info('Done augment.')


augment.__doc__ = get_gen_proc_arg_parser().format_help()
