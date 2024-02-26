# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
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
#  Revision: 4
#
# ***********************************************************************
#

import logging

from caom2 import CalibrationLevel, DataProductType, ReleaseType
from caom2.caom_util import int_32


__all__ = [
    'Hdf5ObsBlueprint',
    'ObsBlueprint',
    '_to_float',
    '_to_int',
    '_to_int_32',
    '_to_str',
]


class classproperty:
    """
    Class property used for CAOM2_ELEMENTS in ObsBleprint
    """

    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class ObsBlueprint:
    """
    Class that represents the blueprint of a CAOM2 Observation that can be used to build an observation.

    The following CAOM2 elements can be specified in the blueprint:
    _CAOM2_ELEMENTS

    The blueprint designates the source of each of these attributes as either FITS keywords with possible default
    values or sets the actual values. The blueprint can be checked by simply displaying it.

    For example:

    # display the default blueprint when WCS axes are not specified
    print(ObsBlueprint())

    # display the default blueprint when WCS axes are specified
    print(ObsBlueprint(position_axis=(1, 2), energy_axis=3, polarization_axis=4, time_axis=5))

    # create a blueprint and customize it
    ob = ObsBlueprint(position_axis=(1, 2), energy_axis=3, polarization_axis=4, time_axis=5))
    ob.set('Observation.algorithm.name', 'exposure')
    ob.add_attribute('Chunk.energy.axis.axis.ctype', ['MYCTYPE'], extension=1)
    ob.add_attribute('Chunk.energy.axis.axis.ctype', 'MYCTYPE2', extension=1)
    ob.set('Chunk.energy.velang', 33, extension=1)
    ob.set_default('Chunk.position.coordsys', 'RA-DEC', extension=1)

    ob.set('Chunk.energy.velang', 44, extension=2)
    print(ob)

    """

    _CAOM2_ELEMENTS = [
        'CompositeObservation.members',
        'DerivedObservation.members',
        'Observation.observationID',
        'Observation.type',
        'Observation.intent',
        'Observation.sequenceNumber',
        'Observation.metaRelease',
        'Observation.metaReadGroups',
        'Observation.metaProducer',
        'Observation.requirements.flag',
        'Observation.algorithm.name',
        'Observation.instrument.name',
        'Observation.instrument.keywords',
        'Observation.proposal.id',
        'Observation.proposal.pi',
        'Observation.proposal.project',
        'Observation.proposal.title',
        'Observation.proposal.keywords',
        'Observation.target.name',
        'Observation.target.type',
        'Observation.target.standard',
        'Observation.target.redshift',
        'Observation.target.keywords',
        'Observation.target.moving',
        'Observation.target.targetID',
        'Observation.target_position.point.cval1',
        'Observation.target_position.point.cval2',
        'Observation.target_position.coordsys',
        'Observation.target_position.equinox',
        'Observation.telescope.name',
        'Observation.telescope.geoLocationX',
        'Observation.telescope.geoLocationY',
        'Observation.telescope.geoLocationZ',
        'Observation.telescope.keywords',
        'Observation.environment.seeing',
        'Observation.environment.humidity',
        'Observation.environment.elevation',
        'Observation.environment.tau',
        'Observation.environment.wavelengthTau',
        'Observation.environment.ambientTemp',
        'Observation.environment.photometric',
        'Plane.productID',
        'Plane.metaRelease',
        'Plane.dataRelease',
        'Plane.dataProductType',
        'Plane.calibrationLevel',
        'Plane.dataQuality',
        'Plane.metaReadGroups',
        'Plane.dataReadGroups',
        'Plane.metaProducer',
        'Plane.provenance.name',
        'Plane.provenance.version',
        'Plane.provenance.project',
        'Plane.provenance.producer',
        'Plane.provenance.runID',
        'Plane.provenance.reference',
        'Plane.provenance.lastExecuted',
        'Plane.provenance.keywords',
        'Plane.provenance.inputs',
        'Plane.metrics.sourceNumberDensity',
        'Plane.metrics.background',
        'Plane.metrics.backgroundStddev',
        'Plane.metrics.fluxDensityLimit',
        'Plane.metrics.magLimit',
        'Plane.metrics.sampleSNR',
        'Plane.observable.ucd',
        'Artifact.productType',
        'Artifact.releaseType',
        'Artifact.contentChecksum',
        'Artifact.contentLength',
        'Artifact.contentType',
        'Artifact.contentRelease',
        'Artifact.contentReadGroups',
        'Artifact.uri',
        'Artifact.metaProducer',
        'Part.name',
        'Part.productType',
        'Part.metaProducer',
        'Chunk',
        'Chunk.naxis',
        'Chunk.observableAxis',
        'Chunk.positionAxis1',
        'Chunk.positionAxis2',
        'Chunk.energyAxis',
        'Chunk.timeAxis',
        'Chunk.polarizationAxis',
        'Chunk.metaProducer',
        'Chunk.observable.dependent.bin',
        'Chunk.observable.dependent.axis.ctype',
        'Chunk.observable.dependent.axis.cunit',
        'Chunk.observable.independent.bin',
        'Chunk.observable.independent.axis.ctype',
        'Chunk.observable.independent.axis.cunit',
        'Chunk.position.coordsys',
        'Chunk.position.equinox',
        'Chunk.position.resolution',
        'Chunk.position.axis.axis1.ctype',
        'Chunk.position.axis.axis1.cunit',
        'Chunk.position.axis.axis2.ctype',
        'Chunk.position.axis.axis2.cunit',
        'Chunk.position.axis.error1.syser',
        'Chunk.position.axis.error1.rnder',
        'Chunk.position.axis.error2.syser',
        'Chunk.position.axis.error2.rnder',
        'Chunk.position.axis.function.cd11',
        'Chunk.position.axis.function.cd12',
        'Chunk.position.axis.function.cd21',
        'Chunk.position.axis.function.cd22',
        'Chunk.position.axis.function.dimension.naxis1',
        'Chunk.position.axis.function.dimension.naxis2',
        'Chunk.position.axis.function.refCoord.coord1.pix',
        'Chunk.position.axis.function.refCoord.coord1.val',
        'Chunk.position.axis.function.refCoord.coord2.pix',
        'Chunk.position.axis.function.refCoord.coord2.val',
        'Chunk.position.axis.range.start.coord1.pix',
        'Chunk.position.axis.range.start.coord1.val',
        'Chunk.position.axis.range.start.coord2.pix',
        'Chunk.position.axis.range.start.coord2.val',
        'Chunk.position.axis.range.end.coord1.pix',
        'Chunk.position.axis.range.end.coord1.val',
        'Chunk.position.axis.range.end.coord2.pix',
        'Chunk.position.axis.range.end.coord2.val',
        'Chunk.energy.specsys',
        'Chunk.energy.ssysobs',
        'Chunk.energy.restfrq',
        'Chunk.energy.restwav',
        'Chunk.energy.velosys',
        'Chunk.energy.zsource',
        'Chunk.energy.ssyssrc',
        'Chunk.energy.velang',
        'Chunk.energy.bandpassName',
        'Chunk.energy.resolvingPower',
        'Chunk.energy.transition',
        'Chunk.energy.transition.species',
        'Chunk.energy.transition.transition',
        'Chunk.energy.axis.axis.ctype',
        'Chunk.energy.axis.axis.cunit',
        'Chunk.energy.axis.bounds.samples',
        'Chunk.energy.axis.error.syser',
        'Chunk.energy.axis.error.rnder',
        'Chunk.energy.axis.function.naxis',
        'Chunk.energy.axis.function.delta',
        'Chunk.energy.axis.function.refCoord.pix',
        'Chunk.energy.axis.function.refCoord.val',
        'Chunk.energy.axis.range.start.pix',
        'Chunk.energy.axis.range.start.val',
        'Chunk.energy.axis.range.end.pix',
        'Chunk.energy.axis.range.end.val',
        'Chunk.polarization.axis.axis.ctype',
        'Chunk.polarization.axis.axis.cunit',
        'Chunk.polarization.axis.bounds.samples',
        'Chunk.polarization.axis.error.syser',
        'Chunk.polarization.axis.error.rnder',
        'Chunk.polarization.axis.function.naxis',
        'Chunk.polarization.axis.function.delta',
        'Chunk.polarization.axis.function.refCoord.pix',
        'Chunk.polarization.axis.function.refCoord.val',
        'Chunk.polarization.axis.range.start.pix',
        'Chunk.polarization.axis.range.start.val',
        'Chunk.polarization.axis.range.end.pix',
        'Chunk.polarization.axis.range.end.val',
        'Chunk.time.exposure',
        'Chunk.time.resolution',
        'Chunk.time.timesys',
        'Chunk.time.trefpos',
        'Chunk.time.mjdref',
        'Chunk.time.axis.axis.ctype',
        'Chunk.time.axis.axis.cunit',
        'Chunk.time.axis.bounds.samples',
        'Chunk.time.axis.error.syser',
        'Chunk.time.axis.error.rnder',
        'Chunk.time.axis.function.naxis',
        'Chunk.time.axis.function.delta',
        'Chunk.time.axis.function.refCoord.pix',
        'Chunk.time.axis.function.refCoord.val',
        'Chunk.time.axis.range.start.pix',
        'Chunk.time.axis.range.start.val',
        'Chunk.time.axis.range.end.pix',
        'Chunk.time.axis.range.end.val',
        'Chunk.observable.axis.axis.ctype',
        'Chunk.observable.axis.axis.cunit',
        'Chunk.observable.axis.function.refCoord.pix',
        'Chunk.custom.axis.axis.ctype',
        'Chunk.custom.axis.axis.cunit',
        'Chunk.custom.axis.bounds.samples',
        'Chunk.custom.axis.error.syser',
        'Chunk.custom.axis.error.rnder',
        'Chunk.custom.axis.function.naxis',
        'Chunk.custom.axis.function.delta',
        'Chunk.custom.axis.function.refCoord.pix',
        'Chunk.custom.axis.function.refCoord.val',
        'Chunk.custom.axis.range.start.pix',
        'Chunk.custom.axis.range.start.val',
        'Chunk.custom.axis.range.end.pix',
        'Chunk.custom.axis.range.end.val',
    ]

    # replace _CAOM2_ELEMENTS in __doc__ with the real elements
    __doc__ = __doc__.replace('_CAOM2_ELEMENTS', '\n'.join(['\t\t{}'.format(elem) for elem in _CAOM2_ELEMENTS]))

    def __init__(
        self,
        position_axes=None,
        energy_axis=None,
        polarization_axis=None,
        time_axis=None,
        obs_axis=None,
        custom_axis=None,
        module=None,
        update=True,
        instantiated_class=None,
    ):
        """
        Ctor
        :param position_axes: tuple of form (int, int) indicating the indexes of position axis
        :param energy_axis: index of energy axis (int)
        :param polarization_axis: index of polarization axis (int)
        :param time_axis: index of time axis (int)
        :param obs_axis: index of observable axis (int)
        :param custom_axis: index of custom axis (int)
        :param module: user-provided code, will be loaded with importlib.import_module if a value is provided.
        """

        if position_axes and isinstance(position_axes, tuple) and (len(position_axes) != 2):
            raise ValueError('Invalid position axis: {}. Must be tuple with 2 elements'.format(str(position_axes)))

        self.logger = logging.getLogger(__name__)

        # this is the default blueprint
        self._plan = {}
        tmp = {
            'Observation.metaRelease': (
                ['DATE', 'DATE-OBS', 'UTCOBS', 'UTCDATE', 'UTC-DATE', 'MJDOBS', 'MJD_OBS'],
                None,
            ),
            'Observation.instrument.name': (['INSTRUME'], None),
            'Observation.type': (['OBSTYPE'], None),
            'Observation.environment.ambientTemp': (['TEMPERAT'], None),
            # set the default for SimpleObservation construction
            'Observation.algorithm.name': (['PROCNAME'], 'exposure'),
            'Observation.instrument.keywords': (['INSTMODE'], None),
            'Observation.proposal.id': (['RUNID'], None),
            'Observation.target.name': (['OBJECT'], None),
            'Observation.telescope.name': (['TELESCOP'], None),
            'Observation.telescope.geoLocationX': (['OBSGEO-X'], None),
            'Observation.telescope.geoLocationY': (['OBSGEO-Y'], None),
            'Observation.telescope.geoLocationZ': (['OBSGEO-Z'], None),
            'Observation.observationID': (['OBSID'], None),
            'Plane.calibrationLevel': ([], CalibrationLevel.RAW_STANDARD),
            'Plane.dataProductType': ([], DataProductType.IMAGE),
            'Plane.metaRelease': (['RELEASE', 'REL_DATE'], None),
            'Plane.dataRelease': (['RELEASE', 'REL_DATE'], None),
            'Plane.productID': (['RUNID'], None),
            'Plane.provenance.name': (['XPRVNAME'], None),
            'Plane.provenance.project': (['ADC_ARCH'], None),
            'Plane.provenance.producer': (['ORIGIN'], None),
            'Plane.provenance.reference': (['XREFER'], None),
            'Plane.provenance.lastExecuted': (['DATE-FTS'], None),
            'Artifact.releaseType': ([], ReleaseType.DATA),
            'Chunk': 'include',
        }
        # using the tmp to make sure that the keywords are valid
        for key in tmp:
            self.set(key, tmp[key])

        self._extensions = {}

        # contains the standard WCS keywords in the FITS file expected by the astropy.WCS package.
        self._wcs_std = {'Chunk.naxis': 'ZNAXIS,NAXIS'}
        self._pos_axes_configed = False
        self._energy_axis_configed = False
        self._time_axis_configed = False
        self._polarization_axis_configed = False
        self._obs_axis_configed = False
        self._custom_axis_configed = False
        if position_axes:
            self.configure_position_axes(position_axes)

        if energy_axis:
            self.configure_energy_axis(energy_axis)

        if polarization_axis:
            self.configure_polarization_axis(polarization_axis)

        if time_axis:
            self.configure_time_axis(time_axis)

        if obs_axis:
            self.configure_observable_axis(obs_axis)

        if custom_axis:
            self.configure_custom_axis(custom_axis)

        if module:
            self._module = module
        else:
            self._module = None
        self._module_instance = instantiated_class
        # if True, existing values are used instead of defaults
        self._update = update
        # a data structure to carry around twelve bits of data at a time:
        # the first item in the set is the ctype index, and the second is whether or not the index means anything,
        # resulting in a call to the blueprint configure_* methods if it's True.
        self._axis_info = {
            'custom': (0, False),
            'dec': (0, False),
            'energy': (0, False),
            'obs': (0, False),
            'polarization': (0, False),
            'ra': (0, False),
            'time': (0, False),
        }

    def configure_custom_axis(self, axis, override=True):
        """
        Set the expected FITS custom keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the custom axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._custom_axis_configed:
            self.logger.debug('Attempt to configure already-configured custom axis.')
            return

        if override:
            self.set('Chunk.custom.axis.axis.ctype', ([f'CTYPE{axis}'], None))
            self.set('Chunk.custom.axis.axis.cunit', ([f'CUNIT{axis}'], None))
            self.set('Chunk.custom.axis.function.naxis', ([f'NAXIS{axis}'], None))
            self.set('Chunk.custom.axis.function.delta', ([f'CDELT{axis}'], None))
            self.set('Chunk.custom.axis.function.refCoord.pix', ([f'CRPIX{axis}'], None))
            self.set('Chunk.custom.axis.function.refCoord.val', ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.custom.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.custom.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.custom.axis.function.naxis'] = f'NAXIS{axis}'
        self._wcs_std['Chunk.custom.axis.function.delta'] = f'CDELT{axis}'
        self._wcs_std['Chunk.custom.axis.function.refCoord.pix'] = f'CRPIX{axis}'
        self._wcs_std['Chunk.custom.axis.function.refCoord.val'] = f'CRVAL{axis}'

        self._custom_axis_configed = True

    def configure_position_axes(self, axes, override=True):
        """
        Set the expected FITS spatial keywords by indices in the blueprint and the wcs_std lookup.

        :param axes: The index expected for the position axes.
        :return:
        """
        if self._pos_axes_configed:
            self.logger.debug('Attempt to configure already-configured position axes.')
            return

        if override:
            self.set('Chunk.position.coordsys', (['RADESYS'], None))
            self.set('Chunk.position.equinox', (['EQUINOX', 'EPOCH'], None))
            self.set('Chunk.position.axis.axis1.ctype', ([f'CTYPE{axes[0]}'], None))
            self.set('Chunk.position.axis.axis1.cunit', ([f'CUNIT{axes[0]}'], None))
            self.set('Chunk.position.axis.axis2.ctype', ([f'CTYPE{axes[1]}'], None))
            self.set('Chunk.position.axis.axis2.cunit', ([f'CUNIT{axes[1]}'], None))
            self.set('Chunk.position.axis.error1.syser', ([f'CSYER{axes[0]}'], None))
            self.set('Chunk.position.axis.error1.rnder', ([f'CRDER{axes[0]}'], None))
            self.set('Chunk.position.axis.error2.syser', ([f'CSYER{axes[1]}'], None))
            self.set('Chunk.position.axis.error2.rnder', ([f'CRDER{axes[1]}'], None))
            self.set('Chunk.position.axis.function.cd11', ([f'CD{axes[0]}_{axes[0]}'], None))
            self.set('Chunk.position.axis.function.cd12', ([f'CD{axes[0]}_{axes[1]}'], None))
            self.set('Chunk.position.axis.function.cd21', ([f'CD{axes[1]}_{axes[0]}'], None))
            self.set('Chunk.position.axis.function.cd22', ([f'CD{axes[1]}_{axes[1]}'], None))
            self.set('Chunk.position.axis.function.dimension.naxis1', ([f'ZNAXIS{axes[0]}', f'NAXIS{axes[0]}'], None))
            self.set('Chunk.position.axis.function.dimension.naxis2', ([f'ZNAXIS{axes[1]}', f'NAXIS{axes[1]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.pix', ([f'CRPIX{axes[0]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.val', ([f'CRVAL{axes[0]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.pix', ([f'CRPIX{axes[1]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.val', ([f'CRVAL{axes[1]}'], None))

        self._wcs_std['Chunk.position.coordsys'] = 'RADESYS'
        self._wcs_std['Chunk.position.equinox'] = 'EQUINOX'

        self._wcs_std['Chunk.position.axis.axis1.ctype'] = f'CTYPE{axes[0]}'
        self._wcs_std['Chunk.position.axis.axis1.cunit'] = f'CUNIT{axes[0]}'
        self._wcs_std['Chunk.position.axis.axis2.ctype'] = f'CTYPE{axes[1]}'
        self._wcs_std['Chunk.position.axis.axis2.cunit'] = f'CUNIT{axes[1]}'
        self._wcs_std['Chunk.position.axis.error1.syser'] = f'CSYER{axes[0]}'
        self._wcs_std['Chunk.position.axis.error1.rnder'] = f'CRDER{axes[0]}'
        self._wcs_std['Chunk.position.axis.error2.syser'] = f'CSYER{axes[1]}'
        self._wcs_std['Chunk.position.axis.error2.rnder'] = f'CRDER{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.cd11'] = f'CD{axes[0]}_{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.cd12'] = f'CD{axes[0]}_{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.cd21'] = f'CD{axes[1]}_{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.cd22'] = f'CD{axes[1]}_{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = f'NAXIS{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = f'NAXIS{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] = f'CRPIX{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] = f'CRVAL{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] = f'CRPIX{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] = f'CRVAL{axes[1]}'

        self._pos_axes_configed = True

    def configure_energy_axis(self, axis, override=True):
        """
        Set the expected FITS energy keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the energy axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._energy_axis_configed:
            self.logger.debug('Attempt to configure already-configured energy axis.')
            return

        if override:
            self.set('Chunk.energy.specsys', (['SPECSYS'], None))
            self.set('Chunk.energy.ssysobs', (['SSYSOBS'], None))
            self.set('Chunk.energy.restfrq', (['RESTFRQ'], None))
            self.set('Chunk.energy.restwav', (['RESTWAV'], None))
            self.set('Chunk.energy.velosys', (['VELOSYS'], None))
            self.set('Chunk.energy.zsource', (['ZSOURCE'], None))
            self.set('Chunk.energy.ssyssrc', (['SSYSSRC'], None))
            self.set('Chunk.energy.velang', (['VELANG'], None))

            self.set('Chunk.energy.bandpassName', ([], None))
            self.set('Chunk.energy.resolvingPower', ([], None))

            self.set('Chunk.energy.axis.axis.ctype', ([f'CTYPE{axis}'], None))
            self.set('Chunk.energy.axis.axis.cunit', ([f'CUNIT{axis}'], None))
            self.set('Chunk.energy.axis.error.syser', ([f'CSYER{axis}'], None))
            self.set('Chunk.energy.axis.error.rnder', ([f'CRDER{axis}'], None))
            self.set('Chunk.energy.axis.function.naxis', ([f'NAXIS{axis}'], None))
            self.set('Chunk.energy.axis.function.delta', ([f'CDELT{axis}'], None))
            self.set('Chunk.energy.axis.function.refCoord.pix', ([f'CRPIX{axis}'], None))
            self.set('Chunk.energy.axis.function.refCoord.val', ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.energy.specsys'] = 'SPECSYS'
        self._wcs_std['Chunk.energy.ssysobs'] = 'SSYSOBS'
        self._wcs_std['Chunk.energy.restfrq'] = 'RESTFRQ'
        self._wcs_std['Chunk.energy.restwav'] = 'RESTWAV'
        self._wcs_std['Chunk.energy.velosys'] = 'VELOSYS'
        self._wcs_std['Chunk.energy.zsource'] = 'ZSOURCE'
        self._wcs_std['Chunk.energy.ssyssrc'] = 'SSYSSRC'
        self._wcs_std['Chunk.energy.velang'] = 'VELANG'

        self._wcs_std['Chunk.energy.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.energy.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.energy.axis.error.syser'] = f'CSYER{axis}'
        self._wcs_std['Chunk.energy.axis.error.rnder'] = f'CRDER{axis}'
        self._wcs_std['Chunk.energy.axis.function.naxis'] = f'NAXIS{axis}'
        self._wcs_std['Chunk.energy.axis.function.delta'] = f'CDELT{axis}'
        self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = f'CRPIX{axis}'
        self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = f'CRVAL{axis}'

        self._energy_axis_configed = True

    def configure_polarization_axis(self, axis, override=True):
        """
        Set the expected FITS polarization keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the polarization axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._polarization_axis_configed:
            self.logger.debug('Attempt to configure already-configured polarization axis.')
            return

        if override:
            self.set('Chunk.polarization.axis.axis.ctype', ([f'CTYPE{axis}'], None))
            self.set('Chunk.polarization.axis.axis.cunit', ([f'CUNIT{axis}'], None))
            self.set('Chunk.polarization.axis.function.naxis', ([f'NAXIS{axis}'], None))
            self.set('Chunk.polarization.axis.function.delta', ([f'CDELT{axis}'], None))
            self.set('Chunk.polarization.axis.function.refCoord.pix', ([f'CRPIX{axis}'], None))
            self.set('Chunk.polarization.axis.function.refCoord.val', ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.polarization.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.polarization.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.polarization.axis.function.naxis'] = f'NAXIS{axis}'
        self._wcs_std['Chunk.polarization.axis.function.delta'] = f'CDELT{axis}'
        self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = f'CRPIX{axis}'
        self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = f'CRVAL{axis}'

        self._polarization_axis_configed = True

    def configure_observable_axis(self, axis, override=True):
        """
        Set the expected FITS observable keywords by index in the blueprint and the wcs_std lookup.
        Note: observable axis is not a standard WCS and it's not used by astropy.wcs so, arguably, it can be removed.
        It is here for now for consistency purposes.
        :param axis: The index expected for the observable axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._obs_axis_configed:
            self.logger.debug('Attempt to configure already-configured observable axis.')
            return

        if override:
            self.set('Chunk.observable.axis.axis.ctype', ([f'CTYPE{axis}'], None))
            self.set('Chunk.observable.axis.axis.cunit', ([f'CUNIT{axis}'], None))
            self.set('Chunk.observable.axis.function.refCoord.pix', ([f'CRPIX{axis}'], None))

        self._wcs_std['Chunk.observable.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.observable.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.observable.axis.function.refCoord.pix'] = f'CRPIX{axis}'

        self._obs_axis_configed = True

    def configure_time_axis(self, axis, override=True):
        """
        Set the expected FITS time keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the time axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._time_axis_configed:
            self.logger.debug('Attempt to configure already-configured time axis.')
            return

        if override:
            self.set('Chunk.time.exposure', (['EXPTIME', 'INTTIME'], None))
            self.set('Chunk.time.timesys', (['TIMESYS'], None))
            self.set('Chunk.time.trefpos', (['TREFPOS'], None))
            self.set('Chunk.time.mjdref', (['MJDREF'], None))
            self.set('Chunk.time.resolution', (['TIMEDEL'], None))
            self.set('Chunk.time.axis.axis.ctype', ([f'CTYPE{axis}'], None))
            self.set('Chunk.time.axis.axis.cunit', ([f'CUNIT{axis}'], None))
            self.set('Chunk.time.axis.error.syser', ([f'CSYER{axis}'], None))
            self.set('Chunk.time.axis.error.rnder', ([f'CRDER{axis}'], None))
            self.set('Chunk.time.axis.function.naxis', ([f'NAXIS{axis}'], None))
            self.set('Chunk.time.axis.function.delta', ([f'CDELT{axis}'], None))
            self.set('Chunk.time.axis.function.refCoord.pix', ([f'CRPIX{axis}'], None))
            self.set('Chunk.time.axis.function.refCoord.val', ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.time.exposure'] = 'EXPTIME'
        self._wcs_std['Chunk.time.resolution'] = 'TIMEDEL'
        self._wcs_std['Chunk.time.timesys'] = 'TIMESYS'
        self._wcs_std['Chunk.time.trefpos'] = 'TREFPOS'
        self._wcs_std['Chunk.time.mjdref'] = 'MJDREF'

        self._wcs_std['Chunk.time.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.time.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.time.axis.error.syser'] = f'CSYER{axis}'
        self._wcs_std['Chunk.time.axis.error.rnder'] = f'CRDER{axis}'
        self._wcs_std['Chunk.time.axis.function.naxis'] = f'NAXIS{axis}'
        self._wcs_std['Chunk.time.axis.function.delta'] = f'CDELT{axis}'
        self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = f'CRPIX{axis}'
        self._wcs_std['Chunk.time.axis.function.refCoord.val'] = f'CRVAL{axis}'

        self._time_axis_configed = True

    def _guess_axis_info(self):
        """Look for info regarding axis types in the blueprint wcs_std. Configure the blueprint according to the
        guesses.
        """
        for ii in self._plan:
            if isinstance(self._plan[ii], tuple):
                for value in self._plan[ii][0]:
                    if (value.startswith('CTYPE')) and value[-1].isdigit():
                        value = value.split('-')[0]
                        self._guess_axis_info_from_ctypes(ii, int(value[-1]))
            else:
                value = self._plan[ii]
                if value is None:
                    continue
                if (value.startswith('CTYPE')) and value[-1].isdigit():
                    value = value.split('-')[0]
                    self._guess_axis_info_from_ctypes(ii, int(value[-1]))

        self._guess_axis_info_from_plan()

    def _guess_axis_info_from_plan(self):
        for ii in self._plan:
            if ii.startswith('Chunk.position') and ii.endswith('axis1.ctype') and not self._axis_info['ra'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'ra')
                self._axis_info['ra'] = (configured_index, True)
            elif ii.startswith('Chunk.position') and ii.endswith('axis2.ctype') and not self._axis_info['dec'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'dec')
                self._axis_info['dec'] = (configured_index, True)
            elif ii.startswith('Chunk.energy') and not self._axis_info['energy'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'energy')
                self._axis_info['energy'] = (configured_index, True)
            elif ii.startswith('Chunk.time') and not self._axis_info['time'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'time')
                self._axis_info['time'] = (configured_index, True)
            elif ii.startswith('Chunk.polarization') and not self._axis_info['polarization'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'polarization')
                self._axis_info['polarization'] = (configured_index, True)
            elif ii.startswith('Chunk.observable') and not self._axis_info['obs'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'obs')
                self._axis_info['obs'] = (configured_index, True)
            elif ii.startswith('Chunk.custom') and not self._axis_info['custom'][1]:
                configured_index = self._get_configured_index(self._axis_info, 'custom')
                self._axis_info['custom'] = (configured_index, True)

        if self._axis_info['ra'][1] and self._axis_info['dec'][1]:
            self.configure_position_axes((self._axis_info['ra'][0], self._axis_info['dec'][0]), False)
        elif self._axis_info['ra'][1] or self._axis_info['dec'][1]:
            raise ValueError(
                'Only one positional axis found '
                '(ra/dec): {}/{}'.format(self._axis_info['ra'][0], self._axis_info['dec'][0])
            )
        else:
            # assume that positional axis are 1 and 2 by default
            if (
                self._axis_info['time'][0] in [1, 2]
                or self._axis_info['energy'][0] in [1, 2]
                or self._axis_info['polarization'][0] in [1, 2]
                or self._axis_info['obs'][0] in [1, 2]
                or self._axis_info['custom'][0] in [1, 2]
            ):
                raise ValueError('Cannot determine the positional axis')
            else:
                self.configure_position_axes((1, 2), False)

        if self._axis_info['time'][1]:
            self.configure_time_axis(self._axis_info['time'][0], False)
        if self._axis_info['energy'][1]:
            self.configure_energy_axis(self._axis_info['energy'][0], False)
        if self._axis_info['polarization'][1]:
            self.configure_polarization_axis(self._axis_info['polarization'][0], False)
        if self._axis_info['obs'][1]:
            self.configure_observable_axis(self._axis_info['obs'][0], False)
        if self._axis_info['custom'][1]:
            self.configure_custom_axis(self._axis_info['custom'][0], False)

    def _guess_axis_info_from_ctypes(self, lookup, counter):
        """
        Check for the presence of blueprint keys in the plan, and whether or
        not they indicate an index in their configuration.

        :param lookup: Blueprint plan key.
        :param counter: Value to set the index to for an axis.
        :param axis_info: local data structure to pass around what is
            configured, and what is it's value.
        """
        if lookup.startswith('Chunk.energy'):
            self._axis_info['energy'] = (counter, True)
        elif lookup.startswith('Chunk.polarization'):
            self._axis_info['polarization'] = (counter, True)
        elif lookup.startswith('Chunk.time'):
            self._axis_info['time'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith('axis1.ctype'):
            self._axis_info['ra'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith('axis2.ctype'):
            self._axis_info['dec'] = (counter, True)
        elif lookup.startswith('Chunk.observable'):
            self._axis_info['obs'] = (counter, True)
        elif lookup.startswith('Chunk.custom'):
            self._axis_info['custom'] = (counter, True)
        else:
            raise ValueError(f'Unrecognized axis type: {lookup}')

    def _get_configured_index(self, axis_info, lookup):
        """Find the next available index value among those that are not set.

        :param axis_info: local data structure to pass around what is
            configured, and what is it's value."""
        DEFAULT_INDICES = {'ra': 1, 'dec': 2, 'energy': 3, 'time': 4, 'polarization': 5, 'obs': 6, 'custom': 7}

        # the logic - if the default index is already used, assign the lowest index that is unused, otherwise use the
        # default index

        max_index = 0
        min_index = 7
        default_index = DEFAULT_INDICES[lookup]
        default_used = False
        for axis in axis_info:
            # do two unrelated things in this for loop
            # 1. determine where to start counting
            if axis_info[axis][1]:
                max_index = max(max_index, axis_info[axis][0])
                min_index = min(min_index, axis_info[axis][0])
            # 2. determine if the default is used
            if axis_info[axis][1] and default_index == axis_info[axis][0]:
                default_used = True

        configured_index = 0
        if default_used:
            if min_index == 1:
                configured_index = max_index + 1
            else:
                configured_index = min(1, min_index)
        else:
            configured_index = default_index
        return configured_index

    def load_from_file(self, file_name):
        """
        Load a blueprint from a file. The expected input format is the same as is output by _serialize. This means
        there's lots of stripping of extra spaces, equals signs, and the word default. Also manage square brackets
        as list construction.

        Accept comments that start with '#'.

        :param file_name: The fully-qualified pathname for the blueprint file on disk.
        """
        ext = 0
        with open(file_name) as file:
            for line in file:
                line = line.split('#')[0]
                if '=' in line:
                    key, value = line.split('=', 1)
                    if 'default' in value:
                        temp = value.split(', default')
                        default = temp[1].replace('=', '').strip()
                        temp_list = [
                            ii.replace('[', '').replace(']', '').replace('\'', '').strip() for ii in temp[0].split(',')
                        ]
                        if 'None' in default:
                            default = None
                        cleaned_up_value = (temp_list, default)
                    else:
                        if value.strip() and value.strip()[0] == '(':
                            cleaned_up_value = tuple(ii.strip() for ii in value.strip().replace(
                                    '(', ''
                                ).replace(')', '').replace('\'', '').split(','))
                        elif '[' in value:
                            temp_list = value.replace('[', '').replace(']', '').replace('\'', '').split(',')
                            temp_list_2 = []
                            for ii in temp_list:
                                temp_list_2.append(ii.strip().strip('\n'))
                            cleaned_up_value = (temp_list_2, None)
                        else:
                            cleaned_up_value = value.strip('\n').strip()
                            if cleaned_up_value == 'None':
                                cleaned_up_value = None
                    self.set(key.strip(), cleaned_up_value, ext)
                elif 'extension' in line:
                    # pattern is 'extension #:'
                    new_ext = _to_int(line.strip('extension').strip('\n').strip(':'))
                    if isinstance(new_ext, int):
                        self.logger.info(f'Add extension {new_ext} to blueprint.')
                        ext = new_ext
        self._guess_axis_info()

    @classproperty
    def CAOM2_ELEMENTS(cls):
        """
        List of valid names of CAOM2 elements.
        :return:
        """
        return list(ObsBlueprint._CAOM2_ELEMENTS)  # return a copy

    @classmethod
    def check_caom2_element(cls, caom2_element):
        """
        Checks that an element is a valid caom2_element in the blueprint. It checks that it's part of the
        ObsBlueprint._CAOM2_ELEMENTS
        :param caom2_element: name CAOM2 element to check
        :raises KeyError
        """
        if caom2_element not in cls._CAOM2_ELEMENTS:
            raise KeyError('{} not a valid CAOM2 element name (mispelling?).'.format(caom2_element))

    @staticmethod
    def check_chunk(caom2_element):
        """
        Checks that an element is a valid Chunk-type caom2_element
        :param caom2_element: name CAOM2 element to check
        :raises ValueError
        """
        if not caom2_element.startswith('Chunk'):
            raise ValueError("Extension number refers to Chunk elements only")

    @staticmethod
    def check_extension(extension):
        if extension is not None and extension < 0:
            raise ValueError(f'Extension count failure. {extension} should be >= 0')

    def __str__(self):
        plan = self._serialize(self._plan)

        extensions = ''
        if self._extensions:
            for key in sorted(self._extensions):
                extensions = extensions + f'\nextension {key}:\n' + self._serialize(self._extensions[key])
        return plan + extensions

    def _serialize(self, src):
        return '\n'.join(
            [
                '{} = {}'.format(
                    key,
                    '{}, default = {}'.format(src[key][0], src[key][1]) if isinstance(src[key], tuple) else src[key],
                )
                for key in ObsBlueprint._CAOM2_ELEMENTS
                if key in src
            ]
        )

    def set(self, caom2_element, value, extension=0):
        """
        Sets the value associated with an element in the CAOM2 model. Value cannot be a tuple.
        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param value: new value of the CAOM2 element
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = value
        else:
            self._plan[caom2_element] = value

    def add_attribute(self, caom2_element, attribute, extension=0):
        """
        Adds an attribute in the list of other attributes associated with an caom2 element.
        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param attribute: name of attribute the element is mapped to
        :param extension: extension number (used only for Chunk elements)
        :raises AttributeError if the caom2 element has already an associated value or KeyError if the caom2 element
        does not exists.
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise AttributeError(f'No extension {extension} in the blueprint')
            else:
                if caom2_element in self._extensions[extension]:
                    if isinstance(self._extensions[extension][caom2_element], tuple):
                        if attribute not in self._extensions[extension][caom2_element][0]:
                            self._extensions[extension][caom2_element][0].insert(0, attribute)
                    else:
                        raise AttributeError(
                            (f'No attributes in extension {extension} ' f'associated with keyword {caom2_element}')
                        )
                else:
                    self._extensions[extension][caom2_element] = ([attribute], None)
        else:
            if caom2_element in self._plan:
                if isinstance(self._plan[caom2_element], tuple):
                    if attribute not in self._plan[caom2_element][0]:
                        self._plan[caom2_element][0].insert(0, attribute)
                else:
                    raise AttributeError(f'No attributes associated with ' f'keyword {caom2_element}')
            else:
                self._plan[caom2_element] = ([attribute], None)

    def add_table_attribute(self, caom2_element, ttype_attribute, extension=0, index=0):
        """
        Adds a FITS BINTABLE TTYPE* lookup, to a list of other FITS attributes associated with an caom2 element.
        This does not co-exist with non-table attributes.

        There is no support for default values for table attributes.

        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param ttype_attribute: name of TTYPE attribute element is mapped to
        :param extension: extension number (used only for Chunk elements)
        :param index: which row values to return. If index is None, all row values will be returned as a
            comma-separated list.
        :raises AttributeError if the caom2 element has already an associated value or KeyError if the caom2
            element does not exists.
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            if extension in self._extensions:
                if caom2_element in self._extensions[extension]:
                    if ObsBlueprint.is_table(self._extensions[extension][caom2_element]):
                        if ttype_attribute not in self._extensions[extension][caom2_element][1]:
                            self._extensions[extension][caom2_element][1].insert(0, ttype_attribute)
                    else:
                        raise AttributeError(
                            ('No TTYPE attributes in extension {} associated ' 'with keyword {}').format(
                                extension, caom2_element
                            )
                        )
                else:
                    self._extensions[extension][caom2_element] = ('BINTABLE', [ttype_attribute], index)
            else:
                self._extensions[extension] = {}
                self._extensions[extension][caom2_element] = ('BINTABLE', [ttype_attribute], index)
        else:
            if caom2_element in self._plan:
                if ObsBlueprint.is_table(self._plan[caom2_element]):
                    if ttype_attribute not in self._plan[caom2_element][1]:
                        self._plan[caom2_element][1].insert(0, ttype_attribute)
                else:
                    raise AttributeError('No TTYPE attributes associated ' 'with keyword {}'.format(caom2_element))
            else:
                self._plan[caom2_element] = ('BINTABLE', [ttype_attribute], None)

    def set_default(self, caom2_element, default, extension=0):
        """
        Sets the default value of a caom2 element that is associated with attributes. If the element does not exist
        or does not have a list of associated attributes, default is set as the associated value of the element.

        If set is called for the same caom2_element after this, the default value will be reset to None.

        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param default: default value
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            if caom2_element in self._extensions[extension] and isinstance(
                self._extensions[extension][caom2_element], tuple
            ):
                self._extensions[extension][caom2_element] = (self._extensions[extension][caom2_element][0], default)
            else:
                # default is the only value
                self._extensions[extension][caom2_element] = default
        else:
            if (caom2_element in self._plan) and isinstance(self._plan[caom2_element], tuple):
                self._plan[caom2_element] = (self._plan[caom2_element][0], default)
            else:
                # override the value
                self._plan[caom2_element] = default

    def delete(self, caom2_element, extension=0):
        """
        Deletes an element from the blueprint
        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :raises exceptions if the element or extension not found
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise ValueError('Extension {} not configured in blueprint'.format(extension))
            if caom2_element in self._extensions[extension]:
                del self._extensions[extension][caom2_element]
            if len(self._extensions[extension]) == 0:
                del self._extensions[extension]
        else:
            if caom2_element in self._plan:
                del self._plan[caom2_element]

    def clear(self, caom2_element, extension=0):
        """
        Clears the value for an element in the blueprint by resetting it to an empty list with no default.

        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :raises exceptions if the element or extension not found
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise ValueError('Extension {} not configured in blueprint'.format(extension))
            if caom2_element in self._extensions[extension]:
                self._extensions[extension][caom2_element] = ([], None)
        else:
            if caom2_element in self._plan:
                self._plan[caom2_element] = ([], None)

    def _get(self, caom2_element, extension=0):
        """
        Returns the source associated with a CAOM2 element
        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :return: Tuple of the form (list_of_associated_attributes, default_value) OR the actual value associated
            with the CAOM2 element
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            if (extension in self._extensions) and (caom2_element in self._extensions[extension]):
                return self._extensions[extension][caom2_element]

        # look in the minimal plan
        if caom2_element not in self._plan:
            return None
        else:
            return self._plan[caom2_element]

    def has_chunk(self, extension):
        """What does the plan say about creating chunks for an extension?

        :return True if there should be a chunk to go along with a part
        """
        value = ''
        if extension is not None and extension in self._extensions:
            if 'Chunk' in self._extensions[extension]:
                value = self._extensions[extension]['Chunk']
        elif 'Chunk' in self._plan:
            if (extension is not None and extension == 0) or (extension is None):
                value = self._plan['Chunk']
        return not value == '{ignore}'

    @staticmethod
    def is_table(value):
        """Hide the blueprint structure from clients - they shouldn't need to know that a value of type tuple
        requires special processing."""
        return ObsBlueprint.needs_lookup(value) and value[0] == 'BINTABLE'

    @staticmethod
    def is_function(value):
        """
        Check if a blueprint value has Python 'function' syntax. The "'/' not in value" clause excludes strings with
        syntax that enables addressing HDF5 arrays.

        :return: True if the value is the name of a function to be executed, False, otherwise
        """
        return (
            not ObsBlueprint.needs_lookup(value)
            and isinstance(value, str)
            and isinstance(value, str)
            and '(' in value
            and ')' in value
            and '/' not in value
        )

    @staticmethod
    def has_default_value(value):
        """"""
        return isinstance(value, tuple) and value[1]

    @staticmethod
    def has_no_value(value):
        """If functions return None, try not to update the WCS with this value."""
        return value is None or (isinstance(value, str) and 'None' in value.strip())

    @staticmethod
    def needs_lookup(value):
        """Hide the blueprint structure from clients - they shouldn't need to know that a value of type tuple
        requires special processing."""
        return isinstance(value, tuple)

    def get_configed_axes_count(self):
        """:return how many axes have been configured to read from WCS"""
        configed_axes = 0
        if self._pos_axes_configed:
            configed_axes += 2
        if self._energy_axis_configed:
            configed_axes += 1
        if self._time_axis_configed:
            configed_axes += 1
        if self._polarization_axis_configed:
            configed_axes += 1
        if self._obs_axis_configed:
            configed_axes += 1
        if self._custom_axis_configed:
            configed_axes += 1
        return configed_axes

    @property
    def update(self):
        return self._update

    @update.setter
    def update(self, value):
        self._update = value


class Hdf5ObsBlueprint(ObsBlueprint):
    """
    Class that specializes the CAOM2 Observation construction based on HDF5 file content.

    The blueprint designates the source of each of these attributes as either HDF5 Dataset or Group values.
    Specific or default values may also be indicated in the same fashion os for an ObsBlueprint. The blueprint can
    be checked by simply displaying it.

    HDF5-specific example:
    # create a blueprint and customize it
    ob = Hdf5ObsBlueprint(position_axes=(1, 2)

    # lookup value starting with // means rooted at base of the hdf5 file
    ob.add_attribute('Observation.target.name', '//header/object/obj_id')

    # lookup value starting with / means rooted at the base of the "find_roots_here" parameter for Hdf5Parser
    # (integer) means return only the value with the index of "integer" from a list
    ob.add_attribute('Chunk.position.axis.function.refCoord.coord1.pix', '/header/wcs/crpix(0)')

    # (integer:integer) means return only the value with the index of "integer" from a list, followed by "integer"
    # from the list in the list
    ob.add_attribute('Chunk.position.axis.function.cd11', '/header/wcs/cd(0:0)')
    print(ob)

    """

    def __init__(
        self,
        position_axes=None,
        energy_axis=None,
        polarization_axis=None,
        time_axis=None,
        obs_axis=None,
        custom_axis=None,
        module=None,
        update=True,
        instantiated_class=None,
    ):
        """
        There are no sensible/known HDF5 defaults for WCS construction, so default to ensuring the blueprint
        executes with mostly values of None.

        Use the attribute _wcs_std, so that the list of WCS keywords used as input is known.
        """
        super().__init__(
            position_axes,
            energy_axis,
            polarization_axis,
            time_axis,
            obs_axis,
            custom_axis,
            module,
            update,
            instantiated_class,
        )
        tmp = {
            'Observation.algorithm.name': ([], 'exposure'),
            'Plane.calibrationLevel': ([], CalibrationLevel.RAW_STANDARD),
            'Plane.dataProductType': ([], DataProductType.IMAGE),
            'Artifact.releaseType': ([], ReleaseType.DATA),
            'Chunk': 'include',
        }
        # using the tmp to make sure that the keywords are valid
        for key in tmp:
            self.set(key, tmp[key])

    def configure_custom_axis(self, axis, override=True):
        """
        Set the expected custom keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the custom axis.
        :param override: Set to False when reading from a file.
        """
        if self._custom_axis_configed:
            self.logger.debug('Attempt to configure already-configured custom axis.')
            return

        if override:
            self.set('Chunk.custom.axis.axis.ctype', ([], None))
            self.set('Chunk.custom.axis.axis.cunit', ([], None))
            self.set('Chunk.custom.axis.function.naxis', ([], 1))
            self.set('Chunk.custom.axis.function.delta', ([], None))
            self.set('Chunk.custom.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.custom.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.custom.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.custom.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.custom.axis.function.naxis'] = ''
        self._wcs_std['Chunk.custom.axis.function.delta'] = ''
        self._wcs_std['Chunk.custom.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.custom.axis.function.refCoord.val'] = ''
        self._custom_axis_configed = True

    def configure_position_axes(self, axes, override=True):
        """
        Set the expected spatial keywords by indices in the blueprint and the wcs_std lookup.

        :param axes: The index expected for the position axes.
        :param override: Set to False when reading from a file.
        """
        if self._pos_axes_configed:
            self.logger.debug('Attempt to configure already-configured position axes.')
            return

        if override:
            self.set('Chunk.position.coordsys', ([], None))
            self.set('Chunk.position.equinox', ([], None))
            self.set('Chunk.position.axis.axis1.ctype', ([], None))
            self.set('Chunk.position.axis.axis1.cunit', ([], None))
            self.set('Chunk.position.axis.axis2.ctype', ([], None))
            self.set('Chunk.position.axis.axis2.cunit', ([], None))
            self.set('Chunk.position.axis.error1.syser', ([], None))
            self.set('Chunk.position.axis.error1.rnder', ([], None))
            self.set('Chunk.position.axis.error2.syser', ([], None))
            self.set('Chunk.position.axis.error2.rnder', ([], None))
            self.set('Chunk.position.axis.function.cd11', ([], None))
            self.set('Chunk.position.axis.function.cd12', ([], None))
            self.set('Chunk.position.axis.function.cd21', ([], None))
            self.set('Chunk.position.axis.function.cd22', ([], None))
            self.set('Chunk.position.axis.function.dimension.naxis1', ([], 1))
            self.set('Chunk.position.axis.function.dimension.naxis2', ([], 1))
            self.set('Chunk.position.axis.function.refCoord.coord1.pix', ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.val', ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.pix', ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.val', ([], None))

        self._wcs_std['Chunk.position.coordsys'] = ''
        self._wcs_std['Chunk.position.equinox'] = ''

        self._wcs_std['Chunk.position.axis.axis1.ctype'] = ''
        self._wcs_std['Chunk.position.axis.axis1.cunit'] = ''
        self._wcs_std['Chunk.position.axis.axis2.ctype'] = ''
        self._wcs_std['Chunk.position.axis.axis2.cunit'] = ''
        self._wcs_std['Chunk.position.axis.error1.syser'] = ''
        self._wcs_std['Chunk.position.axis.error1.rnder'] = ''
        self._wcs_std['Chunk.position.axis.error2.syser'] = ''
        self._wcs_std['Chunk.position.axis.error2.rnder'] = ''
        self._wcs_std['Chunk.position.axis.function.cd11'] = ''
        self._wcs_std['Chunk.position.axis.function.cd12'] = ''
        self._wcs_std['Chunk.position.axis.function.cd21'] = ''
        self._wcs_std['Chunk.position.axis.function.cd22'] = ''
        self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = ''
        self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] = ''

        self._pos_axes_configed = True

    def configure_energy_axis(self, axis, override=True):
        """
        :param axis: The index expected for the energy axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._energy_axis_configed:
            self.logger.debug('Attempt to configure already-configured energy axis.')
            return

        if override:
            self.set('Chunk.energy.specsys', ([], None))
            self.set('Chunk.energy.ssysobs', ([], None))
            self.set('Chunk.energy.restfrq', ([], None))
            self.set('Chunk.energy.restwav', ([], None))
            self.set('Chunk.energy.velosys', ([], None))
            self.set('Chunk.energy.zsource', ([], None))
            self.set('Chunk.energy.ssyssrc', ([], None))
            self.set('Chunk.energy.velang', ([], None))

            self.set('Chunk.energy.bandpassName', ([], None))
            self.set('Chunk.energy.resolvingPower', ([], None))

            self.set('Chunk.energy.axis.axis.ctype', ([], None))
            self.set('Chunk.energy.axis.axis.cunit', ([], None))
            self.set('Chunk.energy.axis.error.syser', ([], None))
            self.set('Chunk.energy.axis.error.rnder', ([], None))
            self.set('Chunk.energy.axis.function.naxis', ([], 1))
            self.set('Chunk.energy.axis.function.delta', ([], None))
            self.set('Chunk.energy.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.energy.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.energy.specsys'] = ''
        self._wcs_std['Chunk.energy.ssysobs'] = ''
        self._wcs_std['Chunk.energy.restfrq'] = ''
        self._wcs_std['Chunk.energy.restwav'] = ''
        self._wcs_std['Chunk.energy.velosys'] = ''
        self._wcs_std['Chunk.energy.zsource'] = ''
        self._wcs_std['Chunk.energy.ssyssrc'] = ''
        self._wcs_std['Chunk.energy.velang'] = ''

        self._wcs_std['Chunk.energy.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.energy.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.energy.axis.error.syser'] = ''
        self._wcs_std['Chunk.energy.axis.error.rnder'] = ''
        self._wcs_std['Chunk.energy.axis.function.naxis'] = ''
        self._wcs_std['Chunk.energy.axis.function.delta'] = ''
        self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = ''
        self._energy_axis_configed = True

    def configure_polarization_axis(self, axis, override=True):
        """
        Set the expected polarization keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the polarization axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._polarization_axis_configed:
            self.logger.debug('Attempt to configure already-configured polarization axis.')
            return

        if override:
            # STOKES is the only value allowed for PolarizationWCS ctype.
            self.set('Chunk.polarization.axis.axis.ctype', ([], 'STOKES'))
            self.set('Chunk.polarization.axis.axis.cunit', ([], None))
            self.set('Chunk.polarization.axis.function.naxis', ([], 1))
            self.set('Chunk.polarization.axis.function.delta', ([], None))
            self.set('Chunk.polarization.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.polarization.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.polarization.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.polarization.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.polarization.axis.function.naxis'] = ''
        self._wcs_std['Chunk.polarization.axis.function.delta'] = ''
        self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = ''

        self._polarization_axis_configed = True

    def configure_observable_axis(self, axis, override=True):
        """
        Set the expected observable keywords by index in the blueprint and the wcs_std lookup.
        Note: observable axis is not a standard WCS and it's not used by astropy.wcs so, arguably, it can be
        removed. It is here for now for consistency purposes.
        :param axis: The index expected for the observable axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._obs_axis_configed:
            self.logger.debug('Attempt to configure already-configured observable axis.')
            return

        if override:
            self.set('Chunk.observable.axis.axis.ctype', ([], None))
            self.set('Chunk.observable.axis.axis.cunit', ([], None))
            self.set('Chunk.observable.axis.function.refCoord.pix', ([], None))

        self._wcs_std['Chunk.observable.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.observable.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.observable.axis.function.refCoord.pix'] = ''

        self._obs_axis_configed = True

    def configure_time_axis(self, axis, override=True):
        """
        Set the expected time keywords by index in the blueprint and the wcs_std lookup.

        :param axis: The index expected for the time axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._time_axis_configed:
            self.logger.debug('Attempt to configure already-configured time axis.')
            return

        if override:
            self.set('Chunk.time.exposure', ([], None))
            self.set('Chunk.time.timesys', ([], None))
            self.set('Chunk.time.trefpos', ([], None))
            self.set('Chunk.time.mjdref', ([], None))
            self.set('Chunk.time.resolution', ([], None))
            self.set('Chunk.time.axis.axis.ctype', ([], None))
            self.set('Chunk.time.axis.axis.cunit', ([], None))
            self.set('Chunk.time.axis.error.syser', ([], None))
            self.set('Chunk.time.axis.error.rnder', ([], None))
            self.set('Chunk.time.axis.function.naxis', ([], 1))
            self.set('Chunk.time.axis.function.delta', ([], None))
            self.set('Chunk.time.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.time.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.time.exposure'] = ''
        self._wcs_std['Chunk.time.resolution'] = ''
        self._wcs_std['Chunk.time.timesys'] = ''
        self._wcs_std['Chunk.time.trefpos'] = ''
        self._wcs_std['Chunk.time.mjdref'] = ''

        self._wcs_std['Chunk.time.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.time.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.time.axis.error.syser'] = ''
        self._wcs_std['Chunk.time.axis.error.rnder'] = ''
        self._wcs_std['Chunk.time.axis.function.naxis'] = ''
        self._wcs_std['Chunk.time.axis.function.delta'] = ''
        self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.time.axis.function.refCoord.val'] = ''

        self._time_axis_configed = True

    def set(self, caom2_element, value, extension=0):
        """
        Sets the value associated with an element in the CAOM2 model. Value cannot be a tuple.
        :param caom2_element: name CAOM2 element (as in ObsBlueprint.CAOM2_ELEMEMTS)
        :param value: new value of the CAOM2 element
        :param extension: extension number (used only for Chunk elements)
        """
        if hasattr(value, 'decode'):
            value = value.decode('utf-8')
        super().set(caom2_element, value, extension)

    def _guess_axis_info(self):
        self._guess_axis_info_from_plan()


def _to_float(value):
    return float(value) if value is not None else None


def _to_int(value):
    return int(value) if value is not None else None


def _to_int_32(value):
    if value is None:
        return None
    elif isinstance(value, str):
        return int_32(value)
    else:
        return value


def _to_str(value):
    if value is None or str(value).strip() == '':
        result = None
    else:
        result = str(value).strip()
    return result
