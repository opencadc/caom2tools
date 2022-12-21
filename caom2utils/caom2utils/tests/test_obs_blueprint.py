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

import pytest
import sys


# @pytest.mark.skip('')
def test_obs_blueprint():
    # test the CAOM2_ELEMENTS property
    assert ObsBlueprint.CAOM2_ELEMENTS == ObsBlueprint._CAOM2_ELEMENTS

    # updating the CAOM2_ELEMENTS property should not update the original
    elems = ObsBlueprint.CAOM2_ELEMENTS
    elems = elems[2:]
    assert elems != ObsBlueprint.CAOM2_ELEMENTS

    # default config (one entry per row...)
    assert str(ObsBlueprint()).count('\n') == 24
    print(ObsBlueprint())

    # default config with WCS info
    assert str(ObsBlueprint(position_axes=(1, 2), energy_axis=3,
                            polarization_axis=4, time_axis=5,
                            obs_axis=6, custom_axis=7)).count('\n') == 90

    ob = ObsBlueprint()
    ob.configure_position_axes(axes=(1, 2))
    ob.configure_energy_axis(axis=3)
    ob.configure_polarization_axis(axis=4)
    ob.configure_time_axis(axis=5)
    ob.configure_observable_axis(axis=6)
    ob.configure_custom_axis(axis=7)

    # test that configuring something that's already configured doesn't break
    # anything
    ob.configure_position_axes(axes=(1, 2))
    ob.configure_energy_axis(axis=3)
    ob.configure_polarization_axis(axis=4)
    ob.configure_time_axis(axis=5)
    ob.configure_observable_axis(axis=6)
    ob.configure_custom_axis(axis=7)

    # set attribute
    ob.set('Observation.instrument.name', 'NIRI')
    assert ob._plan['Observation.instrument.name'] == 'NIRI'
    assert 'Observation.instrument.name = NIRI' in str(ob)

    # set default
    ob.clear('Observation.instrument.keywords')
    ob.add_attribute('Observation.instrument.keywords', 'INSTMODE')
    assert "Observation.instrument.keywords = ['INSTMODE'], default = None" \
           in str(ob)
    ob.set_default('Observation.instrument.keywords', 'TEST')
    assert ob._plan['Observation.instrument.keywords'][1] == 'TEST'
    assert ob._plan['Observation.instrument.keywords'][0] == ['INSTMODE']
    assert "Observation.instrument.keywords = ['INSTMODE'], default = TEST" \
           in str(ob)

    # set fits attribute
    ob.add_attribute('Observation.proposal.id', 'PROP')
    ob.add_attribute('Observation.proposal.id', 'PROP2')
    ob.set_default('Observation.proposal.id', 'NOPROP')
    assert ob._plan['Observation.proposal.id'][0] == ['PROP2', 'PROP', 'RUNID']
    assert ob._plan['Observation.proposal.id'][1] == 'NOPROP'
    assert ("Observation.proposal.id = ['PROP2', 'PROP', 'RUNID'], "
                "default = NOPROP") in str(ob)

    # set in extension
    ob.set('Chunk.energy.velang', 33, extension=1)
    extension1_str = str(ob)[str(ob).index('extension 1'):]
    assert 'Chunk.energy.velang = 33' in extension1_str

    # set fits attribute in extension
    ob.add_attribute('Chunk.energy.axis.axis.ctype', 'MYCTYPE',
                          extension=1)
    ob.add_attribute('Chunk.energy.axis.axis.ctype', 'MYCTYPE2',
                          extension=1)
    ob.set_default('Chunk.energy.axis.axis.ctype', 'NOCTYPE', extension=1)
    extension1_str = str(ob)[str(ob).index('extension 1'):]
    assert ("Chunk.energy.axis.axis.ctype = ['MYCTYPE2', 'MYCTYPE'], "
            "default = NOCTYPE") in extension1_str

    # set in a different extension
    ob.set('Chunk.energy.velang', 44, extension=2)
    extension2_str = str(ob)[str(ob).index('extension 2'):]
    assert 'Chunk.energy.velang = 44' in extension2_str

    # test get
    assert ob._get('Observation.instrument.name') == 'NIRI'
    assert ob._get('Observation.instrument.keywords')[0] == ['INSTMODE']
    assert ob._get('Observation.instrument.keywords')[1] == 'TEST'
    assert ob._get('Chunk.energy.velang', extension=2) == 44
    assert ob._get('Chunk.energy.velang', extension=1) == 33
    assert ob._get('Chunk.energy.axis.axis.ctype', extension=1)[0] ==\
        ['MYCTYPE2', 'MYCTYPE']
    assert ob._get('Chunk.energy.axis.axis.ctype', extension=1)[1] == 'NOCTYPE'
    # test get when keyword not present in extension and the default is used
    assert ob._get('Chunk.energy.specsys', extension=33)[0] == ['SPECSYS']
    assert ob._get('Chunk.energy.specsys', extension=44)[1] is None
    # get element not set
    assert 'Observation.target.redshift' not in ob._plan
    assert ob._get('Observation.target.redshift') is None

    # bintable :)
    ob.add_table_attribute('CompositeObservation.members', 'FICS', extension=0)
    result = ob._get('CompositeObservation.members', extension=0)
    assert result is not None
    assert len(result) == 3, len(result)

    # delete attributes
    assert len(ob._plan) != 0
    assert len(ob._extensions) != 0
    for i in ObsBlueprint.CAOM2_ELEMENTS:
        ob.delete(i)
    assert len(ob._plan) == 0

    # clear
    ob.clear('Chunk.energy.axis.axis.ctype', extension=1)

    # delete attributes from extensions
    ob.delete('Chunk.energy.velang', extension=1)
    ob.delete('Chunk.energy.axis.axis.ctype', extension=1)
    ob.delete('CompositeObservation.members', extension=0)
    ob.delete('Chunk.energy.velang', extension=2)
    assert len(ob._extensions) == 0

    # set defaults in extension
    ob.set_default('Chunk.energy.axis.axis.ctype', 'NOCTYPE', extension=3)
    extension3_str = str(ob)[str(ob).index('extension 3'):]
    assert "Chunk.energy.axis.axis.ctype = NOCTYPE" in extension3_str
    assert len(ob._extensions) == 1

    # testing error cases
    ob = ObsBlueprint()
    ob.configure_position_axes(axes=(1, 2))
    ob.configure_energy_axis(axis=3)
    ob.configure_polarization_axis(axis=4)
    ob.configure_time_axis(axis=5)

    # non CAOM2 element name
    with pytest.raises(KeyError):
        ob.set('Nonexistent', 33)
    with pytest.raises(KeyError):
        ob.add_attribute('Nonexistent', 33)
    with pytest.raises(KeyError):
        ob.set_default('Nonexistent', 33)
    with pytest.raises(KeyError):
        ob._get('Nonexistent')
    with pytest.raises(KeyError):
        ob.delete('Nonexistent')
    with pytest.raises(KeyError):
        ob.clear('Nonexistent')

    # repeat with extension specified
    with pytest.raises(KeyError):
        ob.set('Chunk.Nonexistent', 33, extension=1)
    with pytest.raises(KeyError):
        ob.add_attribute('Chunk.Nonexistent', 33, extension=1)
    with pytest.raises(KeyError):
        ob.set_default('Chunk.Nonexistent', 33, extension=1)
    with pytest.raises(KeyError):
        ob._get('Nonexistent', extension=33)
    with pytest.raises(KeyError):
        ob.delete('Chunk.Nonexistent', extension=1)
    with pytest.raises(KeyError):
        ob.clear('Chunk.Nonexistent', extension=1)

    # CAOM2 element not Chunk with extension specified
    with pytest.raises(ValueError):
        ob.set('Observation.observationID', 33, extension=1)
    with pytest.raises(ValueError):
        ob.set_default('Observation.observationID', 33, extension=1)
    with pytest.raises(ValueError):
        ob.add_attribute('Observation.observationID', 'AA', extension=1)
    with pytest.raises(ValueError):
        ob.delete('Observation.observationID', extension=1)
    with pytest.raises(ValueError):
        ob.clear('Observation.observationID', extension=1)

    # CAOM2 element Chunk with extension
    with pytest.raises(ValueError):
        ob.clear('Chunk.energy.axis.axis.cunit', extension=33)

    # add FITS attribute to element that does not exist
    assert 'Chunk.energy.transition' not in ob._plan
    ob.set('Chunk.energy.transition', 'Name')
    with pytest.raises(AttributeError):
        ob.add_attribute('Chunk.energy.transition', 'BP')

    # call set_fits_attribute with argument other than list
    with pytest.raises(AttributeError):
        ob.add_attribute('Chunk.energy.transition', 33)

    # delete element from a non-existent extension
    with pytest.raises(ValueError):
        ob.delete('Chunk.energy.transition', extension=66)

    # adding the same thing twice does nothing - the test values are defaults
    result = ob._get('Observation.metaRelease')
    initial_result_length = (len(result[0]))
    ob.add_attribute('Observation.metaRelease', 'DATE-OBS')
    result = ob._get('Observation.metaRelease')
    add_result_length = (len(result[0]))
    assert initial_result_length == add_result_length
    # in an extension
    result = ob._get('Chunk.energy.specsys', extension=1)
    initial_result_length = (len(result[0]))
    ob.add_attribute('Chunk.energy.specsys', 'SPECSYS')
    result = ob._get('Chunk.energy.specsys', extension=1)
    add_result_length = (len(result[0]))
    assert initial_result_length == add_result_length, result


def test_load_from_file_configure():
    ob = ObsBlueprint()
    assert not ob._pos_axes_configed, \
        'Failure to initialize configure_position_axes'
    assert not ob._energy_axis_configed, \
        'Failure to initialize configure_energy_axis'
    assert not ob._custom_axis_configed, 'custom config'
    assert not ob._obs_axis_configed, 'obs config'
    assert not ob._polarization_axis_configed, 'pol config'
    assert not ob._time_axis_configed, 'time config'
    ob.add_attribute('Chunk.position.axis.axis1.ctype', 'CTYPE1')
    ob.add_attribute('Chunk.position.axis.axis2.ctype', 'CTYPE2')
    ob.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    ob._guess_axis_info()
    assert ob._pos_axes_configed, 'Failure to call configure_position_axes'
    assert ob._energy_axis_configed, 'Failure to call configure_energy_axis'
    assert ob._wcs_std['Chunk.energy.axis.axis.ctype'] == 'CTYPE3', \
        ob._wcs_std['Chunk.energy.axis.axis.ctype']

    ob = ObsBlueprint()
    ob.add_attribute('Chunk.position.axis.axis1.ctype', 'CTYPE3')
    ob.add_attribute('Chunk.position.axis.axis2.ctype', 'CTYPE4')
    ob.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    ob._guess_axis_info()
    assert ob._pos_axes_configed, 'Failure to call configure_position_axes'
    assert ob._energy_axis_configed, 'Failure to call configure_energy_axis'
    assert ob._wcs_std['Chunk.energy.axis.axis.ctype'] == 'CTYPE1', \
        ob._wcs_std['Chunk.energy.axis.axis.ctype']

    ob = ObsBlueprint()
    ob.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.energy.axis.axis.ctype'] == 'CTYPE3', \
        ob._wcs_std['Chunk.energy.axis.axis.ctype']

    ob = ObsBlueprint()
    ob.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.polarization.axis.axis.ctype'] == 'CTYPE5', \
        ob._wcs_std['Chunk.polarization.axis.axis.ctype']
    assert ob._polarization_axis_configed, 'pol config'

    ob = ObsBlueprint()
    ob.set('Chunk.observable.axis.axis.ctype', 'COUNT')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.observable.axis.axis.ctype'] == 'CTYPE6', \
        ob._wcs_std['Chunk.observable.axis.axis.ctype']
    assert ob._obs_axis_configed, 'obs config'

    ob = ObsBlueprint()
    ob.set('Chunk.custom.axis.axis.ctype', 'FARDEP')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.custom.axis.axis.ctype'] == 'CTYPE7', \
        ob._wcs_std['Chunk.custom.axis.axis.ctype']
    assert ob._custom_axis_configed, 'custom config'

    ob = ObsBlueprint()
    ob.set('Chunk.time.axis.axis.ctype', 'TIME')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.time.axis.axis.ctype'] == 'CTYPE4', \
        ob._wcs_std['Chunk.time.axis.axis.ctype']
    assert ob._time_axis_configed, 'time config'

    # should get the position axes by default
    ob = ObsBlueprint()
    ob._guess_axis_info()
    assert ob._pos_axes_configed, 'pos config'
    assert not ob._energy_axis_configed, 'energy config'
    assert not ob._custom_axis_configed, 'custom config'
    assert not ob._obs_axis_configed, 'obs config'
    assert not ob._polarization_axis_configed, 'pol config'
    assert not ob._time_axis_configed, 'time config'

    # for the truly creative instrument scientist
    ob = ObsBlueprint()
    ob.add_attribute('Chunk.polarization.axis.axis.ctype', 'CTYPE1')
    ob.add_attribute('Chunk.custom.axis.axis.ctype', 'CTYPE2')
    ob.add_attribute('Chunk.position.axis.axis1.ctype', 'CTYPE3')
    ob.add_attribute('Chunk.position.axis.axis2.ctype', 'CTYPE4')
    ob.add_attribute('Chunk.time.axis.axis.ctype', 'CTYPE5')
    ob.add_attribute('Chunk.energy.axis.axis.ctype', 'CTYPE6')
    ob.add_attribute('Chunk.observable.axis.axis.ctype', 'CTYPE7')
    ob._guess_axis_info()
    assert ob._wcs_std['Chunk.polarization.axis.axis.ctype'] == 'CTYPE1', \
        ob._wcs_std['Chunk.polarization.axis.axis.ctype']
    assert ob._wcs_std['Chunk.custom.axis.axis.ctype'] == 'CTYPE2', \
        ob._wcs_std['Chunk.custom.axis.axis.ctype']
    assert ob._wcs_std['Chunk.position.axis.axis1.ctype'] == 'CTYPE3', \
        ob._wcs_std['Chunk.position.axis.axis1.ctype']
    assert ob._wcs_std['Chunk.position.axis.axis2.ctype'] == 'CTYPE4', \
        ob._wcs_std['Chunk.position.axis.axis2.ctype']
    assert ob._wcs_std['Chunk.time.axis.axis.ctype'] == 'CTYPE5', \
        ob._wcs_std['Chunk.time.axis.axis.ctype']
    assert ob._wcs_std['Chunk.energy.axis.axis.ctype'] == 'CTYPE6', \
        ob._wcs_std['Chunk.energy.axis.axis.ctype']
    assert ob._wcs_std['Chunk.observable.axis.axis.ctype'] == 'CTYPE7', \
        ob._wcs_std['Chunk.observable.axis.axis.ctype']

    with pytest.raises(ValueError):
        ob = ObsBlueprint()
        ob.add_attribute('Chunk.polarization.axis.axis.ctype', 'CTYPE1')
        ob._guess_axis_info()


def test_has_chunk():
    # the CFHT case
    ob = ObsBlueprint()
    ob.configure_position_axes((1, 2))
    ob.set('Chunk', '{ignore}')
    ob.set('Chunk.position.axis.axis1.ctype', 'RA---SIN', 1)
    assert not ob.has_chunk(0)
    assert ob.has_chunk(1)

    # the OMM case
    ob = ObsBlueprint()
    ob.configure_position_axes((1, 2))
    ob.set('Chunk', '{ignore}', 1)
    assert ob.has_chunk(0)
    assert not ob.has_chunk(1)


def test_ctor_failure():
    with pytest.raises(ValueError):
        ObsBlueprint(position_axes=(1, 2, 3))
