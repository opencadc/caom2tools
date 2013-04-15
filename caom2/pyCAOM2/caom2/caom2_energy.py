# -*- coding: latin-1 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2010.                            (c) 2010.
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
#***********************************************************************
#

"""Defines the caom2.Energy object"""

from caom2_object import Caom2Object
from caom2_entity import AbstractCaom2Entity
from caom2_energy_transition import EnergyTransition
from caom2_enums import EnergyBand
from types.caom2_interval import Interval
import util.caom2_util as util

class Energy(Caom2Object):
    """ Energy """

    def __init__(self):
        """
        Initialize an Energy instance.

        Arguments:
        None
        """
        self._value = None
        self._bounds = None
        self._dimension = None
        self._resolving_power = None
        self._sample_size = None
        self._bandpass_name = None
        self._em_band = None
        self._transition = None

        self._print_attributes = ['value',
                                  'dimension',
                                  'resolving_power',
                                  'sample_size',
                                  'bandpass_name',
                                  'em_band',
                                  'transition',
                                  'bounds']


    # Properties

    @property
    def value(self):
        """ Energy value """
        return self._value

    @value.setter
    def value(self, value):
        if value is not None:
            assert isinstance(value, float), (
                    "energy value is not a float: {0}".format(value))
        self._value = value

    @property
    def bounds(self):
        """ Energy bounds """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        if value is not None:
            assert isinstance(value, Interval), (
                    "energy bounds is not an Interval: {0}".format(value))
        self._bounds = value

    @property
    def dimension(self):
        """DIMENSION (NUMBER OF PIXELS) ALONG ENERGY AXIS."""
        return self._dimension

    @dimension.setter
    def dimension(self, value):
        if value is not None:
            assert isinstance(value, long), (
                    "energy dimension is not a long: {0}".format(value))
        self._dimension = value

    @property
    def resolving_power(self):
        """ Resolving power """
        return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        if value is not None:
            assert isinstance(value, float), (
                    "resolving power is not a float: {0}".format(value))
        self._resolving_power = value

    @property
    def sample_size(self):
        """ Sample size """
        return self._sample_size

    @sample_size.setter
    def sample_size(self, value):
        if value is not None:
            assert isinstance(value, float), (
                    "sample size is not a float: {0}".format(value))
        self._sample_size = value

    @property
    def bandpass_name(self):
        """ Bandpass name """
        return self._bandpass_name

    @bandpass_name.setter
    def bandpass_name(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "bandpass name is not str: {0}".format(value))
        self._bandpass_name = value

    @property
    def em_band(self):
        """ EM Band """
        return self._em_band

    @em_band.setter
    def em_band(self, value):
        if value is not None:
            assert isinstance(value, EnergyBand), (
                    "em_Band is not an EnergyBand: {0}".format(value))
        self._em_band = value

    @property
    def transition(self):
        """ Energy transition """
        return self._transition

    @transition.setter
    def transition(self, value):
        if value is not None:
            assert isinstance(value, EnergyTransition), (
                    "transition is not an EnergyTransition: {0}".format(value))
        self._transition = value
