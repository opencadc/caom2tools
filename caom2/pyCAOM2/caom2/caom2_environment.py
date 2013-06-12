#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
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

"""
Defines the caom2.Environment class, populates the caom2_environment values
"""

import util.caom2_util as util


class Environment(object):
    """A CAOM2 Environment Object.

    This object contains the various values that can be set in the environment
    table entry. Normally each Observation object will have an associate
    Environment Object."""

    def __init__(self):
        """
        Initializes an Environment instance

        Arguments:
        """
        self._seeing = None
        self._humidity = None
        self._elevation = None
        self._tau = None
        self._wavelength_tau = None
        self._ambient_temp = None
        self._photometric = None
        self._print_attributes = ['seeing', 'photometric',
                                  'tau', 'wavelength_tau', 'elevation',
                                  'ambient_temp', 'humidity']

    # Properties
    @property
    def seeing(self):
        """atmospheric seeing (FWHM) in arcsec.

        units: arcseconds
        """
        return self._seeing

    @seeing.setter
    def seeing(self, value):
        util.typeCheck(value, float, 'seeing')
        util.valueCheck(value, 0, 360 * 3600.0, 'seeing')
        self._seeing = value

    @property
    def humidity(self):
        """Relative humidity, expressed as a fraction between 0 and 1.

        units: fraction
        """
        return self._humidity

    @humidity.setter
    def humidity(self, value):
        # set the humidity to value, which  must be +ve fraction less than 1
        util.typeCheck(value, float, 'humidity')
        util.valueCheck(value, 0, 1, 'humidity')
        self._humidity = value

    @property
    def elevation(self):
        """ Elevation above horizon (0 to 90 deg) at which tau was measured.

        units: deg
        """
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        util.typeCheck(value, float, 'elevation')
        util.valueCheck(value, 0, 90, 'elevation')
        self._elevation = value

    @property
    def tau(self):
        """The tau at zennith at the time of observation.

        units:  fraction

        This tau can be used, in combination with the elevation,
        to determine the actual tau."""
        return self._tau

    @tau.setter
    def tau(self, value):
        util.typeCheck(value, float, 'tau')
        util.valueCheck(value, 0, 1, 'tau')
        self._tau = value

    @property
    def wavelength_tau(self):
        """Wavelength at which tau was measured
        (normally 225GHz converted to wavelength).

        units: meters

        """
        return self._wavelength_tau

    @wavelength_tau.setter
    def wavelength_tau(self, value):
        util.typeCheck(value, float, 'wavelength_tau')
        util.valueCheck(value, 0, 1E3, 'wavelength_tau')
        self._wavelength_tau = value

    @property
    def ambient_temp(self):
        """The ambient air temperature at time of observation.

        unit: Celsius degrees
        """
        return self._ambient_temp

    @ambient_temp.setter
    def ambient_temp(self, value):
        util.typeCheck(value, float, 'ambient_temp')
        util.valueCheck(value, -100, 100, 'ambient_temp')
        self._ambient_temp = value

    @property
    def photometric(self):
        """A boolean flag (True/False) indicating if
        the observational conditions were photometric."""
        return self._photometric

    @photometric.setter
    def photometric(self, value):
        util.typeCheck(value, bool, 'photometric')
        self._photometric = value
