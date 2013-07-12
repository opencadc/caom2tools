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

"""defines the SpectralWCS class"""

from caom2.caom2_energy_transition import EnergyTransition
from caom2_coord_axis1d import CoordAxis1D
from caom2.util import caom2_util as util
from caom2.caom2_object import Caom2Object


class SpectralWCS(Caom2Object):
    """A transformation that maps pixel coordinates to spectral ones.

    Note that a 2D image has implicit 'spectral' (and temporal)
    dimension that is one pixel in size.  Basically pixels are really
    multidimensional voxels.

    Due to FITS standards this pixel starts at 0.5 and runs to 1.5
    with the centre being at 1.0

    """

    def __init__(self,
                 axis,
                 specsys,
                 ssysobs=None,
                 ssyssrc=None,
                 restfrq=None,
                 restwav=None,
                 velosys=None,
                 zsource=None,
                 velang=None,
                 bandpass_name=None,
                 transition=None,
                 resolving_power=None
                 ):
        """The SpectralWCS can be defined in a number of different ways, to
        define one you must provide a CoordAxis1D object that maps the
        pixels to WCS values and a reference specsys.  After that the
        user can add what ever parts seam useful.  More info is more
        helpful for searching.

        """

        self.axis = axis
        self.specsys = specsys
        self.ssysobs = ssysobs
        self.ssyssrc = ssyssrc
        self.restfrq = restfrq
        self.restwav = restwav
        self.velosys = velosys
        self.zsource = zsource
        self.velang = velang
        self.bandpass_name = bandpass_name
        self.transition = transition
        self.resolving_power = resolving_power

    @property
    def axis(self):
        """A 1D coordinate axis object that contains the pix/wcs
        transformation values.

        eg.  CoordAxis1D(Axis('wave','flux'),...)

        """
        return self._axis

    @axis.setter
    def axis(self, value):
        util.typeCheck(value, CoordAxis1D, 'axis', override=False)
        self._axis = value

    @property
    def specsys(self):
        """describes the reference frame in use for the spectral-axis
        coordinate(s).

        eg. BARYCENT

        type: str
        """
        return self._specsys

    @specsys.setter
    def specsys(self, value):
        util.typeCheck(value, str, 'specsys', override=False)
        self._specsys = value

    @property
    def ssysobs(self):
        """describes the spectral reference frame that is constant over the
        range of the non-spectral world coordinates

        For example, for a large image the the wavelength at the edges
        is different from the centres.  This refernce frame is one where they
        are not different.

        Nominally 'TOPOCENT'

        type: str
        """
        return self._ssysobs

    @ssysobs.setter
    def ssysobs(self, value):
        util.typeCheck(value, str, 'ssysobs')
        self._ssysobs = value

    @property
    def ssyssrc(self):
        """The reference frame in which zsource is expressed.

        eg. BARYCENT
        type: string
        """
        return self._ssyssrc

    @ssyssrc.setter
    def ssyssrc(self, value):
        util.typeCheck(value, str, 'ssyssrc')
        self._ssyssrc = value

    @property
    def restfrq(self):
        """The frequency of the spectal feature being observed.

        unit: Hz
        type: float
        """
        return self._restfrq

    @restfrq.setter
    def restfrq(self, value):
        util.typeCheck(value, float, 'restfrq')
        self._restfrq = value

    @property
    def restwav(self):
        """The wavelength of spectral feature being observed,
        not the wavelength observed but the wavelength of the
        feature when at rest..

        unit: m
        type: float
        """
        return self._restwav

    @restwav.setter
    def restwav(self, value):
        util.typeCheck(value, float, 'restwav')
        self._restwav = value

    @property
    def velosys(self):
        """Relative radial velocity between the observer and the selected
        standard of rest in the direction of the celestial reference
        coordinate.

        eg. 26000 m/s


        unit: m/s
        type: float
        """
        return self._velosys

    @velosys.setter
    def velosys(self, value):
        util.typeCheck(value, float, 'velosys')
        self._velosys = value

    @property
    def zsource(self):
        """The redshift of the source emitting the photons.

        almost always None

        unit: z
        type: float
        """
        return self._zsource

    @zsource.setter
    def zsource(self, value):
        util.typeCheck(value, float, 'zsource')
        util.valueCheck(value, -0.5, 1200, 'zsource')
        self._zsource = value

    @property
    def velang(self):
        """I don't know what this is... angle of the velocity ??? """
        return self._velang

    @velang.setter
    def velang(self, value):
        util.typeCheck(value, float, 'velang')
        self._velang = value

    @property
    def bandpass_name(self):
        """string the represent the bandpass of the observation.

        eg. r'
        type: str
        """
        return self._bandpass_name

    @bandpass_name.setter
    def bandpass_name(self, value):
        util.typeCheck(value, str, 'bandpass_name')
        self._bandpass_name = value

    @property
    def transition(self):
        """which molecular transition has been observed.

        type: EnergyTransition object (see caom2.EnergyTransition for help)
        """
        return self._transition

    @transition.setter
    def transition(self, value):
        util.typeCheck(value, EnergyTransition, "transition")
        self._transition = value

    @property
    def resolving_power(self):
        """The R value of the spectal coverage.

        Normally this is something like dlamda/lamda

        unit:  RATIO
        type: float
        """
        return self._resolving_power

    @resolving_power.setter
    def resolving_power(self, value):
        util.typeCheck(value, float, 'resolving_power')
        util.valueCheck(value, 0, 1E7, 'resolving_power')
        self._resolving_power = value
