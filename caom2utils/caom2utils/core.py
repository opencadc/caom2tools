# -*- coding: utf-8 -*-
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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

""" Defines utilities for working with fits files """

ENERGY_KEYWORDS = [
    'FREQ',
    'ENER',
    'WAVN',
    'VRAD',
    'WAVE',
    'VOPT',
    'ZOPT',
    'AWAV',
    'VELO',
    'BETA']


from astropy.io import fits
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk
from caom2 import SpectralWCS, CoordAxis1D


def augment_artifact(artifact, file, collection=None):
    hdulist = fits.open(file)
    hdulist.close()
    parts = len(hdulist)
    print(repr(hdulist[0].header))

    if not artifact:
        assert not collection
        artifact = Artifact('ad:{}/{}'.format(collection, file),
                            ProductType.SCIENCE, ReleaseType.DATA) #TODO

    for p in range(parts):
        if str(p) not in artifact.parts.keys():
            artifact.parts.add(Part(str(p)))
        part = artifact.parts[str(p)]
        if not part.chunks:
            part.chunks.append(Chunk())
        chunk = part.chunks[0]
        augment_energy(chunk, hdulist[0].header)
        print('*******Chunk - {}'.format(chunk))


def augment_energy(chunk, header):
    specsys = header.get('SPECSYS', None)
    naxis = header.get('NAXIS', 0) or header.get('ZAXIS', 0)
    if not chunk.energy:
        chunk.energy = SpectralWCS(naxis, specsys)


    # determine the index of energy axis
    index = None
    for i in range(1, naxis + 1):
        if header['CTYPE{}'.format(i)] is not None:
           if header['CTYPE{}'.format(i)].split('-')[0].strip() in ENERGY_KEYWORDS:
               index = i
               break

    if not index:
        return

    chunk.ctype = header['CTYPE{}'.format(index).split('-')[0].strip()]
    #chunk.naxis = header['NAXIS{}'.format(index)]


    # Artifact.productType = artifact.productType
    # Artifact.releaseType = artifact.releaseType
    #
    # Part.name = part.name
    # Part.productType = part.productType
    #
    # Chunk.naxis = ZNAXIS, NAXIS
    # Chunk.observableAxis = chunk.observableAxis
    # Chunk.positionAxis1 = getPositionAxis()
    # Chunk.positionAxis2 = getPositionAxis()
    # Chunk.energyAxis = getEnergyAxis()
    # Chunk.timeAxis = getTimeAxis()
    # Chunk.polarizationAxis = getPolarizationAxis()
    #
    # Chunk.observable.dependent.bin = observable.dependent.bin
    # Chunk.observable.dependent.axis.ctype = observable.dependent.ctype
    # Chunk.observable.dependent.axis.cunit = observable.dependent.cunit
    # Chunk.observable.independent.bin = observable.independent.bin
    # Chunk.observable.independent.axis.ctype = observable.independent.ctype
    # Chunk.observable.independent.axis.cunit = observable.independent.cunit
    # Chunk.energy.specsys = SPECSYS
    # Chunk.energy.ssysobs = SSYSOBS
    # Chunk.energy.restfrq = RESTFRQ
    # Chunk.energy.restwav = RESTWAV
    # Chunk.energy.velosys = VELOSYS
    # Chunk.energy.zsource = ZSOURCE
    # Chunk.energy.ssyssrc = SSYSSRC
    # Chunk.energy.velang = VELANG
    # Chunk.energy.bandpassName = bandpassName
    # Chunk.energy.resolvingPower = resolvingPower
    # Chunk.energy.transition.species = energy.transition.species
    # Chunk.energy.transition.transition = energy.transition.transition
    # Chunk.energy.axis.axis.ctype = CTYPE{energyAxis}
    # Chunk.energy.axis.axis.cunit = CUNIT{energyAxis}
    # Chunk.energy.axis.bounds.samples = energy.samples
    # Chunk.energy.axis.error.syser = CSYER{energyAxis}
    # Chunk.energy.axis.error.rnder = CRDER{energyAxis}
    # Chunk.energy.axis.function.naxis = NAXIS{energyAxis}
    # Chunk.energy.axis.function.delta = CDELT{energyAxis}
    # Chunk.energy.axis.function.refCoord.pix = CRPIX{energyAxis}
    # Chunk.energy.axis.function.refCoord.val = CRVAL{energyAxis}
    # Chunk.energy.axis.range.start.pix = energy.range.start.pix
    # Chunk.energy.axis.range.start.val = energy.range.start.val
    # Chunk.energy.axis.range.end.pix = energy.range.end.pix
    # Chunk.energy.axis.range.end.val = energy.range.end.val