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

"""defines the ObservableAxis class

"""

from caom2_slice import Slice
from caom2.util import caom2_util as util
from caom2.caom2_object import Caom2Object

class ObservableAxis(Caom2Object):
    """The slice through the data structure that provides the thing being
    described by this Axis.

    this data structure is used when a file contains data that has set
    of measured values (the dependent variable) and might also have
    the coordinate at which those values are measured in another Slice
    of the file.

    The Slice refers to a column in a FITS image. The bin is the
    column index and ctype/cunit for the axis describe what the column
    contiains.

    eg.

    NAXIS=2
    NAXIS1=3
    NAXIS2=N
    l1 f1 s1
    l2 f2 s2
    l3 f3 s3
    .
    .
    lN fN sN

    where l? is the wavelength at which a measure of flux f? has been
    made and l is the first column of a FITS data structure that is
    3,N in size.  s? is a third slice that would be used to define
    another observable. When defining the s? obserable the independent
    variable must be defined for that ObservableAxis too.

    The l/f obserable would be recorded as

    dependent=Slice(Axis('wave','nm'),bin=1)
    independent=Slice(Axis('flux','Jy'),bin=2)
    Chunk.observable_axis=ObservableAxis(dependent, independent)
    

    """
    def __init__(self, dependent, independent=None):


        self.dependent = dependent
        self.independent = independent

    @property
    def dependent(self):
        """The dependent (y) variable slice. 

        A slice provides the bin and the type/unit of the observable axis
        """
        
        return self._dependent

    @dependent.setter
    def dependent(self,value):
        util.typeCheck(value, Slice, 'dependent', override=False)
        self._dependent=value

    @property
    def independent(self):
        """The dependent (y) variable slice. 

        A slice that provides the pixel value and the type/unit for
        conversion

        """
        return self._independent

    @independent.setter
    def independent(self, value):
        util.typeCheck(value, Slice, "independent")
        self._independent = value
