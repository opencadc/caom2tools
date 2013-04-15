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

"""defines the caom2.Plane class"""


from caom2_entity import AbstractCaom2Entity
from caom2_artifact import Artifact
from caom2_metrics import Metrics
from caom2_provenance import Provenance
from caom2_enums import CalibrationLevel
from caom2_enums import DataProductType
from util.caom2_util import TypedOrderedDict
from util.caom2_util import typeCheck
from util.caom2_util import valueCheck
from datetime import datetime
from caom2.util import caom2_util as util

class Plane(AbstractCaom2Entity):
    """ Plane class """

    def __init__(self, product_id,
                 artifacts=None,
                 meta_release=None,
                 data_release=None,
                 data_product_type=None,
                 calibration_level=None,
                 provenance=None,
                 metrics=None):
        """
        Initialize a Plane instance

        Arguments:
        product_id : product ID
        """
        super(Plane, self).__init__()
        self.product_id = product_id
        if artifacts is None:
            artifacts = TypedOrderedDict((Artifact),)
        self._artifacts = artifacts


        self.meta_release = None
        self.data_release = None
        self.data_product_type = None
        self.calibration_level = None
        self.provenance = None
        self.metrics = None

        # computed fields
        # agregated from the Chunks during ingestion
        self._position = None
        self._energy = None
        self._time = None
        self._polarization = None

        self._print_attributes = ['product_id', 'meta_release', 'data_release',
                                  'data_product_type', 'calibration_level', 'provenance',
                                  'metrics', 
                                  #'position', 'time', 'energy','polarization'
                                  'artifacts']

    def _key(self):
        return (self.product_id)

    def __hash__(self):
        return hash(self._key())

    # Properties
    @property
    def product_id(self):
        """A string that identifies the data product, within a given
        observation, that is stored in this plane.

        eg: '1234567p'
        type: str
        """
        return self._product_id

    @product_id.setter
    def product_id(self, value):
        util.typeCheck(value, str, 'product_id', override=False)
        self._product_id = value

    @property
    def key(self):
        """ The dictionary key for a plane is product ID """
        return self._product_id

    @property
    def artifacts(self):
        """A TypeList of artifacts that are part of this plane.
        
        individual artifacts are constructed and then added to the plane.

        eg. Plane.artifacts.add(Artifact('ad:CFHT/1234567p')), see the
        Arifact help for more info on making an Aritfact object return
        """
        return self._artifacts

    @artifacts.setter
    def artifacts(self,value):
        util.typeCheck(value,TypedOrderedDict,'artifacts', override=False)
        self._artifacts = value

    @property
    def meta_release(self):
        """The date when the metadata describing this observation become
        public.

        eg.  Plane.meta_releaes=datetime.datetime(2012,1,1,0,0,0)
        indicates that the metadata become public at 00:00:00 on
        January 1st 2012.

        if there is no meta_release period, set value to None

        unit: calendar date
        type: datetime
        """
        return self._meta_release

    @meta_release.setter
    def meta_release(self, value):
        util.typeCheck(value, datetime, 'meta_release')
        util.valueCheck(value, 
                        datetime(1800,1,1,0,0,0), 
                        datetime(2050,1,1,0,0,0), 
                        'meta_release')
        self._meta_release = value

    @property
    def data_release(self):
        """The date when the data contained in this plane become public.

        eg.  Plane.data_releaes=datetime.datetime(2012,1,1,0,0,0)
        indicates that ar 00:00:00 on January 1st 2012 the data
        assocaited with this Plane will become publicly accessible.

        If there is no data_release period, set value to None.

        unit: calendar date
        type: datetime
        """
        return self._data_release

    @data_release.setter
    def data_release(self, value):
        util.typeCheck(value, datetime, 'data_release')
        util.valueCheck(value, 
                        datetime(1800,1,1,0,0,0), 
                        datetime(2050,1,1,0,0,0), 
                        'data_release')
        self._data_release = value

    @property
    def data_product_type(self):
        """The type of file structure that this plane contains.

        eg.
        Plane.data_product_type = 'EVENTLIST'

        see DataProductType.names() for allowed values

        """
        return self._data_product_type

    @data_product_type.setter
    def data_product_type(self, value):
        util.typeCheck(value, DataProductType, 'data_product_type')
        self._data_product_type = value

    @property
    def calibration_level(self):
        """a string that represents the level of calibrattion (aka processing)
        the data contained in this plane have reecieved.  The string
        is converted to an integer during storage.

        eg. Plane.calibration_level = "RAW_STANDARD"
        type: str

        Must be one of CalibrationLevel.names()

	"""
        return self._calibration_level

    @calibration_level.setter
    def calibration_level(self, value):
        util.typeCheck(value, CalibrationLevel, "calibration_level")
        self._calibration_level = value

    @property
    def provenance(self):
        """The process that created the data referred to by this Plane.

        eg.  Plane.provenance=caom2.Provenance("Elixir")
        """
        return self._provenance

    @provenance.setter
    def provenance(self, value):
        util.typeCheck(value, Provenance, "provenance")
        self._provenance = value

    @property
    def metrics(self):
        """reference to an object that contains metrics of this plane.

        eg. Plane.metrics = caom2.Metrics()
        """
        return self._metrics

    @metrics.setter
    def metrics(self, value):
        util.typeCheck(value, Metrics, 'metrics')
        self._metrics = value

    #@property
    #def observable(self):
    #    """ """
    #    return self._observable

    @property
    def position(self):
        """A caom2 Position object that is developed from
        the agregation of the Chunks that are children of 
        the Plane. 
        
        agregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._position

    @property
    def energy(self):
        """A caom2 Energy object that is developed from
        the agregation of the Chunks that are children of 
        the Plane. 
        
        agregation happens during ingest and is not part
        of the python module at this time.
        """
        """ Energy """
        return self._energy

    @property
    def time(self):
        """A caom2 Time object that is developed from
        the agregation of the Chunks that are children of 
        the Plane. 
        
        agregation happens during ingest and is not part
        of the python module at this time.
        """
        """ Time """
        return self._time

    @property
    def polarization(self):
        """A caom2 Polarization object that is developed from
        the agregation of the Chunks that are children of 
        the Plane. 
        
        agregation happens during ingest and is not part
        of the python module at this time.
        """
        return self._polarization

    # Compute derived fields

    def compute_position(self):
        raise NotImplementedError("Agregation of position has not been implemenetd in this module")

    def compute_energy(self):
        raise NotImplementedError("Agregation of energy has not been implemenetd in this module")

    def compute_time(self):
        raise NotImplementedError("Agregation of time has not been implemenetd in this module")

    def compute_polarization(self):
        raise NotImplementedError("Agregation of polarization has not been implemenetd in this module")
