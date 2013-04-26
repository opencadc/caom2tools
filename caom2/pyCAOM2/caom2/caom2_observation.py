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

"""definition of the  caom2.Observation object."""

from caom2_entity import AbstractCaom2Entity
from caom2_algorithm import Algorithm
from caom2_environment import Environment
from caom2_instrument import Instrument
from caom2_plane import Plane
from caom2_proposal import Proposal
from caom2_target import Target
from caom2_telescope import Telescope
from caom2_observation_uri import ObservationURI
from caom2_enums import ObservationIntentType
from util.caom2_util import TypedOrderedDict
from util.caom2_util import validate_path_component
from util import caom2_util as util
from datetime import datetime


class Observation(AbstractCaom2Entity):
    """
    Observation

    An observation describes a set of empirical data.

    The CAOM2 observation is described at:
    http://www.cadc.hia.nrc.gc.ca/caom2

    The Observation object is the top level container for the
    meta-data assocaited with an observation.  The variaous attributes
    described here provide ways of storing the complete set of
    associated meta-data as well as links to the actual observation
    data.

    Observation -> Target
                -> Plane(s) 
                -> Instrument
                -> Telescope
                -> Proposal
                -> Environment

    Plane -> Artifact(s)

    Artifact -> Part(s)

    Part -> Chunk(s)

    The actual 'chunks' of observational data are reference by Chunk
    objects.  Information about the Spatial/Frequency/Time aspects of
    an Observation are expressed at the Chunk level.                 

    The Chunk contains refernces to caom2 objects that fully describe
    the circumstances of that chunk of observation.  Often a 'Chunk'
    is a single extension in a FITS image. But can also be a
    particular column of values from a FITS Table or some other data
    object.

    Chunk -> SpatialWCS
          -> TemporalWCS
          -> SpectralWCS
          -> PolarizationWCS
          -> (Observable) 
    """

    def __init__(self, 
                 collection, 
                 observation_id, 
                 algorithm,
                 sequence_number=None,
                 intent=None,
                 obs_type=None,
                 proposal=None,
                 telescope=None,
                 instrument=None,
                 target=None,
                 meta_release=None,
                 planes=None,
                 environment=None
                 ):
        """
        Initializes an Observation instance

        Arguments: collection : where the observation is from
        (eg. 'HST')

        observation_id : a unique identifier within that collection
                         (eg. '111')

        algorithm : the algorithm used to create the observation.  For
                    a telescope observation this is always 'exposure'
        """
        super(Observation, self).__init__()
        
        self.collection = collection
        self.observation_id = observation_id
        self.algorithm = algorithm

        self.uri = ObservationURI.get_observation_uri(collection,
                                                      observation_id)

        self.sequence_number = sequence_number
        self.intent = intent
        self.obs_type = obs_type
        self.proposal = proposal
        self.telescope = telescope
        self.instrument = instrument
        self.target = target
        self.environment = environment
        self.meta_release = meta_release
        if planes is None:
            planes = TypedOrderedDict((Plane),)
        self.planes = planes

        ## hold the list of attributes we should pprint
        self._print_attributes = ['collection','observation_id',
                                  'algorithm','sequence_number', 
                                  'intent', 'obs_type', 'proposal', 
                                  'telescope','instrument', 'target', 
                                  'environment', 'meta_release',
                                  'planes']


    # Properties
    @property
    def collection(self):
        """The name of the collection of observations, normally a telescope
        name.
        
        type: str
        """
        return self._collection

    @collection.setter
    def collection(self, value):
        util.typeCheck(value, str, 'collection', override=False)
        self._collection = value

    @property
    def observation_id(self):
        """A string that uniquely identifies this obseravtion within the given
        collection.

        type: str
        """
        return self._observation_id

    @observation_id.setter
    def observation_id(self, value):
        util.typeCheck(value, str, 'observation_id', override=False)
        self._observation_id = value

    @property
    def uri(self):
        """A URI for this observation referenced in the caom system.

        This attribute is auto geneqrated from the other metadata.
        type: str
        """
        return self._uri

    @uri.setter
    def uri(self, value):
        util.typeCheck(value, ObservationURI, 'uri')
        self._uri = value

    @property
    def planes(self):
        """A typed ordered dictionary containing plane objects associated with
        this observations.

        see caom2.Plane for details about a creating a plane.

        type: TypedOrderDict((Plane),)

        eg. Observation.planes.add(Plane("SCIENCE"))
        """
        return self._planes

    @planes.setter
    def planes(self, value):
        util.typeCheck(value, TypedOrderedDict,'planes')
        self._planes = value

    @property
    def algorithm(self):
        """The process that was used to select this observation.

        normally 'exposure' for an observation or the name of a group
        process for a composite observation"""
        return self._algorithm

    @algorithm.setter
    def algorithm(self, value):
        if isinstance(value, str):
            value = Algorithm(value)
        util.typeCheck(value, Algorithm, 'algorithm')
        self._algorithm = value

    @property
    def intent(self):
        """The original intent of having this data.
        
        type: ObservationIntentType
        
        see ObservationIntentType.names() for allowed values

        """
        return self._intent

    @intent.setter
    def intent(self, value):
        if isinstance(value, str):
            value = ObservationIntentType(value)
        util.typeCheck(value, ObservationIntentType, 'intent')
        self._intent = value

    @property
    def sequence_number(self):
        """An integer counter that reprsents this observation within a
        collection.  eg.  EXPNUM type: int
        """
        return self._sequence_number

    @sequence_number.setter
    def sequence_number(self, value):
        util.typeCheck(value, int, 'sequence_number')
        self._sequence_number = value

    @property
    def obs_type(self):
        """The OBSTYPE of the observation being recorded.

        eg. OBJECT, FLAT, BIAS
        type: str
        """
        return self._obs_type

    @obs_type.setter
    def obs_type(self, value):
        util.typeCheck(value, str, 'obs_type')
        self._obs_type = value

    @property
    def proposal(self):
        """Refence to a Proposal object that describe the science proposal
        that lead to this observation.

        can be None
        see caom2.Proposal for help on building a Proposal object
        type: caom2.Proposal
        """
        return self._proposal

    @proposal.setter
    def proposal(self, value):
        util.typeCheck(value, Proposal, "proposal")
        self._proposal = value

    @property
    def telescope(self):
        """Reference to a Telescope object associated with this observation.

        can be None
        type: caom2.Telescope
        """
        return self._telescope

    @telescope.setter
    def telescope(self, value):
        util.typeCheck(value, Telescope, 'telescope')
        self._telescope = value

    @property
    def instrument(self):
        """Reference to an Instrument object associated with this observation.

        can be None
        type: caom2.Instrument
        """
        return self._instrument

    @instrument.setter
    def instrument(self, value):
        if isinstance(value, str):
            value = Instrument(str)
        util.typeCheck(value, Instrument, "instrument")
        self._instrument = value

    @property
    def target(self):
        """Reference to a Target object associted with this observation.

        can be None
        type: caom2.Target
        """
        return self._target

    @target.setter
    def target(self, value):
        if isinstance(value, str):
            value = Target(str)
        util.typeCheck(value, Target, 'target')
        self._target = value

    @property
    def environment(self):
        """Reference to an Environment object associated with this
        observation.
        
        can be None
        type: caom2.Environment
        """
        return self._environment

    @environment.setter
    def environment(self, value):
        util.typeCheck(value, Environment,'environment')
        self._environment = value

    @property
    def meta_release(self):
        """A datetime value indicating when the meta-data of this observation
        is publicly accessible.

        This only controls access to the information about the
        observation. Access to the observational data is controlled
        via the Plane.data_release attribute (see the planes
        attribute).

        eg. '2012/11/28 12:00:00'
        type: datatime
        """
        return self._meta_release

    @meta_release.setter
    def meta_release(self, value):
        util.typeCheck(value, datetime, 'meta_release')
        self._meta_release = value
