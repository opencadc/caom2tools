caom2
=====

.. image:: https://img.shields.io/pypi/v/caom2.svg   
    :target: https://pypi.python.org/pypi/caom2

Common Archive Observation Model - data engineering tools

caom2 module

The caom2 module is a library implementing the Common Archive
Observation Model (CAOM-2.3) for manipulating CAOM observations and
reading and writing XML documents.

http://www.opencadc.org/caom2/

To create a minimal Simple Observation
--------------------------------------

.. code:: python

        import sys
        from caom2 import SimpleObservation, TypedOrderedDict, Plane, Artifact,\
                          Part, Chunk, ObservationWriter, ProductType,\
                          ReleaseType, TypedList

        observation = SimpleObservation('collection', 'observationID')

        observation.planes = TypedOrderedDict(Plane)
        plane = Plane('productID')
        observation.planes['productID'] = plane

        plane.artifacts = TypedOrderedDict(Artifact)
        artifact = Artifact('uri:foo/bar', ProductType.SCIENCE, ReleaseType.META)
        plane.artifacts['uri:foo/bar'] = artifact

        artifact.parts = TypedOrderedDict(Part)
        part = Part('name')
        artifact.parts['name'] = part

        part.chunks = TypedList(Chunk)
        chunk = Chunk()
        part.chunks.append(chunk)

        writer = ObservationWriter()
        writer.write(observation, sys.stdout)

The output:

.. code:: xml

    <?xml version='1.0' encoding='UTF-8'?>
    <caom2:Observation xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:type="caom2:SimpleObservation" caom2:id="00000000-0000-0000-3d6a-420eab45bf2e" caom2:lastModified="2016-11-24T08:40:54.003">
      <caom2:collection>collection</caom2:collection>
      <caom2:observationID>observationID</caom2:observationID>
      <caom2:algorithm>
        <caom2:name>exposure</caom2:name>
      </caom2:algorithm>
      <caom2:planes>
        <caom2:plane caom2:id="00000000-0000-0000-1fe6-420eab45c0bc" caom2:lastModified="2016-11-24T08:40:54.003">
          <caom2:productID>productID</caom2:productID>
          <caom2:artifacts>
            <caom2:artifact caom2:id="00000000-0000-0000-7adf-420eab45c170" caom2:lastModified="2016-11-24T08:40:54.004">
              <caom2:uri>uri:foo/bar</caom2:uri>
              <caom2:productType>science</caom2:productType>
              <caom2:releaseType>meta</caom2:releaseType>
              <caom2:parts>
                <caom2:part caom2:id="00000000-0000-0000-16f0-420eab45c246" caom2:lastModified="2016-11-24T08:40:54.004">
                  <caom2:name>name</caom2:name>
                  <caom2:chunks>
                    <caom2:chunk caom2:id="00000000-0000-0000-0c8f-420eab45c2a1" caom2:lastModified="2016-11-24T08:40:54.004"/>
                  </caom2:chunks>
                </caom2:part>
              </caom2:parts>
            </caom2:artifact>
          </caom2:artifacts>
        </caom2:plane>
      </caom2:planes>
    </caom2:Observation>

To create a complete Observation
--------------------------------------

.. code:: python

    # make it compatible with Python 2 and 3
    from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
    from datetime import datetime
    import sys
    from caom2 import SimpleObservation, Plane, Artifact, Part, Chunk,\
                      TypedOrderedDict, ObservationWriter, ProductType, \
                      ReleaseType, TypedList, Target, TargetPosition, \
                      TargetType, ObservationIntentType, Instrument, \
                      Telescope, Environment, DataProductType, Provenance, \
                      CalibrationLevel, Metrics, Proposal, Point, Slice, Axis,\
                      ObservableAxis, CoordAxis1D, CoordAxis2D, SpatialWCS,\
                      SpectralWCS, EnergyTransition, TemporalWCS, CoordFunction1D,\
                      RefCoord, PolarizationWCS

    observation = SimpleObservation('collection', 'observationID')
    observation.obs_type = 'flat'
    observation.intent = ObservationIntentType.SCIENCE
    observation.meta_release = datetime(2016, 11, 22, 11, 53, 44, 0)

    observation.proposal = Proposal('proposal id')
    observation.proposal.pi_name = 'pi name'
    observation.proposal.project = 'proposal project'
    observation.proposal.title = 'proposal title'
    observation.proposal.keywords.update({'proposal', 'key', 'words'})

    observation.target = Target('target name')
    observation.target.target_type = TargetType.OBJECT
    observation.target.standard = False
    observation.target.redshift = 1.5
    observation.target.keywords.update({'target', 'key', 'words'})

    point = Point(1.0, 2.0)
    observation.target_position = TargetPosition(point, 'coordsys')
    observation.target_position.equinox = 3.0

    observation.telescope = Telescope('telescope name')
    observation.telescope.geo_location_x = 1.0
    observation.telescope.geo_location_y = 2.0
    observation.telescope.geo_location_z = 3.0
    observation.telescope.keywords.update({'telescope', 'key', 'words'})

    observation.instrument = Instrument('instrument name')
    observation.instrument.keywords.update({'instrument', 'key', 'words'})

    observation.env = Environment()
    observation.env.seeing = 0.08
    observation.env.humidity = 0.35
    observation.env.elevation = 2.7
    observation.env.tau = 0.7
    observation.env.wavelength_tau = 450e-6
    observation.env.ambient_temp = 20.0
    observation.env.photometric = True

    observation.planes = TypedOrderedDict(Plane)
    plane = Plane('productID')
    observation.planes['productID'] = plane

    plane.meta_release = datetime(2016, 11, 22, 12, 26, 21, 0)
    plane.data_release = datetime(2018, 01, 01, 00, 00, 00, 0)
    plane.data_product_type = DataProductType.IMAGE
    plane.calibration_level = CalibrationLevel.PRODUCT

    plane.provenance = provenance = Provenance('name')
    plane.provenance.version = 'version'
    plane.provenance.product = 'product'
    plane.provenance.producer = 'producer'
    plane.provenance.run_id = 'run_id'
    plane.provenance.reference = 'http://foo/bar'
    plane.provenance.last_executed = datetime(2016, 11, 22, 12, 28, 16, 0)
    plane.provenance.keywords.update({'provenance', 'key', 'words'})

    plane.metrics = Metrics()
    plane.metrics.source_number_density = 1.0
    plane.metrics.background = 2.0
    plane.metrics.background_std_dev = 3.0
    plane.metrics.flux_density_limit = 4.0
    plane.metrics.mag_limit = 5.0

    plane.artifacts = TypedOrderedDict(Artifact)
    artifact = Artifact('uri:foo/bar', ProductType.SCIENCE, ReleaseType.META)
    plane.artifacts['uri:foo/bar'] = artifact

    artifact.content_type = 'application/fits'
    artifact.content_length = 12345L

    artifact.parts = TypedOrderedDict(Part)
    part = Part('name')
    artifact.parts['name'] = part
    part.product_type = ProductType.SCIENCE

    part.chunks = TypedList(Chunk)
    chunk = Chunk()
    part.chunks.append(chunk)

    chunk.product_type = ProductType.SCIENCE
    chunk.naxis = 5
    chunk.observable_axis = 1
    chunk.position_axis_1 = 1
    chunk.position_axis_2 = 2
    chunk.energy_axis = 3
    chunk.time_axis = 4
    chunk.polarization_axis = 5

    observable_axis = Slice(Axis('observable_ctype', 'observable_cunit'), 1L)
    chunk.observable = ObservableAxis(observable_axis)

    position_axis = CoordAxis2D(Axis('position_ctype_1', 'position_cunit_1'),
                                Axis('position_ctype_2', 'position_cunit_2'))
    chunk.position = SpatialWCS(position_axis)
    chunk.position.coordsys = 'position coordsys'
    chunk.position.equinox = 2000.0
    chunk.position.resolution = 0.5

    energy_axis = CoordAxis1D(Axis('energy_ctype', 'energy_cunit'))
    chunk.energy = SpectralWCS(energy_axis, 'specsys')
    chunk.energy.ssysobs = 'ssysobs'
    chunk.energy.ssyssrc = 'ssyssrc'
    chunk.energy.restfrq = 1.0
    chunk.energy.restwav = 2.0
    chunk.energy.velosys = 3.0
    chunk.energy.zsource = 4.0
    chunk.energy.velang = 5.0
    chunk.energy.bandpassName = 'bandpass name'
    chunk.energy.resolvingPower = 6.0
    chunk.energy.transition = EnergyTransition('H', '21cm')

    time_axis = CoordAxis1D(Axis('time_ctype', 'time_cunit'))
    chunk.time = TemporalWCS(time_axis)
    chunk.time.exposure = 1.0
    chunk.time.resolution = 2.0
    chunk.time.timesys = 'UTC'
    chunk.time.trefpos = 'TOPOCENTER'
    chunk.time.mjdref = 3.0

    polarization_axis = CoordAxis1D(Axis('STOKES'))
    polarization_axis.function = CoordFunction1D(4L, 1.0, RefCoord(1.0, 1.0))
    chunk.polarization = PolarizationWCS(polarization_axis)

    writer = ObservationWriter()
    writer.write(observation, sys.stdout)

The output:

.. code:: xml

	<?xml version='1.0' encoding='UTF-8'?>
	<caom2:Observation xmlns:caom2="http://www.opencadc.org/caom2/xml/v2.3" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:type="caom2:SimpleObservation" caom2:id="00000000-0000-0000-21ae-41feaaab49f6" caom2:lastModified="2016-11-23T13:35:24.404">
	  <caom2:collection>collection</caom2:collection>
	  <caom2:observationID>observationID</caom2:observationID>
	  <caom2:metaRelease>2016-11-22T11:53:44.000</caom2:metaRelease>
	  <caom2:algorithm>
		<caom2:name>exposure</caom2:name>
	  </caom2:algorithm>
	  <caom2:type>flat</caom2:type>
	  <caom2:intent>science</caom2:intent>
	  <caom2:proposal>
		<caom2:id>proposal id</caom2:id>
		<caom2:pi>pi name</caom2:pi>
		<caom2:project>proposal project</caom2:project>
		<caom2:title>proposal title</caom2:title>
		<caom2:keywords>proposal words key</caom2:keywords>
	  </caom2:proposal>
	  <caom2:target>
		<caom2:name>target name</caom2:name>
		<caom2:type>object</caom2:type>
		<caom2:standard>false</caom2:standard>
		<caom2:redshift>1.5</caom2:redshift>
		<caom2:keywords>words key target</caom2:keywords>
	  </caom2:target>
	  <caom2:targetPosition>
		<caom2:coordsys>coordsys</caom2:coordsys>
		<caom2:equinox>3.0</caom2:equinox>
		<caom2:coordinates>
		  <caom2:cval1>1.0</caom2:cval1>
		  <caom2:cval2>2.0</caom2:cval2>
		</caom2:coordinates>
	  </caom2:targetPosition>
	  <caom2:telescope>
		<caom2:name>telescope name</caom2:name>
		<caom2:geoLocationX>1.0</caom2:geoLocationX>
		<caom2:geoLocationY>2.0</caom2:geoLocationY>
		<caom2:geoLocationZ>3.0</caom2:geoLocationZ>
		<caom2:keywords>words key telescope</caom2:keywords>
	  </caom2:telescope>
	  <caom2:instrument>
		<caom2:name>instrument name</caom2:name>
		<caom2:keywords>instrument words key</caom2:keywords>
	  </caom2:instrument>
	  <caom2:planes>
		<caom2:plane caom2:id="00000000-0000-0000-f768-41feaaab4bbc" caom2:lastModified="2016-11-23T13:35:24.404">
		  <caom2:productID>productID</caom2:productID>
		  <caom2:metaRelease>2016-11-22T12:26:21.000</caom2:metaRelease>
		  <caom2:dataRelease>2018-01-01T00:00:00.000</caom2:dataRelease>
		  <caom2:dataProductType>image</caom2:dataProductType>
		  <caom2:calibrationLevel>3</caom2:calibrationLevel>
		  <caom2:provenance>
			<caom2:name>name</caom2:name>
			<caom2:version>version</caom2:version>
			<caom2:producer>producer</caom2:producer>
			<caom2:runID>run_id</caom2:runID>
			<caom2:reference>http://foo/bar</caom2:reference>
			<caom2:lastExecuted>2016-11-22T12:28:16.000</caom2:lastExecuted>
			<caom2:keywords>provenance words key</caom2:keywords>
		  </caom2:provenance>
		  <caom2:metrics>
			<caom2:sourceNumberDensity>1.0</caom2:sourceNumberDensity>
			<caom2:background>2.0</caom2:background>
			<caom2:backgroundStddev>3.0</caom2:backgroundStddev>
			<caom2:fluxDensityLimit>4.0</caom2:fluxDensityLimit>
			<caom2:magLimit>5.0</caom2:magLimit>
		  </caom2:metrics>
		  <caom2:artifacts>
			<caom2:artifact caom2:id="00000000-0000-0000-d905-41feaaab4ca0" caom2:lastModified="2016-11-23T13:35:24.404">
			  <caom2:uri>uri:foo/bar</caom2:uri>
			  <caom2:productType>science</caom2:productType>
			  <caom2:releaseType>meta</caom2:releaseType>
			  <caom2:contentType>application/fits</caom2:contentType>
			  <caom2:contentLength>12345</caom2:contentLength>
			  <caom2:parts>
				<caom2:part caom2:id="00000000-0000-0000-909d-41feaaab4d2d" caom2:lastModified="2016-11-23T13:35:24.405">
				  <caom2:name>name</caom2:name>
				  <caom2:productType>science</caom2:productType>
				  <caom2:chunks>
					<caom2:chunk caom2:id="00000000-0000-0000-2ef1-41feaaab4d74" caom2:lastModified="2016-11-23T13:35:24.405">
					  <caom2:productType>science</caom2:productType>
					  <caom2:naxis>5</caom2:naxis>
					  <caom2:observableAxis>1</caom2:observableAxis>
					  <caom2:positionAxis1>1</caom2:positionAxis1>
					  <caom2:positionAxis2>2</caom2:positionAxis2>
					  <caom2:energyAxis>3</caom2:energyAxis>
					  <caom2:timeAxis>4</caom2:timeAxis>
					  <caom2:polarizationAxis>5</caom2:polarizationAxis>
					  <caom2:observable>
						<caom2:dependent>
						  <caom2:axis>
							<caom2:ctype>observable_ctype</caom2:ctype>
							<caom2:cunit>observable_cunit</caom2:cunit>
						  </caom2:axis>
						  <caom2:bin>1</caom2:bin>
						</caom2:dependent>
					  </caom2:observable>
					  <caom2:position>
						<caom2:axis>
						  <caom2:axis1>
							<caom2:ctype>position_ctype_1</caom2:ctype>
							<caom2:cunit>position_cunit_1</caom2:cunit>
						  </caom2:axis1>
						  <caom2:axis2>
							<caom2:ctype>position_ctype_2</caom2:ctype>
							<caom2:cunit>position_cunit_2</caom2:cunit>
						  </caom2:axis2>
						</caom2:axis>
						<caom2:coordsys>position coordsys</caom2:coordsys>
						<caom2:equinox>2000.0</caom2:equinox>
						<caom2:resolution>0.5</caom2:resolution>
					  </caom2:position>
					  <caom2:energy>
						<caom2:axis>
						  <caom2:axis>
							<caom2:ctype>energy_ctype</caom2:ctype>
							<caom2:cunit>energy_cunit</caom2:cunit>
						  </caom2:axis>
						</caom2:axis>
						<caom2:specsys>specsys</caom2:specsys>
						<caom2:ssysobs>ssysobs</caom2:ssysobs>
						<caom2:ssyssrc>ssyssrc</caom2:ssyssrc>
						<caom2:restfrq>1.0</caom2:restfrq>
						<caom2:restwav>2.0</caom2:restwav>
						<caom2:velosys>3.0</caom2:velosys>
						<caom2:zsource>4.0</caom2:zsource>
						<caom2:velang>5.0</caom2:velang>
						<caom2:transition>
						  <caom2:species>H</caom2:species>
						  <caom2:transition>21cm</caom2:transition>
						</caom2:transition>
					  </caom2:energy>
					  <caom2:time>
						<caom2:axis>
						  <caom2:axis>
							<caom2:ctype>time_ctype</caom2:ctype>
							<caom2:cunit>time_cunit</caom2:cunit>
						  </caom2:axis>
						</caom2:axis>
						<caom2:timesys>UTC</caom2:timesys>
						<caom2:trefpos>TOPOCENTER</caom2:trefpos>
						<caom2:mjdref>3.0</caom2:mjdref>
						<caom2:exposure>1.0</caom2:exposure>
						<caom2:resolution>2.0</caom2:resolution>
					  </caom2:time>
					  <caom2:polarization>
						<caom2:axis>
						  <caom2:axis>
							<caom2:ctype>STOKES</caom2:ctype>
						  </caom2:axis>
						  <caom2:function>
							<caom2:naxis>4</caom2:naxis>
							<caom2:delta>1.0</caom2:delta>
							<caom2:refCoord>
							  <caom2:pix>1.0</caom2:pix>
							  <caom2:val>1.0</caom2:val>
							</caom2:refCoord>
						  </caom2:function>
						</caom2:axis>
					  </caom2:polarization>
					</caom2:chunk>
				  </caom2:chunks>
				</caom2:part>
			  </caom2:parts>
			</caom2:artifact>
		  </caom2:artifacts>
		</caom2:plane>
	  </caom2:planes>
	</caom2:Observation>
