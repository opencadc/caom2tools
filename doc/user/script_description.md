# How to describe and load a CAOM2 Observation using Python scripts

Ensure the pre-conditions described [here](../README.md).

The method `caom2utils.fits2caom2.augment` uses the concept of a blueprint to capture the description of a CAOM2 Observation as a 
mapping of a Telescope Data Model (`COLLECTION`) to the CAOM2 data model. This describes how to extend that application to customize the mapping for a `COLLECTION`.

`augment` works by creating or augmenting a CAOM2 Observation record, which can then be loaded via the CADC service.

`augment` creates the Observation record using information contained in a FITS file. The python module `fits2caom2`, from the python package `caom2utils`,
examines the FITS file and uses a blueprint, embodied in an instance of the ObsBlueprint class, to define default values, override values, and mappings to augment the FITS header. The keywords and values in the augmented FITS header are then used to fill in corresponding CAOM2 entities and attributes.

There are two alternate ways to provide input file metadata to the caom2gen application:
* have the file located on disk, and use the --local parameter
* have the file located in a CADC archive. The artifact URI portion of the lineage parameter will be used to resolve the archive and file name.

## Observation Blueprints

The blueprint is one way to capture the mapping of the telescope data model (TDM) to the CAOM2 data model. The blueprint can identify:
* what information to obtain from the FITS header, 
* defaults in case the FITS header is incomplete,
* hard-coded value when the FITS header should be ignored, or doesn't have information, and
* python functions which will be loaded and executed at run-time to augment FITS keyword values.

The blueprint is a set of key-value pairs, where the values have three possible representations.

The three representations are: defaults, overrides, and FITS keyword mappings.

There is a sample blueprint in [this file](https://github.com/opencadc-metadata-curation/collection2caom2/blob/master/test_obs.blueprint).

The keys are the long-form names for the CAOM2 model elements and attributes. The complete set of valid keys can be found by executing the following:

    python
    > import caom2utils
    > pydoc caom2utils.fits2caom2.ObsBlueprint

### Changing What a Blueprint Looks Like, By Extension

A blueprint may be provided by one of two ways: as a file on disk, or programmatically.

#### File Blueprint Usage

    Observation.observationID = ['OBSID'], default = TEST_OBS
    Plane.dataRelease = 2017-08-31T00:00:00
    Chunk.position.coordsys = ['RADECSYS,RADESYS']

   * Observation.observationID provides a default value of `test_obs`, which is used if the keyword `OBSID` does not exist in the FITS file.
   * Plane.dataRelease provides an override value, which is always used.
   * Chunk.position.coordsys provides a list of FITS keywords to try. If the first value is not in the FITS header, the second one is queried. If neither of them exist, there were be no value for Chunk.position.coordsys in the CAOM2 observation.

#### Programmatic Blueprint Usage 

An example of this implementation is in (https://github.com/opencadc-metadata-curation/vlass2caom2)

    bp = ObsBlueprint(position_axes=(1,2), time_axis=3, energy_axis=4, polarization_axis=5, observable_axis=6)
    bp.set_default('Observation.observationID', 'TEST_OBS')
    bp.set('Plane.dataRelease', '2017-08-31T00:00:00')
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADESYS')

   * Observation.observationID provides a default value of `TEST_OBS`, which is used if the keyword `OBSID` does not exist in the FITS file.
   * Plane.dataRelease provides an override value, which is always used when setting the plane-level data release data in the CAOM2 instance.
   * Chunk.position.coordsys provides a list of FITS keywords to try.

To make WCS content available in the blueprint, as opposed to setting the indices in the ObsBlueprint constructor, as many of the following functions for which there is metadata in a FITS file, may be called on a blueprint instance:

    bp = ObsBlueprint()
    bp.configure_position_axes((1, 2)) 
    bp.configure_energy_axis(3) 
    bp.configure_time_axis(4)
    bp.configure_polarization_axis(5) 
    bp.configure_observable_axis(6)

## Putting It All Together

The following script is an end-to-end example of describing and loading a CAOM2 Observation given a FITS file, and a blueprint.

    ```
    from cadcutils import net
    from caom2 import obs_reader_writer, DataProductType, CalibrationLevel
    from caom2repo import CAOM2RepoClient
    from caom2utils import fits2caom2

    # configure and create the CADC service client
    netrc_fqn = '/usr/src/app/.netrc'
    subject = net.Subject(netrc=netrc_fqn)
    # remove the resource_id parameter to use production resources
    repo_client = CAOM2RepoClient(subject, resource_id='ivo://cadc.nrc.ca/sc2repo')

    # describe the Observation by setting up the mapping between the
    # COLLECTION and CAOM2
    bp = fits2caom2.ObsBlueprint()
    bp.configure_position_axes((1, 2))
    bp.set_default('Observation.observationID', 'TEST_OBS')
    bp.set('Plane.dataRelease', '2017-08-31T00:00:00')
    bp.set('Plane.calibrationLevel', CalibrationLevel.RAW_STANDARD)
    bp.set('Plane.dataProductType', DataProductType.IMAGE)
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')

    # apply the mapping to the FITS file
    kwargs = {}
    uri = 'ad:COLLECTION/TEST_FILE.FITS'
    blueprints = {uri: bp}
    result = fits2caom2.augment(blueprints=blueprints,
                                no_validate=False,
                                dump_config=False,
                                plugin=None,
                                out_obs_xml='./TEST_OBS.XML',
                                in_obs_xml=None,
                                collection='COLLECTION',
                                observation='TEST_OBS',
                                product_id='TEST_PRODUCT_ID',
                                uri=uri,
                                netrc=netrc_fqn,
                                file_name='file:///usr/src/app/TEST_FILE.FITS',
                                verbose=False,
                                debug=False,
                                quiet=False,
                                caom_namespace=obs_reader_writer.CAOM23_NAMESPACE,
                                **kwargs)

    # load the observation into memory
    reader = obs_reader_writer.ObservationReader(False)
    observation = reader.read('./TEST_OBS.XML')

    # create the observation record with the service
    #
    # use 'update' if the observation has already been loaded to the CADC service. The service
    # generates a parallel set of database keys that must be honoured. The id values must be
    # consistent when doing 'create' and 'update' calls, or a "This observation already exists"
    # error will occur.
    repo_client.create(observation)
    ```

## More Information

If you want:

* to add direct CAOM2 model manipulation to your script, see [here](https://github.com/opencadc/caom2tools/tree/master/caom2) for an introduction to the possibilities.

* a description of the latest version of the CAOM2 model, see [here](http://www.opencadc.org/caom2/).

* a description of the operational version of the CADC service, see [here](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/ams/).

* examples of the model and the client embedded in end-to-end workflows, see [here](https://github.com/opencadc-metadata-curation). Each application in this repository uses the tactic of programatically creating a unique blueprint for each file that is ingested, and then creating or updating the resulting CAOM2 Observation.
