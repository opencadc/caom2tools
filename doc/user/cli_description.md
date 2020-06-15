# How to describe and load a CAOM2 Observation using the Command Line

1. Ensure the pre-conditions as described [here](https://github.com/SharonGoliath/caom2tools/tree/s2505/doc#preconditions)

1. Use the file [test_obs.blueprint](https://github.com/opencadc-metadata-curation/collection2caom2/blob/master/test_obs.blueprint) as the initial version of the blueprint file.

1. Run caom2gen.

    ```
    caom2gen --out test_obs.xml --observation COLLECTION test_obs --blueprint ./test_obs.blueprint 
    --local ./test_file.fits --lineage test_file/ad:COLLECTION/test_file.fits
    ```

1. There should be a file named `test_obs.xml` in the working directory.

1. Run caom2-repo. There will be no output if the command succeeds.

    ```
    caom2-repo create --netrc ./netrc --resource-id ivo://cadc.nrc.ca/sc2repo test_obs.xml
    ```

1. Everything after this is making refinements to the mapping between file content and CAOM2 instance members. This means issuing `caom2-repo update` commands, instead of `caom2-repo create` commands, to make changes on the server. However, caom2-repo is particular about its ids, so after the first successful execution of caom2-repo create, do this:

    ```
    caom2-repo read -n --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION test_fits > test_obs_read.xml
    ```

1. There should be a file named `test_obs_read.xml` on disk. It will be different from `test_obs.xml`, because the service generates a parallel set of keys that must be honoured. If you look in this file, you will see that it is different than the `test_obs.xml` used for the create. For each of the observation, plane, artifact, part, and chunk elements there are id, lastModified, maxLastModified, metaChecksum, and accMetaChecksum values. In particular, the id values must be consistent when doing `caom2-repo update` calls, or a "This observation already exists" error will occur.

1. After you've generated this output file, use the following commands to repeatedly make and view changes to the mapping between the `COLLECTION` data and the CAOM2 instance:

    ```
    caom2gen -o test_obs.xml --in test_obs_read.xml  --blueprint ./test_obs.blueprint --local ./test_file.fits 
    --lineage test_file/ad:COLLECTION/test_file.fits
    caom2-repo update --netrc ./netrc --resource-id ivo://cadc.nrc.ca/sc2repo test_obs.xml
    ```

1. In your browser, go to http://sc2.canfar.net/search, enter `test_fits` into the `Observation ID` search field, click search, then click the `test_fits` link in the `Obs. ID`. This will display the details of the CAOM2 instance for `test_Fits` in a new tab.

1. Modify the blueprint to change mappings between the `COLLECTION` data model and the CAOM2 data model. If more complicated metadata mappings are required, investigate the use of the --module and --plugin parameters to [caom2gen](https://github.com/opencadc/caom2tools/tree/master/caom2utils). There are addition `caom2gen` parameters described here as well.

1. Should entries ever need to be deleted from the CAOM2 repository, replace `COLLECTION` with the appropriate value, and replace `test_fits` with the observation ID that is being deleted:

    ```
    caom2-repo delete --netrc ./netrc --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION test_fits
    ```
