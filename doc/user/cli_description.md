# How to describe and load a CAOM2 Observation using the Command Line

1. Ensure the pre-conditions as described [here](https://github.com/SharonGoliath/caom2tools/tree/s2505/doc#preconditions)

1. Use the file [test_obs.blueprint](https://github.com/opencadc-metadata-curation/collection2caom2/blob/master/test_obs.blueprint) as the initial version of the blueprint file. For more information on the concept of blueprints, and their use, see [here](https://github.com/SharonGoliath/caom2tools/blob/s2505/doc/user/script_description.md#observation-blueprints).

1. Run caom2gen. The value provided for the `--local` parameter must be a fully qualified path name.

    ```
    caom2gen --out TEST_OBS.XML --observation COLLECTION TEST_OBS --blueprint ./test_obs.blueprint 
    --local /fully/qualified/path/TEST_FILE.FITS --lineage test_file/ad:COLLECTION/TEST_FILE.FITS
    ```

1. There should be a file named `TEST_OBS.XML` in the working directory.

1. Run caom2-repo. There will be no output if the command succeeds.

    ```
    caom2-repo create --netrc ./.netrc --resource-id ivo://cadc.nrc.ca/sc2repo TEST_OBS.XML
    ```

1. Everything after this is making refinements to the mapping between file content and CAOM2 instance members. This means issuing `caom2-repo update` commands, instead of `caom2-repo create` commands, to make changes on the server. However, caom2-repo is particular about its ids, so after the first successful execution of caom2-repo create, do this:

    ```
    caom2-repo read --netrc ./.netrc --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION TEST_OBS > TEST_OBS_READ.XML
    ```

1. There should be a file named `TEST_OBS_READ.XML` on disk. It will be different than the `TEST_OBS.XML` used for `create`, because the service generates a parallel set of keys that must be honoured. For each of the observation, plane, artifact, part, and chunk elements there are `id`, `lastModified`, `maxLastModified`, `metaChecksum`, and `accMetaChecksum` values. In particular, the `id` values must be consistent when doing `caom2-repo update` calls, or a "This observation already exists" error will occur.

1. After you've generated this output file, use the following commands to iteratively make and view changes to the mapping between the `COLLECTION` data and the CAOM2 instance:

    ```
    caom2gen -o TEST_OBS.XML --in TEST_OBS_READ.XML  --blueprint ./test_obs.blueprint --local /fully/qualified/path/TEST_FILE.FITS 
    --lineage TEST_FILE/ad:COLLECTION/TEST_FILE.FITS
    caom2-repo update --netrc ./.netrc --resource-id ivo://cadc.nrc.ca/sc2repo TEST_OBS.XML
    ```

1. In your browser, go to http://sc2.canfar.net/search, enter `TEST_OBS` into the `Observation ID` search field, click search, then click the `TEST_OBS` link in the `Obs. ID` column of the `Results` tab. This will display the details of the CAOM2 instance for `TEST_OBS` in a new tab.

1. Modify the blueprint to change mappings between the `COLLECTION` data model and the CAOM2 data model. If more complicated metadata mappings are required, investigate the use of the `--module` and `--plugin` parameters to [caom2gen](https://github.com/opencadc/caom2tools/tree/master/caom2utils). There are additional `caom2gen` parameters described here as well.

1. Should entries ever need to be deleted from the CAOM2 repository, replace `COLLECTION` with the appropriate value, and replace `TEST_OBS` with the observation ID that is being deleted:

    ```
    caom2-repo delete --netrc ./.netrc --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION TEST_OBS
    ```
