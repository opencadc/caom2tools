# How to describe and load a CAOM2 Observation using the Command Line

1. This description assumes:
    1. a working knowledge of python. [Prefer python3, please](https://pythonclock.org/),
    1. a linux-type environment,
    1. a working directory location, where all files discussed are placed, and
    1. that you have a [CADC account](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/auth/request.html), which is configured by CADC to have read and write access to a CAOM `COLLECTION`.

1. This description uses the parameters `test_file.fits`, `test_obs.xml` and `COLLECTION`. Replace these values appropriately when executing the commands.

1. Copy the file `test_file.fits` in the working directory. This metadata in this file will be described in the CAOM Observation created during this example.

1. The example will cause an instance to be created in the [CAOM2 sandbox](http://sc2.canfar.net/search/).  If you visit prior to the creation of the first CAOM instance for a collection, that collection will not show in the Collection data train. Even after successful creation of a CAOM instance for a collection, it can take up to one day for the collection to be displayed on the UI.

    To use production CADC services, remove `resource-id` parameters.

1. Install the following python dependencies:

    ```
    pip install caom2repo
    pip install caom2utils
    ```

1. Get credentials organized. The examples assume the use of a [.netrc file](https://www.systutorials.com/docs/linux/man/5-netrc/). The .netrc file content should include the following, with cadcusername and cadcpassword replaced with your CADC username and password values:

    ````
    machine www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca login canfarusername password canfarpassword
    machine www.canfar.net login canfarusername password canfarpassword
    machine ws-cadc.canfar.net login canfarusername password canfarpassword
    machine sc2.canfar.net login canfarusername password canfarpassword
    ````

1. The caom2-repo client also supports username/password and X509 certificates. If you want to use X509 
certificates use the --cert parameter instead of the -n parameter in all the commands. Obtaining and using a CADC certificate is described TBD.

1. Test the install. Commands are case-sensitive.

    ```
    caom2-repo read --netrc ./netrc --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION abc
    ```

    This will report an error:

    ```
    Client Error: Not Found for url: http://sc2.canfar.net/sc2repo/auth-observations/COLLECTION/abc.
    ```

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
    caom2-repo update -n --resource-id ivo://cadc.nrc.ca/sc2repo test_obs.xml
    ```

1. In your browser, go to http://sc2.canfar.net/search, enter `test_fits` into the `Observation ID` search field, click search, then click the `test_fits` link in the `Obs. ID`. This will display the details of the CAOM2 instance for `test_Fits` in a new tab.

1. Modify the blueprint to change mappings between the `COLLECTION` data model and the CAOM2 data model. If more complicated metadata mappings are required, investigate the use of the --module and --plugin parameters to [caom2gen](https://github.com/opencadc/caom2tools/tree/master/caom2utils).

1. [How To Use caom2gen](https://github.com/opencadc-metadata-curation/collection2caom2/wiki/How-To-Use-caom2gen).

1. Should entries ever need to be deleted from the CAOM2 repository, replace `COLLECTION` with the appropriate value, and replace `test_fits` with the observation ID that is being deleted:

    ```
    caom2-repo delete -n --resource-id ivo://cadc.nrc.ca/sc2repo COLLECTION test_fits
    ```
