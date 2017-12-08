caom2utils
=========

.. image:: https://img.shields.io/pypi/v/caom2utils.svg   
    :target: https://pypi.python.org/pypi/caom2utils

Utilities to facililate working with the CAOM2 model.


# Design Decisions

Is it ok to set 'lazy_load_hdus' to False in the astropy fits.open call? What's the maximum number of HDUs in the files CADC archives?


```
cvodb=>  select a.uri,count(*) as numparts from caom2.Artifact a join caom2.Part p on a.artifactID=p.artifactID where a.contentType='application/fits' group by a.artifactID having count(*) > 10 order by numparts desc limit 10;
                     uri                     | numparts
---------------------------------------------+----------
 ad:JCMT/jcmth20100515_00048_02_rimg_pro_000 |     2344
 ad:JCMT/jcmth20100515_00048_02_rsp_pro_000  |      579
 ad:HSTCA/oc6n01010_x1d                      |      334
 ad:HSTCA/oc6nf1010_x1d                      |      301
 ad:HSTCA/n8r7c2010_jit                      |      259
 ad:HSTCA/n8r7c2010_jif                      |      259
 ad:HSTCA/n8r7d2010_jif                      |      259
 ad:HSTCA/n8r7d2010_jit                      |      259
 ad:HSTCA/n8r7d3010_jif                      |      254
 ad:HSTCA/n8r7d3010_jit                      |      254
(10 rows)
```

Ran a test with the first JCMT file and the first HSTCA file on mach37 (2GB RAM) with no OOM issue.
