caom2utils
=========

.. image:: https://img.shields.io/pypi/v/caom2utils.svg   
    :target: https://pypi.python.org/pypi/caom2utils

Utilities to facililate working with the CAOM2 model.


Design Decisions
================

1. Is it ok to set 'lazy_load_hdus' to False in the astropy fits.open call? What's the maximum number of HDUs in the files CADC archives?::

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


Ran a test with the first JCMT file and the first HSTCA file on mach37 (2GB RAM) with no OOM issue.

Tests
=====

CFHT
----

* 2087482o  MegaPrime raw

* 2087482p  MegaPrime processed

* 1916216o Espadons raw

* 1916216i  Espadons processed  - weird format

* 1916216p Espadons polarization composite  - weird format

* 2216850f  Espadons flat  - calibration, not on-sky

* 2216860b Espadons bias   - calibration, not on-sky

* 1709071o  WIRCam raw

* 1709071g  WIRCam raw guide cube - little time-series cutouts around some stars::

   curl --location-trusted -g -o 1709071g.fits 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/CFHT/1709071g.fits?cutout=[0][1:1,1:1,1:1]&cutout=[1][1:1,1:1,1:1]&cutout=[2][1:1,1:1,1:1]&cutout=[3][1:1,1:1,1:1]&cutout=[4][1:1,1:1,1:1]&cutout=[5][1:1,1:1,1:1]'

   def test_update():
      hdulist = fits.open(sample_file_time_axes)
      header = hdulist[0].header

      for ii in range(6):
          print(ii)
          if ii == 5:
              break
          header = hdulist[ii+1].header
          header['CTYPE3'] = 'TIME'
          header['CUNIT3'] = 'd'
          header['CSYER3'] = 1e-07
          header['CRDER3'] = 1e-07
          header['CRPIX3'] = 0.5
          header['CRVAL3'] = 56789.4298069
          header['CDELT3'] = 2.31481e-07
          header['NAXIS3'] = 1
          header['EXPTIME'] = 0.02
          header['TIMEDEL'] = 0.02
          header['TIMESYS'] = 'UTC'

      hdulist.writeto(sample_file_time_axes, overwrite=True)


* 1709071p  WIRCam processed

* 2136164o SITELLE raw

* 2136164p SITELLE processed spectral cube

CGPS
----

* CGPS_MA1_HI_line_image.fits

    def test_update():
        hdulist = fits.open(sample_file_4axes_obs)
        header = hdulist[0].header
        header['RUNID'] = 'HI-line'
        header['PROCNAME'] = 'exposure'
        header['INSTRUME'] = 'DRAO-ST'
        header['OBSGEO-X'] = -2100330.87517
        header['OBSGEO-Y'] = -3694247.82445
        header['OBSGEO-Z'] = 4741018.33097
        header['DATE-FTS'] = '2000-10-16'
        header['OBSID'] = 'MA1_DRAO-ST'
        header['XREFER'] = 'http://dx.doi.org/10.1086/375301'
        header['XPRVNAME'] = 'CGPS MOSAIC'
        hdulist.writeto(sample_file_4axes_obs, overwrite=True)


