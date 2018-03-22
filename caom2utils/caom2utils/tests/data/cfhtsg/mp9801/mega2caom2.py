from astropy.io import fits
from re import findall, search, compile

def image2ID(image) :
    id = ''
    if (type(image) is str) :
        id=image.replace('.fits','')
    return id

def getInputs (image) :
    ### open the image or just use the header provided
    if (type(image) is str) :
        hdul= fits.open(image)
    else :
        hdul=image

    history= str(hdul[0].get('HISTORY'))

    ### find all the exposure numbers:
    exposures=findall(r'\d{6,7}',history)

    ### add the caom:CFHT/ bits:
    for i,value in enumerate(exposures) :
        exposures[i]='caom:CFHT/'+exposures[i]
        
        ### join them:
    inputs=','.join(exposures)  
        
    ### complain about the syntax for join:
    ### just saying this has got the dumbest syntax for join in any language I've seen
    return inputs


def getMagnitudeLimit (image) :
### open the image or just use the header provided
    if (type(image) is str) :
        hdul= fits.open(image)
    else :
        hdul=image

### get the magnitude limit keywords 
    # ml_5siga = float(hdul[0].header['ML_5SIGA'])
    # ml_5sig2 = float(hdul[0].header['ML_5SIG2'])
    # ml_5odet = float(hdul[0].header['ML_50DET'])
    ml_5siga = float(hdul[0].get('ML_5SIGA'))
    ml_5sig2 = float(hdul[0].get('ML_5SIG2'))
    ml_5odet = float(hdul[0].get('ML_50DET'))
    if (type(image) is str) :
        hdul.close
### return the first one that doesn't suck
    if (ml_5siga >0.000 and ml_5siga<30.000) : return ml_5siga
    if (ml_5sig2 >0.000 and ml_5sig2<30.000) : return ml_5sig2
    if (ml_5odet >0.000 and ml_5odet<30.000) : return ml_5odet

### return 0 as a default
    return 0.000



def getFilterWidth(image) :
### given a string with a CFHT MegaCam filter name somewhere in it somewhere, 
### return the width of the filter in Angstroms

    filter=search(r"([A-Z]+\.MP\d{4})",image).group(0)

    widthdict = { 'U.MP9301'       :  867.0,
                  'G.MP9401'       : 1624.0,
                  'R.MP9601'       : 1413.0,
                  'I.MP9701'       : 1726.0,
                  'I.MP9702'       : 1728.0,
                  'Z.MP9801'       : 1848.0,
                  
                  'U.MP9302'       :  934.0,
                  'G.MP9402'       : 1570.0,
                  'R.MP9602'       : 1524.0,
                  'I.MP9703'       : 1622.0,
                  'Z.MP9901'       : 1558.0,
                  'GRI.MP9605'     : 4241.0,
                  
                  'HA.MP9603'      :  122.0,
                  'HAOFF.MP9604'   :  122.0,
                  'OIII.MP9501'    :  118.0,
                  'OIIIOFF.MP9502' :  118.0,
                  'CAHK.MP9303'    :  110.0
              }
    return widthdict[filter]
    

def getFilterCentre(image) :
### given a string with a CFHT MegaCam filter name somewhere in it somewhere, 
### return the centre of the filter in Angstroms
    filter=search(r"([A-Z]+\.MP\d{4})",image).group(0)
    centredict = { 'U.MP9301'       : 3754.5, 
                   'G.MP9401'       : 4890.0, 
                   'R.MP9601'       : 6248.5, 
                   'I.MP9701'       : 7763.0, 
                   'I.MP9702'       : 7623.0, 
                   'Z.MP9801'       : 9083.0, 
                   
                   'U.MP9302'       : 3523.0,
                   'G.MP9402'       : 4754.0,
                   'R.MP9602'       : 6414.0,
                   'I.MP9703'       : 7767.0,
                   'Z.MP9901'       : 9228.0,
                   'GRI.MP9605'     : 6100.0,
                   
                   'HA.MP9603'      : 6591.0,
                   'HAOFF.MP9604'   : 6721.0,
                   'OIII.MP9501'    : 5007.0,
                   'OIIIOFF.MP9502' : 5107.0,
                   'CAHK.MP9303'    : 3953.0
               }
    return centredict[filter]
