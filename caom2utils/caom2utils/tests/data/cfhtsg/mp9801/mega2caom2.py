from re import findall, search


def uri2observationID(uri):
    id = search(r"(MegaPipe.\d\d\d.\d\d\d)", uri).group(0)
    return id


def uri2planeproductID(uri):
    id = search(r"(MegaPipe.\d\d\d.\d\d\d.[A-Z]+\.MP\d{4})", uri).group(0)
    return id


def getInputs(hdul):
    history = str(hdul[0].get('HISTORY'))

    # find all the exposure numbers:
    exposures = findall(r'\d{6,7}', history)

    # add the caom:CFHT/ bits:
    for i, value in enumerate(exposures):
        exposures[i] = 'caom:CFHT/' + exposures[i]

        # join them:
    inputs = ' '.join(exposures)

    # complain about the syntax for join:
    # just saying this has got the dumbest syntax for join in any
    # language I've seen
    return inputs


def getProvInputs(hdul):
    history = str(hdul[0].get('HISTORY'))

    # find all the exposure numbers:
    exposures = findall(r'\d{6,7}', history)

    # add the caom:CFHT/ bits:
    for i, value in enumerate(exposures):
        exposures[i] = 'caom:CFHT/' + exposures[i] + '/' + exposures[i] + 'p'

        # join them:
    inputs = ' '.join(exposures)

    # complain about the syntax for join:
    # just saying this has got the dumbest syntax for join in any
    # language I've seen
    return inputs


def getMagnitudeLimit(hdul):
    # get the magnitude limit keywords
    ml_5siga = float(hdul[0].get('ML_5SIGA'))
    ml_5sig2 = float(hdul[0].get('ML_5SIG2'))
    ml_5odet = float(hdul[0].get('ML_50DET'))
    # if (type(image) is str) :
    #     hdul.close
    # return the first one that doesn't suck
    if (ml_5siga > 0.000 and ml_5siga < 30.000):
        return ml_5siga
    if (ml_5sig2 > 0.000 and ml_5sig2 < 30.000):
        return ml_5sig2
    if (ml_5odet > 0.000 and ml_5odet < 30.000):
        return ml_5odet

    # return 0 as a default
    return 0.000


def getFilterWidth(uri):
    # given a string with a CFHT MegaCam filter name somewhere in it somewhere,
    # return the width of the filter in Angstroms
    filter = search(r"([A-Z]+\.MP\d{4})", uri).group(0)

    widthdict = {'U.MP9301': 867.0,
                 'G.MP9401': 1624.0,
                 'R.MP9601': 1413.0,
                 'I.MP9701': 1726.0,
                 'I.MP9702': 1728.0,
                 'Z.MP9801': 1848.0,

                 'U.MP9302': 934.0,
                 'G.MP9402': 1570.0,
                 'R.MP9602': 1524.0,
                 'I.MP9703': 1622.0,
                 'Z.MP9901': 1558.0,
                 'GRI.MP9605': 4241.0,

                 'HA.MP9603': 122.0,
                 'HAOFF.MP9604': 122.0,
                 'OIII.MP9501': 118.0,
                 'OIIIOFF.MP9502': 118.0,
                 'CAHK.MP9303': 110.0
                 }
    return widthdict[filter]


def getFilterCentre(uri):
    # given a string with a CFHT MegaCam filter name somewhere in it somewhere,
    # return the centre of the filter in Angstroms
    filter = search(r"([A-Z]+\.MP\d{4})", uri).group(0)
    centredict = {'U.MP9301': 3754.5,
                  'G.MP9401': 4890.0,
                  'R.MP9601': 6248.5,
                  'I.MP9701': 7763.0,
                  'I.MP9702': 7623.0,
                  'Z.MP9801': 9083.0,

                  'U.MP9302': 3523.0,
                  'G.MP9402': 4754.0,
                  'R.MP9602': 6414.0,
                  'I.MP9703': 7767.0,
                  'Z.MP9901': 9228.0,
                  'GRI.MP9605': 6100.0,

                  'HA.MP9603': 6591.0,
                  'HAOFF.MP9604': 6721.0,
                  'OIII.MP9501': 5007.0,
                  'OIIIOFF.MP9502': 5107.0,
                  'CAHK.MP9303': 3953.0
                  }
    return centredict[filter]
