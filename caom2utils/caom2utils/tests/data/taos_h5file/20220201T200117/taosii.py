from astropy import units
from astropy.coordinates import SkyCoord


def get_exposure(base):
    b = base.get('base')
    mjdrunstart = b['header']['timeseries']['mjdrunstart']
    mjdrunend = b['header']['timeseries']['mjdrunend']
    result = 0.0
    if mjdrunstart is not None and mjdrunend is not None:
        result = mjdrunend - mjdrunstart
    return result


def get_target_position_cval1(base):
    ra, dec_ignore = _get_target_position(base)
    return ra


def get_target_position_cval2(base):
    ra_ignore, dec = _get_target_position(base)
    return dec


def _get_target_position(base):
    import logging
    b = base.get('base')
    try:
        ra = b['header']['object']['obj_ra']
        dec = b['header']['object']['obj_dec']
        # logging.error(f'{ra} {dec}')
        result = SkyCoord(
            ra.decode('utf-8'),
            dec.decode('utf-8'),
            frame='icrs',
            unit=(units.hourangle, units.deg),
        )
        return result.ra.degree, result.dec.degree
    except Exception as e:
        import logging
        import traceback
        logging.error(e)
        logging.error(traceback.format_exc())
        raise e


def get_time_axis_range_end(base):
    b = base.get('base')
    x = b['header']['timeseries']['numepochs']
    return x - 1
