from astropy import units
from astropy.coordinates import SkyCoord


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
        logging.error(f'{ra} {dec}')
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
