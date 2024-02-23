import logging
from astropy.time import Time


def get_exposure(uri):
    return 1200.0


def _get_energy_chunk_range_end_val(header):
    result = None
    camera = header.get('CAMERA')
    if camera == 'RED':
        result = 1060.0  # nm
    elif camera == 'BLUE':
        result = 530.0  # nm
    return result


def _get_energy_chunk_range_start_val(header):
    result = None
    camera = header.get('CAMERA')
    if camera == 'RED':
        result = 520.0  # nm
    elif camera == 'BLUE':
        result = 347.0  # nm
    return result


def get_time_delta(header):
    result = None
    date_obs = header.get('DATE-OBS')
    ut_end = header.get('UTEND')
    ut_start = header.get('UTSTART')
    if date_obs and ut_end and ut_start:
        temp_start = f'{date_obs} {ut_start}'
        temp_end = f'{date_obs} {ut_end}'
        start = Time(temp_start)
        end = Time(temp_end)
        if start and end:
            start.format = 'mjd'
            end.format = 'mjd'
            result = (end - start).value
        else:
            logging.debug(f'Cannot convert {temp_start} or {temp_end} to MJD for {header.get("EXPID")}')
    else:
        logging.error(
            f'Missing one of DATE-OBS {date_obs}, UTSTART {ut_start}, or UTEND {ut_end} in  {header.get("EXPID")}'
        )
    return result


def get_time_function_val(header):
    result = None
    date_obs = header.get('DATE-OBS')
    ut_start = header.get('UTSTART')
    if date_obs and ut_start:
        temp_start = f'{date_obs} {ut_start}'
        start = Time(temp_start)
        if start:
            start.format = 'mjd'
            result = start.value
        else:
            logging.debug(f'Cannot convert {temp_start} to MJD')
    else:
        logging.error(f'Missing one of DATE-OBS {date_obs} or UTSTART {ut_start} in  {header.get("EXPID")}')
    return result
