from astropy.time import Time
from datetime import datetime


def _get_datetime(base):
    b = base.get('base').attrs
    result = None
    d = b.get('OBS_DATE')
    t = b.get('OBS_TIME')
    if d is not None and t is not None:
        dt = f'{d} {t}'
        result = Time(datetime.strptime(dt, '%Y-%m-%d %H:%M:%S.%f'))
        result.format = 'mjd'
        result = result.value
    return result


def _get_energy_resolving_power(base):
    b = base.get('base').attrs
    result = None
    # Laurie Rousseau-Nepton - 11-08-22
    # Resolving Power could be given at the central wavelength of the filter.
    # The formula is R = 1/lambda[nm]* (2*(STEP[nm]*(NAXIS3-zpd_index))/1.2067
    step = b.get('STEP')
    zpd_index = b.get('zpd_index')
    naxis_3 = b.get('step_nb')
    filter_max = b.get('filter_nm_max')
    filter_min = b.get('filter_nm_min')
    wl = None
    if filter_max is not None and filter_min is not None:
        wl = (filter_min + filter_max) / 2
    if (
        step is not None
        and zpd_index is not None
        and naxis_3 is not None
        and wl is not None
    ):
        result = 1 / wl * 2 * (step * (naxis_3 - zpd_index)) / 1.2067
    return result


def _get_exposure(base):
    b = base.get('base').attrs
    # Laurie Rousseau-Nepton - 11-08-22
    # Int. Time could be the total (multiplied by the cube spectral dimension
    # f.attrs.get(‘NAXIS3’)
    result = None
    exposure = b.get('exposure_time')
    naxis_3 = b.get('step_nb')
    if exposure is not None and naxis_3 is not None:
        result = exposure * naxis_3
    return result


def _get_fwhm(base):
    b = base.get('base').attrs
    minimum = b.get('filter_nm_min')
    maximum = b.get('filter_nm_max')
    result = None
    if minimum is not None and maximum is not None:
        result = (maximum - minimum) / 2
    return result
