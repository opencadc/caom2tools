from caom2pipe.manage_composable import convert_to_days, make_datetime, to_float

def get_exposure(uri):
    return 45.003


def get_time_delta(uri):
    result = None
    exptime = get_exposure(0)
    if exptime is not None:
        result = convert_to_days(to_float(exptime))
    return result
