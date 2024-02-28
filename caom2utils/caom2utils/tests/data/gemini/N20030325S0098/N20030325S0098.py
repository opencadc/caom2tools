def get_exposure(uri):
    return 45.003


def get_time_delta(uri):
    result = None
    exptime = get_exposure(0)
    if exptime:
        result = exptime / (24.0 * 3600.0)
    return result
