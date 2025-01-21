def get_proposal_id(uri):
    return 'PROPOSAL ID'


def get_energy_chunk_range_end_val(parameters):
    header = parameters.get('header')
    if (camera_arm := header.get('ARM')) == 'RED':
        result = 920.0  # nm
    elif camera_arm == 'BLUE':
        result = 663.0  # nm
    return result


def get_energy_chunk_range_start_val(parameters):
    header = parameters.get('header')
    result = None
    if (camera_arm := header.get('ARM')) == 'RED':
        result = 649.0  # nm
    elif camera_arm == 'BLUE':
        result = 499.0  # nm
    return result


def get_exposure(header):
    return 4.0


def get_time_delta(parameters):
    header = parameters.get('header')
    result = None
    exptime = get_exposure(header)
    if exptime:
        result = 12.0
    return result


def get_time_function_val(parameters):
    header = parameters.get('header')
    return header.get('MJD')

