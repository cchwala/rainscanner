import wradlib as wrl
import numpy as np
import xarray as xr
import pandas as pd
import tarfile

from tqdm import tqdm


def read_azi_tgz_files_to_xarray_dataset(fn_list,
                                         elevation,
                                         r=None,
                                         az=None,
                                         check_N_az=None,
                                         radar_location=None):
    """ Read azi tgz files and parse them into a xarray DataSet

    Parameters
    ----------

    fn_list : list
        List of paths to .azi files
    elevation : float
        Radar beam elevation in degree
    r : array
        Radar beam ranges in meter
    az : array
        Radar beam azimuth in degree from north
    check_N_az : int
        Check number of az values for each file. Skip file if there is a
        mismatch
    radar_location : tuple of floats
        Radar location of the form (latitude, longitude, altitude)

    Returns
    -------

    xarray.DataSet with radar data and metadata

    """

    # Read the first azi file to get the metadata required for deriving
    # latitude, longitude and altitude
    if check_N_az is None:
        with tarfile.open(fn_list[0]) as tar:
            f = tar.extractfile(tar.getmembers()[0])
            r, az, radar_location = _get_r_az_loc(f)
    # If check_N_az is given, iterate over azi files in first tgz file
    # until az is found with the correct length
    else:
        with tarfile.open(fn_list[0]) as tar:
            for file_in_tar in tar.getmembers():
                f = tar.extractfile(file_in_tar)
                r, az, radar_location = _get_r_az_loc(f)
                if len(az) == check_N_az:
                    break
                else:
                    'Trying to get metadata from %s' % f

    # Overwrite r, az, location with arguments if supplied
    if r is not None:
        r = r
    if az is not None:
        az = az
    if radar_location is not None:
        radar_location = radar_location

    # Build 2D grids for r and az
    r_grid, az_grid = np.meshgrid(r, az)

    lons, lats, alts = wrl.georef.polar2lonlatalt_n(r_grid,
                                                    az_grid,
                                                    elevation,
                                                    radar_location)

    data_list = []
    metadata_list = []
    for fn in fn_list:
        with tarfile.open(fn) as tar:
            for tarinfo in tqdm(tar, desc=('Reading ' + fn)):
                f = tar.extractfile(tarinfo)
                temp_data, temp_metadata = read_azi_file(f)
                if check_N_az is not None:
                    if temp_data.shape[0] != check_N_az:
                        print 'N_az = %d instead of %d. --> Skipping %s' % (
                            temp_data.shape[0], check_N_az, fn)
                        continue
                    if len(temp_metadata['az']) != check_N_az:
                        print 'N_az = %d instead of %d. --> Skipping %s' % (
                            len(temp_metadata['az']), check_N_az, fn)
                        continue
                data_list.append(temp_data)
                metadata_list.append(temp_metadata)

    time_list = [pd.to_datetime(metadata['date'] + ' ' + metadata['time'])
                 for metadata in metadata_list]

    ds = xr.Dataset(coords={'r': ('r', r),
                            'az': ('az', az),
                            'time': ('time', time_list),
                            'latitude': (['az', 'r'], lats),
                            'longitude': (['az', 'r'], lons)},
                    data_vars={'ZH_raw': (['time', 'az', 'r'], data_list)})
    return ds


def read_azi_files_to_xarray_dataset(fn_list,
                                     elevation,
                                     tgz_file=False,
                                     r=None,
                                     az=None,
                                     check_N_az=None,
                                     radar_location=None):
    """ Read .azi files and parse them into a xarray DataSet

    Parameters
    ----------

    fn_list : list
        List of paths to .azi files
    elevation : float
        Radar beam elevation in degree
    r : array
        Radar beam ranges in meter
    az : array
        Radar beam azimuth in degree from north
    check_N_az : int
        Check number of az values for each file. Skip file if there is a
        mismatch
    radar_location : tuple of floats
        Radar location of the form (latitude, longitude, altitude)

    Returns
    -------

    xarray.DataSet with radar data and metadata

    """

    # Read the first azi file to get the metadata required for deriving
    # latitude, longitude and altitude
    if check_N_az is None:
        r, az, radar_location = _get_r_az_loc(fn_list[0])
    # If check_N_az is given, iterate over azi files
    # until az is found with the correct length
    else:
        for fn in fn_list:
            r, az, radar_location = _get_r_az_loc(fn)
            if len(az) == check_N_az:
                break
            else:
                'Trying to get metadata from %s' % fn

    # Overwrite r, az, location with arguments if supplied
    if r is not None:
        r = r
    if az is not None:
        az = az
    if radar_location is not None:
        radar_location = radar_location

    # Build 2D grids for r and az
    r_grid, az_grid = np.meshgrid(r, az)

    lons, lats, alts = wrl.georef.polar2lonlatalt_n(r_grid,
                                                    az_grid,
                                                    elevation,
                                                    radar_location)

    data_list = []
    metadata_list = []
    for fn in tqdm(fn_list):
        temp_data, temp_metadata = read_azi_file(fn)
        if check_N_az is not None:
            if temp_data.shape[0] != check_N_az:
                print 'N_az = %d instead of %d. --> Skipping %s' % (
                    temp_data.shape[0], check_N_az, fn)
                continue
            if len(temp_metadata['az']) != check_N_az:
                print 'N_az = %d instead of %d. --> Skipping %s' % (
                    len(temp_metadata['az']), check_N_az, fn)
                continue
        data_list.append(temp_data)
        metadata_list.append(temp_metadata)

    time_list = [pd.to_datetime(metadata['date'] + ' ' + metadata['time'])
                 for metadata in metadata_list]

    ds = xr.Dataset(coords={'r': ('r', r),
                            'az': ('az', az),
                            'time': ('time', time_list),
                            'latitude': (['az', 'r'], lats),
                            'longitude': (['az', 'r'], lons)},
                    data_vars={'ZH_raw': (['time', 'az', 'r'], data_list)})
    return ds


def read_azi_file(file_name_or_handle):
    """ Read one RAINSCANNER rainbow .azi files

    Parameters
    ----------

    file_name_or_handle : string or file handle
        Path to azi file of open file handle

    Returns
    -------

    ?????

    """

    metadata = {}

    # load rainbow file contents to dict
    rb_dict = wrl.io.read_Rainbow(file_name_or_handle)

    slice_dict = rb_dict['volume']['scan']['slice']

    # get azimuthal data
    az_raw = slice_dict['slicedata']['rayinfo']['data']
    az_depth = float(slice_dict['slicedata']['rayinfo']['@depth'])
    az_range = float(slice_dict['slicedata']['rayinfo']['@rays'])
    angle_step = float(slice_dict['anglestep'])
    az_uncorrected = (az_raw * az_range / 2**az_depth) * angle_step

    # roll (rotate) raw data matrix to adjust for azimuth offset
    # (align data with 0 degree azimuth to the north)
    theta1 = np.argmin(az_uncorrected)
    az = np.roll(az_uncorrected, -theta1, axis=0)

    # create range array
    stoprange = float(slice_dict['stoprange'])
    rangestep = float(slice_dict['rangestep'])
    r = np.arange(0, stoprange, rangestep)

    metadata['r'] = r
    metadata['az'] = az

    # get reflectivity data
    data = slice_dict['slicedata']['rawdata'].pop('data')

    # roll (rotate) raw data matrix to adjust for azimuth offset
    # (align data with 0 degree azimuth to the north)
    data = np.roll(data, -theta1, axis=0)

    # scale data
    data_depth = float(slice_dict['slicedata']['rawdata']['@depth'])
    data_min = float(slice_dict['slicedata']['rawdata']['@min'])
    data_max = float(slice_dict['slicedata']['rawdata']['@max'])
    data = data_min + data * (data_max - data_min) / 2**data_depth

    # calculate the value of original nodata=255
    nodata_x = data_min + 255 * (data_max - data_min) / 2**data_depth

    # change nodata_x to new nodata= -99.99
    # nodata = 255
    nodata = -99.99
    data[data == nodata_x] = nodata

    # Parse metadata from `slice`
    for key, value in slice_dict.iteritems():
        if key == 'slicedata':
            metadata[u'time'] = value['@time']
            metadata[u'date'] = value['@date']
            metadata[u'rayinfo'] = value['rayinfo']
            metadata[u'datainfo'] = value['rawdata']
        else:
            metadata[key] = value

    # Parse more general metadata from lower level
    metadata['id'] = rb_dict['volume']['sensorinfo']['@id']
    metadata['name'] = rb_dict['volume']['sensorinfo']['@name']
    metadata['latitude'] = float(rb_dict['volume']['sensorinfo']['lat'])
    metadata['longitude'] = float(rb_dict['volume']['sensorinfo']['lon'])
    metadata['altitude'] = float(rb_dict['volume']['sensorinfo']['alt'])
    metadata['wavelength'] = float(rb_dict['volume']['sensorinfo']['wavelen'])
    metadata['beamwidth'] = float(rb_dict['volume']['sensorinfo']['beamwidth'])

    return data, metadata


def _get_r_az_loc(fn):
    """ Return r, az and radar_location for one file """
    temp_data, temp_metadata = read_azi_file(fn)

    r = temp_metadata['r'] * 1e3
    az = temp_metadata['az']
    radar_location = (temp_metadata['longitude'],
                      temp_metadata['latitude'],
                      temp_metadata['altitude'])
    return r, az, radar_location
