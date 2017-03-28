import wradlib as wrl
import numpy as np
import xarray as xr
import pandas as pd

from tqdm import tqdm


def read_azi_files_to_xarray_dataset(fn_list,
                                     elevation,
                                     r=None,
                                     az=None,
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
    radar_location : tuple of floats
        Radar location of the form (latitude, longitude, altitude)

    Returns
    -------

    xarray.DataSet with radar data and metadata

    """

    # Read the first azi file to get the metadata required for deriving
    # latitude, longitude and altitude
    temp_data, temp_metadata = read_azi_file(fn_list[0])

    if r is None:
        r = temp_metadata['r'] * 1e3
    if az is None:
        az = temp_metadata['az']
    if radar_location is None:
        radar_location = (temp_metadata['longitude'],
                          temp_metadata['latitude'],
                          temp_metadata['altitude'])

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


def read_azi_file(fn):
    """ Read one RAINSCANNER rainbow .azi files

    Parameters
    ----------

    fn : string
        Path to azi file

    Returns
    -------

    ?????

    """

    metadata = {}

    # load rainbow file contents to dict
    rb_dict = wrl.io.read_Rainbow(fn)

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
