####

import os
#import glob
from glob import glob
from collections import namedtuple

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import wradlib
import wradlib.adjust as adjust
import wradlib.verify as verify
import wradlib.util as util

import datetime
#from datetime import timedelta
import time

import netCDF4 as nc
from netCDF4 import Dataset
from netCDF4 import num2date, date2num

# For progress bars
import ipywidgets
from IPython.core.display import display

# To supress warnings
import warnings
warnings.filterwarnings('ignore')

def get_monthly_data(raw_data_root_path, year_str, month_str):
    fn_list = []
    
    root_path=raw_data_root_path
    month_path = os.path.join(root_path, 
                              year_str, 
                              #(year_str + '_no_clutter_v7-0'),
                              (year_str + '-' + month_str))
    print  month_path
    for day in range(1,32):
        day_dir_str = '%s-%s-%02d' % (year_str, month_str, day)
        day_path = os.path.join(month_path, day_dir_str)
        #day_npz_path = os.path.join(day_path, day_dir_str) # + '_npz')
        fn_list = fn_list + glob(os.path.join(day_path, '*.azi'))
    fn_list.sort()
    return fn_list

def get_dayly_data(year_str, month_str, day):
    fn_list = []
    
    root_path=raw_data_root_path
    month_path = os.path.join(root_path, 
                              (year_str + '_no_clutter_v7-0'),
                              (year_str + '-' + month_str))
    print  month_path
    #for day in range(1,32):
    day_dir_str = '%s-%s-%02d' % (year_str, month_str, day)
    #day_dir_str = '%s-%s-%s-' % (year_str, month_str, day)
    print  day_dir_str
    day_path = os.path.join(month_path, day_dir_str)
    day_npz_path = os.path.join(day_path, day_dir_str + '_npz')
    fn_list = fn_list + glob(os.path.join(day_npz_path, '*.npz'))
    fn_list.sort()
    return fn_list, day_dir_str

def read_npz_data(fn, range_max):
    z_tmp = np.load(fn)
    XML_data_tmp = z_tmp['XML_data']
    data_tmp = z_tmp['data'][:,:range_max]
    r = z_tmp['r'][:range_max]
    azi = z_tmp['azi']
    
    file_name=fn.rpartition('\\')
    file_name=file_name[2].rpartition('.')
    file_name_RS=file_name[0]
    #print (file_name)
    time_file_RS=int(float(file_name_RS[:12]))
    time_temp = datetime.datetime.strptime(str(time_file_RS),'%Y%m%d%H%M')
    time_RS = date2num(time_temp,units='hours since 2000-01-01 00:50:00.0',calendar='standard')
    
    npz_data = namedtuple('npz_data', ['date', 'dBZ', 'r', 'azi', 'XML']) 
    
    return npz_data(date=time_RS, dBZ=data_tmp, r=r, azi=azi, XML=XML_data_tmp)

def read_fn_list(fn_list, range_max):
    data_list = []
    pbar = ipywidgets.FloatProgress(min=0, max=len(fn_list))
    display(pbar)
    print "pbar"

    for i, fn in enumerate(fn_list):
        pbar.value = i
        #try:
        npz_data = read_npz_data(fn, range_max)
        data_list.append(npz_data)
        #except MemoryError:
            #print 'MemoryError for %s' % fn
    return data_list

def read_azi_list(fn_list, range_max):
    XML_data=[]
    data=[]
    r=[]
    azi=[] 
    data_list = []
    XML_list = []
    time_list = []
    
    pbar = ipywidgets.FloatProgress(min=0, max=len(fn_list)-1)
    display(pbar)
    #print "pbar"

    for i, fn in enumerate(fn_list):
        pbar.value = i
        #try:
        #azi_data = read_npz_data(fn, range_max)
        XML_data, azi_data, r, azi = rainbow(fn, XML_data, data, r, azi)
        data_list.append(azi_data)
        XML_list.append(XML_data)
        time = XML_data['date']+' '+XML_data['endtime_round']
        time_list.append(time)
        #except MemoryError:
            #print 'MemoryError for %s' % fn
    return data_list, XML_list, time_list

def rainbow(file_name, XML_data, data, r, azi):

    #file_name="C:\\Users\\bialas-j\\Documents\\standard-gap.azi\\2015\\2015-01\\2015-01-31\\2015013121200000dBuZ.azi"

    # load rainbow file contents to dict
    rbdict = wradlib.io.read_Rainbow(file_name)

    # get azimuthal data
    azi = rbdict['volume']['scan']['slice']['slicedata']['rayinfo']['data']
    azidepth =  float(rbdict['volume']['scan']['slice']['slicedata']['rayinfo']['@depth'])
    azirange = float(rbdict['volume']['scan']['slice']['slicedata']['rayinfo']['@rays'])
    anglestep = float(rbdict['volume']['scan']['slice']['anglestep'])
    azi =  (azi * azirange / 2**azidepth) * anglestep
    ## roll (rotate) raw data matrix to adjust for azimuth offset (aligne data with N-S direction)
    theta1=np.argmin(azi)
    azi=np.roll(azi, -theta1, axis=0)
    #azi=np.arange(0,360,2)
    astart=0.0
    # create range array
    stoprange = float(rbdict['volume']['scan']['slice']['stoprange'])
    rangestep = float(rbdict['volume']['scan']['slice']['rangestep'])
    r = np.arange(0,stoprange,rangestep)
    #r=np.arange(0,500,1)

    # get reflectivity data
    data = rbdict['volume']['scan']['slice']['slicedata']['rawdata']['data']
    ## roll (rotate) raw data matrix to adjust for azimuth offset (aligne data with N-S direction)
    data=np.roll(data, -theta1, axis=0)
    datadepth = float(rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@depth'])
    datamin = float(rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@min'])
    datamax = float(rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@max'])
    
    data = datamin + data * (datamax - datamin) / 2**datadepth
    ## calculate the value of original nodata=255 
    nodata_x = datamin + 255 * (datamax - datamin) / 2**datadepth
    ## change nodata_x to new nodata= -99.99
    #nodata = 255
    nodata = -99.99
    data[data == nodata_x] = nodata

    # get annotation data
    unit = rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@type']
    time = rbdict['volume']['scan']['slice']['slicedata']['@time']
    date = rbdict['volume']['scan']['slice']['slicedata']['@date']
    endtime = rbdict['volume']['@datetime']
    endtime = endtime.rpartition('T')
    endtime = endtime[2]
    starttime = time
    endtime_round = rbdict['volume']['scan']['@time']
    enddate = rbdict['volume']['scan']['@date']
    #highprf = rbdict['volume']['scan']['pargroup']['highprf']
    highprf = float(rbdict['volume']['scan']['slice']['highprf'])
    #lowprf = rbdict['volume']['scan']['pargroup']['lowprf']
    rpm = float(rbdict['volume']['scan']['slice']['antspeed'])
#    lon = rbdict['volume']['sensorinfo']['lon']
#    lat = rbdict['volume']['sensorinfo']['lat']
#    wavelength = rbdict['volume']['sensorinfo']['wavelen']
#    beamwidth = rbdict['volume']['sensorinfo']['beamwidth']
    #sensortype = rbdict['volume']['sensorinfo']['@type']
    #sensorname = rbdict['volume']['sensorinfo']['@name']
#    sensor_id = rbdict['volume']['sensorinfo']['@id']
    max_dynz=float(datamax)
    min_dynz=float(datamin)
    rays_XML = int(rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@rays'])
    bins_XML = int(rbdict['volume']['scan']['slice']['slicedata']['rawdata']['@bins'])
    noise_power_dbz = float(rbdict['volume']['scan']['slice']['noise_power_dbz'])
    log_XML = float(rbdict['volume']['scan']['slice']['log'])
    
    lon=11.028760
    lat=47.728040
    #alt=rbdict['volume']['sensorinfo']['alt']
    alt=955.000000
    wavelength=0.05
    beamwidth=1.0
    #elevation=float(rbdict['volume']['scan']['slice']['posangle'])
    elevation=99.99  ## include check for correct value vs. date time
    sensor_id = "GAP"
    quantity="TH"
    gain = 1.0
    offset=0.0
    undetected = min_dynz  ## -31.5
    rscale = 100.0
    rstart = 0.0
    object_ODIM = "SCAN"
    source = "WMO:0,PLC:GAP"
    #version = "H5rad2.1"
    version = "H5rad 2.3.1"
    product = "SCAN"
    software = "RAINBOW"

    XML_data={"max_dynz":max_dynz,"min_dynz":min_dynz, "rays_XML":rays_XML, "bins_XML":bins_XML,
              "date":date, "enddate": enddate, "time":time, "endtime":endtime,"endtime_round":endtime_round,"lat":lat, "lon":lon, "altitude":alt, "elevation":elevation,
              "system":sensor_id, "quantity":quantity, "gain":gain, "nodata":nodata, "undetected":undetected, "rscale":rscale,
              "source":source, "wavelength":wavelength, "beamwidth":beamwidth, "rscale":rscale, "rstart":rstart, "product":product,
              "software":software, "object":object_ODIM, "offset":offset, "highprf":highprf, "rpm":rpm, "astart":astart}
    #print(XML_data)
    
    return XML_data, data, r, azi

def Reference(io_path_Ref):
    z=np.load(io_path_Ref)
    Ref_XML_data=z['XML_clutter']
    Ref_data=z['Ref_dBZ_mean']
    #print(Ref_data)
    #imshow(Ref_data)
    #print "Ref_data shape:", np.shape(Ref_data)
    #Ref_data=wradlib.trafo.idecibel(Ref_data)
 
    return Ref_data