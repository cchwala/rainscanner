####

import os
#import glob
from glob import glob
from collections import namedtuple

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import netCDF4 as nc
from netCDF4 import Dataset
from netCDF4 import num2date, date2num

def make_dir(raw_data_root_path, year, month):
    io_path_res = os.path.join(raw_data_root_path, (year+'_raw_data'))
    io_path_res_lev1 = os.path.join(io_path_res, (year +'-'+month))
    io_path_res_lev2 = os.path.join(io_path_res_lev1, ('Results_'+month))
    #print io_path_res
    #print io_path_res_lev1
    #print io_path_res_lev2
    
    dir_results = os.listdir(raw_data_root_path)
    #print dir_results
    if str(year+'_raw_data') in dir_results:   ##check if dir "path_Results" allready exists
        print "existing directory:", str(year+'_raw_data')
    else:
        os.mkdir(io_path_res)
        print "new directory:", str(year+'_raw_data')
        
    dir_results = os.listdir(io_path_res)
    #print dir_results
    if str(year +'-'+month) in dir_results:   ##check if dir "path_Results" allready exists
        print "existing directory:", str(year +'-'+month)
    else:
        os.mkdir(io_path_res_lev1)
        print "new directory:", str(year +'-'+month)
    
    dir_results = os.listdir(io_path_res_lev1)
    #print dir_results
    if str('Results_'+month) in dir_results:   ##check if dir "path_Results" allready exists
        print "existing directory:", str('Results_'+month)
    else:
        os.mkdir(io_path_res_lev2)
        print "new directory:", str('Results_'+month)
        
    return io_path_res, io_path_res_lev1, io_path_res_lev2

def scan_errors(a, b, c, thr_empty, file_scheck, filename_scheck, filename_scheck_txt, io_path_res_lev2):
    a0 = ["Version:  "+"preliminary"," ","Sector error"," "]
    b0 = [" "," ","Empty scan, thresh.="+str(thr_empty)," "]
    #c0 = [" "," ","Date mismach, max.delay[s]="+str(max_delay)," "]
    c0 = [" "," ","Date mismach disabled"," "]
    
    dir_results_files = os.listdir(io_path_res_lev2)
    
    if (file_scheck+'.npz') in dir_results_files:    ##check if file "Log_sanity_checks" allready exists
        print "File existing:", str(file_scheck)
        temp_dat = np.load(filename_scheck)
        a0 = temp_dat['A']
        b0 = temp_dat['B']
        c0 = temp_dat['C']
    else:
        #init_s_check_file(filename_scheck, Version_IFU, thr2, max_delay)
        print "File new:", str(file_scheck)
    
    #a, b, c = load_s_check_file(filename_scheck)
    

    A = np.hstack((a0, a))
    B = np.hstack((b0, b))
    C = np.hstack((c0, c))
    
    np.savez_compressed(filename_scheck, A=A, B=B, C=C)  ## save .npz
    
    a_l = np.array(A).tolist()
    b_l = np.array(B).tolist()
    c_l = np.array(C).tolist()
    
    max_length = max(len(a_l), len(b_l), len(c_l))
    
    row_a = (a_l+["                   -"]*(max_length))[:max_length]
    row_b = (b_l+["                   -"]*(max_length))[:max_length]
    row_c = (c_l+["                   -"]*(max_length))[:max_length]

    file = open(filename_scheck_txt, "w")
    for index in range(max_length):
        file.write(str(row_a[index]) + "                 " + str(row_b[index]) + "                  " + str(row_c[index])+"\n")
    
    file.close()

    