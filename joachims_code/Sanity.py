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

import scipy
#from scipy import ndimage
from scipy import signal, fftpack
from scipy.interpolate import UnivariateSpline



def rays_correct(XML_data, correct_ray_number):
    '''
    Check if the number of rays equals a desired number
    '''

    number_is_correct = (np.abs(XML_data["rays_XML"] - correct_ray_number) == 0)
    
    return number_is_correct


def scan_empty(data, thr_empty):
    ''' 
    Check for "empty" scans by checking dBZ values 
    at a known clutter location.

    data = radar dBz matrix
    range_bin_range = range which is used for the check. 

    '''
    min_clutter_dBZ_threshold=thr_empty
    
    rangebin_range1=np.arange(98,108)  ## Mount Hoernle near Bad Kohlgrub
    max_in_rangebin_range1 = np.max(data[:, rangebin_range1])
    
    rangebin_range2=np.arange(230,240)  ## Mount Krottenkopf (Estergebirge)
    max_in_rangebin_range2 = np.max(data[:, rangebin_range2])
    
    rangebin_range3=np.arange(340,350)  ## Mount Zugspitze
    max_in_rangebin_range3 = np.max(data[:, rangebin_range3])
    
    rangebin_range4=np.arange(20,30)  ## near radar
    max_in_rangebin_range4 = np.max(data[:, rangebin_range4])
    
    scan_is_empty = (max_in_rangebin_range1 + max_in_rangebin_range2 + max_in_rangebin_range3+ max_in_rangebin_range4)< min_clutter_dBZ_threshold
    
    return scan_is_empty, thr_empty


def clean_up(data, XML_data, time, thr_empty, ray_number, data_list_clean, XML_list_clean, time_list_clean, a, b):
    
    #### check_1 for sector errors
    number_is_correct = rays_correct(XML_data, ray_number)

    #### check_2 for "empty" scans
    ## thr2: dB threshold to remove "empty" scans (@ 3 range slices)
    scan_is_empty, thr_empty = scan_empty(data, thr_empty)

    #### check_3 for date mismach of the scan
    #date_is_correct, total_seconds = date_mismach(XML_data, max_delay)

    if number_is_correct==False:
        # a.append(time)  ## applicable for list only
        a=np.hstack((a, time))
    elif scan_is_empty==True:
        #b.append(file_name)
        b=np.hstack((b, time))
#            elif date_is_correct==False:
#                #c.append(time)
#                c=np.hstack((c, time))
    else:
        ## add data to XML_dictionary
        w0 = {"thr_empty":thr_empty}
        XML_data.update(w0)
        #w1 = {"max_time_delay":total_seconds}
        #XML_data.update(w1)
        data_list_clean.append(data)
        XML_list_clean.append(XML_data)
        time_list_clean.append(time)
    
    return a, b, data_list_clean, XML_list_clean, time_list_clean

def rfft_xcorr(x, y):
    ### Algoryth by user "eryksun" http://stackoverflow.com/questions/4688715/find-time-shift-between-two-similar-waveforms 
    ### It uses rfft and zero pads the inputs to a power of 2 large enough to ensure linear (i.e. non-circular) correlation:
    M = len(x) + len(y) - 1
    N = 2 ** int(np.ceil(np.log2(M)))
    X = np.fft.rfft(x, N)
    Y = np.fft.rfft(y, N)
    cxy = np.fft.irfft(X * np.conj(Y))
    cxy = np.hstack((cxy[:len(x)], cxy[N-len(y)+1:]))
    return cxy

def match(x, ref):
    ### To match a reference signal, compute rfft_xcorr(x, ref) and search for the peak.
    cxy = rfft_xcorr(x, ref)
    index = np.argmax(cxy)
    if index < len(x):
        return index
    else: # negative lag
        return index - len(cxy)  

def check_azimuth(Ref_data, data, ranges, delta_spl, step, azi_start, k):
    ## roll (rotate) raw data matrix to adjust for azimuth offset (aligne data with N-S direction) 
    
    ## rotate data for test
    #corr=9
    #data=np.roll(data, corr, axis=0)
    #data=np.roll(Ref_data, corr, axis=0) ## for test with Ref_data only
    #imshow(data)

    ## noise filtering
    #data= ndimage.median_filter(data, 1)
    #data= ndimage.gaussian_filter(data, 1)
    #imshow(data)
    
    thetas=[]  ## azimuth offset calclated for different ranges
    for range in ranges:
        ## define data slice for convolution (from 0 - 499)
        range1 = range-delta_spl
        range2 = range1+delta_spl*2
        #print "delta_spl:", delta_spl
        #print "range2 - range1 =", (range2 - range1)
        
        ## prepair Ref_data
        Ref_data_temp = Ref_data[:][:,range1:range2]
        #print "Ref_data_temp shape:", Ref_data_temp.shape
        #Ref_data_temp=Ref_data_temp.flatten()
        Ref_data_temp = Ref_data_temp.flatten('F')
        #print(Ref_data_temp)
        #print("---------------")
        #### convert to linear units
        Ref_data_temp = wradlib.trafo.idecibel(Ref_data_temp)
        #### make spline to the Ref_data
        x0 = np.arange(0, len(Ref_data_temp), 1)
        xs = np.arange(0, len(Ref_data_temp), step)
        Ref_data_temp_spl = UnivariateSpline(x0, Ref_data_temp, k=k, s=None) #default s=None
        #plt.plot(xs, Ref_data_temp_spl(xs), 'g', lw=1)
        #print "Ref_data_temp_spl shape:", Ref_data_temp_spl(x0).shape
        
        ## prepair scan data
        data_temp = data[:][:,range1:range2]
        data_temp = data_temp.flatten('F')
        data_temp = wradlib.trafo.idecibel(data_temp)
        
        #### make spline to the data
        data_temp_spl = UnivariateSpline(x0, data_temp, k=k, s=None) # default s=None
        #plt.plot(xs, data_temp_spl(xs), 'b', lw=1)
        #print "data_temp_spl shape:", data_temp_spl(x0).shape
    
        ## calculate rotation angle    
        #crr=match(data_temp, Ref_data_temp)
        crr = match(data_temp_spl(xs), Ref_data_temp_spl(xs))*step 
        signe = (crr+0.000000001)/abs(crr+0.000000001)
        #print('signe=',signe)
        res1 = divmod(crr,180)
        #print(res1[1])
        res2 = divmod(crr,-180)
        #print(res2[1])
        theta2x = min(abs(res1[1]),abs(res2[1]))*signe
        
        thetas = np.hstack((thetas,theta2x))

    ## check for correct azimuths
    becon_0 = np.round(thetas[0], decimals=1)
    becon_1 = np.round(thetas[1], decimals=1)
    becon_2 = np.round(thetas[2], decimals=1)
    becon_3 = np.round(thetas[3], decimals=1)
    becon_4 = np.round(thetas[4], decimals=1)
    #print "becon_x:", becon_0, becon_1, becon_2, becon_3, becon_4

    becon = thetas
    becon = sorted(becon)
    #print becon
    corr_azi_old = azi_start #corr_azi

    if abs(becon[1]-becon[2])<1.01 and abs(becon[3]-becon[2])<1.01:
        corr_azi_new = int(np.round((becon[1]+becon[2]+becon[3])/3., decimals=0))
        
        ## check to supress noise
        if abs(corr_azi_old-corr_azi_new) < 2:  # ignor noise (+/-1)
            corr_azi = corr_azi_old 
            #print corr_azi
        else: 
            corr_azi = corr_azi_new
    else:
        corr_azi = corr_azi_old
        
    #print "corr_azi=", corr_azi
    
    return corr_azi, thetas 


def aligne_azimuth(Ref_data, data, XML_data, ranges, delta_spl, step, azi_start, k, 
                   data_list_aligned, XML_list_aligned, thetas_list_thetas):
    #### adjust azimuth to reference
    corr_azi, thetas = check_azimuth(Ref_data, data, ranges, delta_spl, step, azi_start, k)

    w2 = {"corr_azi":corr_azi}
    XML_data.update(w2)

    ### aligne data for ODIM_H5
    #corr_azi=0  ## for test reason only!
    data = np.roll(data, -corr_azi, axis=0)
    
    data_list_aligned.append(data)
    XML_list_aligned.append(XML_data)
    thetas_list_thetas.append(thetas)
    
    return data_list_aligned, XML_list_aligned, thetas_list_thetas
 