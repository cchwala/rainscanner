
# coding: utf-8

# In[ ]:
import numpy as np
import wradlib
import pylab as plt
    
def polar_plot(data, r, az, title, units, subplot):
    vmax = np.nanmax(data)
    vmin = np.nanmin(data)
    cgax, caax, paax, pm = wradlib.vis.plot_cg_ppi(data, r=r, az=az, vmax=vmax, vmin=vmin, subplot=subplot)
    #cgax, caax, paax, pm = wradlib.vis.plot_cg_ppi(data, r=r, az=az, vmax=vmax, vmin=vmin)

    #pm.set_cmap('spectral')
    pm.set_cmap('jet')
    title = title
    t = plt.title(title, fontsize=10)
    t.set_y(1.1)
    #cbar = plt.colorbar(pm, shrink=0.75, pad=0.085)
    #cbar = plt.colorbar(pm, shrink=shrink, pad=pad, fraction=fraction)
    cbar = plt.colorbar(pm)
    #caax.set_xlabel('x_range')
    #caax.set_ylabel('y_range')
    caax.set_xlabel('km')
    caax.set_ylabel('km')
    #plt.text(1.0, 1.05, 'azimuth', transform=caax.transAxes, va='bottom',
    #        ha='right')
    cbar.set_label(units)
    plt.tight_layout()

def histogram(data_show2, Max_b, Min_b, width_of_bin, title, units, fig, subplot):
    
    Max = np.nanmax(data_show2)
    Min = np.nanmin(data_show2)
    Mean = np.nanmean(data_show2)
    Median = np.nanmedian(data_show2)

    #nr_of_bins = int((Max-Min)/width_of_bin)
    nr_of_bins = int((Max_b-Min_b)/width_of_bin)
    #hist, bins = np.histogram(data_show2, bins=nr_of_bins, range=(Min_,Max))
    hist, bins = np.histogram(data_show2, bins=nr_of_bins, range=(Min_b,Max_b))
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:])/2
    
    title_ext = '\n'+'Max = '+str(Max)[:5]+" ;"+'  Min = '+str(Min)[:5]+" ;"+'  Mean = '+str(Mean)[:5]+" ;"+'  Median = '+str(Median)[:5]+" "+'\n'+"width_of_bin = "+str(width_of_bin)+'\n'
    
    #fig, ax = plt.subplots()
    ax = fig.add_subplot(subplot) #, aspect='equal')
    ax.bar(center[1:], hist[1:], align='center', width=width)
    ax.set_title(title+title_ext, fontsize=10)
    ax.set_xlabel(units)
    
    #ext_filename3 = "Histogram " + ext_filename2
    #saveto=path_png + '\\' + ext_filename3 + '.png'
    #fig.savefig(saveto, dpi=130)
    #pl.close()
    
    return Max, Min, Mean, Median