
# coding: utf-8

# In[100]:


from IPython.display import HTML, display
HTML('''<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" value="Click here to toggle on/off the raw code."></form>''')


# In[1]:


import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import math
from scipy import stats
import glob
from scipy.interpolate import interp1d
import mplcursors
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib.widgets import Slider, Button, RadioButtons
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets

import IPython
from astropy import constants as const
from astropy import units as u


# In[2]:


ML_disk, ML_bulge, DM_frac, SB_bary, Ltot, HImass, Vflat, index_stop = np.loadtxt('data for good galaxies.csv', delimiter=',', skiprows = 1, usecols = (1,2,3,5,14,20,22,24),unpack = True)
Lbulge  = np.loadtxt('sparc global measurements.csv', delimiter=',', unpack = True, skiprows = 2, usecols = (20)) 
sparc = np.load('sparc_dict.npy').item()
name = np.loadtxt('starkman data.csv', delimiter=',', skiprows = 2, usecols = (0), dtype = 'object', unpack = True)
Mbary = np.loadtxt('starkman data.csv', delimiter=',', skiprows = 2, usecols = (19), unpack = True)
rad_start, rad_stop = np.loadtxt('bracketing curves.csv', delimiter=',', skiprows = 1, usecols = (3,4),unpack = True)

Ltot, Lbulge, HImass , Mbary = Ltot*1e9, Lbulge*1e9, HImass*1e9, Mbary*1e9


# In[3]:


Vflat_arr = np.empty(len(name))
VflatML_arr = np.empty(len(name))


# In[11]:


def Vflat_recurse(index_next, Vbar_sum, vobs):
    
    Vbar_sum_new = Vbar_sum + vobs[-1*index_next]
    Vbar_new = Vbar_sum_new/index_next

    if ((index_next == len(vobs)-1) | (np.isnan(vobs[-1*(index_next+1)]))):
        return Vbar_new
    else:
        if (np.abs(vobs[-1*(index_next+1)] - Vbar_new)/Vbar_new <= 0.05):
            return Vflat_recurse(index_next+1, Vbar_sum_new,vobs)
        else:
            return Vbar_new

for x in range(len(name)):
    vbary = np.sqrt(0.5 * sparc[name[x]+"_Vdisk"]**2 + 0.7 * sparc[name[x]+"_Vbul"]**2 + sparc[name[x]+"_Vgas"]**2 )
    vobs = sparc[name[x]+'_Vobs']
    sparc[name[x] + '_VDM'] = np.sqrt(vobs**2 - vbary**2)
    
    #print(sparc[name[x] + '_VDM'])
    
    Vflat_orig = 0.5 * (vobs[-1] + vobs[-2])
    Vbar_sum = vobs[-1] + vobs[-2]
    
    if (np.abs(vobs[-3] - Vflat_orig)/Vflat_orig <= 0.05):
        Vflat_arr[x] = Vflat_recurse(3, Vbar_sum, vobs)
    else:
        Vflat_arr[x] = 0


# In[12]:


plt.figure(figsize=(8,6))
plt.scatter(Vflat_arr, Vflat_arr/Vflat)
plt.xlabel('my Vflat [km/s]')
plt.ylabel(r"$\mathrm{\frac{my \ \ Vflat}{Fredrico's \ \ Vflat}}$", rotation = 0, fontsize = 15, horizontalalignment = 'right')
plt.yticks(rotation=0)
#print(Vflat_arr)


# # now changing M/L

# In[13]:


def changing_ML(mult_factor, gal_name, x):
    ML_disknew = 0.5 * mult_factor
    ML_bulgenew = 0.7 * mult_factor
    
    vbaryorig = np.sqrt(.5 * sparc[gal_name+"_Vdisk"]**2 + .7 * sparc[gal_name+"_Vbul"]**2 + sparc[gal_name+"_Vgas"]**2)
    vobsorig = np.sqrt(vbaryorig**2 + sparc[gal_name+"_VDM"]**2)
    
    vbarynew = np.sqrt(ML_disknew * sparc[gal_name+"_Vdisk"]**2 + ML_bulgenew * sparc[gal_name+"_Vbul"]**2 + sparc[gal_name+"_Vgas"]**2)
    vobsnew = np.sqrt(vbarynew**2 + sparc[gal_name+"_VDM"]**2)
    
    
    Vflat_orig = 0.5 * (vobsnew[-1] + vobsnew[-2])
    Vbar_sum = vobsnew[-1] + vobsnew[-2]
    
    if (np.abs(vobsnew[-3] - Vflat_orig)/Vflat_orig <= 0.05):
        VflatML_arr[x] = Vflat_recurse(3, Vbar_sum, vobsnew)
    else:
        VflatML_arr[x] = 0
        
    #print(name[x], VflatML_arr[x])
    
#     plt.plot(sparc[gal_name+"_rad"], vobsnew, label = 'new vobs')
#     plt.plot(sparc[gal_name+"_rad"], vobsorig, ls = ':', color = 'red', label = 'original vobs')
#     plt.title(gal_name)
#     plt.legend()
#     #plt.figtext(.1, .05, '{:.3f}, {}'.format(VflatML_arr[x], vobsnew))
#     plt.figtext(.1, .05, 'Vflat = {:.1f} km/s'.format(VflatML_arr[x]))
#     #plt.savefig('test galaxies/'+str(gal_name) + '.png', dpi=300)
#     plt.close()


def graphing_ML(mult_factor):
    for x in range(len(name)):
        changing_ML(mult_factor, name[x], x)

    nonzero = np.where(VflatML_arr > 0)

    plt.figure(figsize=(10,8))
    plt.scatter(Vflat[nonzero], VflatML_arr[nonzero] / Vflat[nonzero])
    plt.xlabel("Fredrico's Vflat [km/s]")
    plt.ylabel(r"$\mathrm{\frac{my \ Vflat}{Fredrico's \ \ Vflat}}$", rotation = 0, fontsize = 15, horizontalalignment = 'right')
    plt.yticks(rotation=0)
    plt.ylim(.8, 1.8)

#interact(graphing_ML, mult_factor=widgets.FloatSlider(value=2, min=0.01, max=5,step=0.01, continuous_update=False)) 


# # changing M/L by shifting bulge and disk individually

# In[14]:


def changing_ML2(ML_disknew, ML_bulgenew, gal_name):
    
    vbaryorig = np.sqrt(.5 * sparc[gal_name+"_Vdisk"]**2 + .7 * sparc[gal_name+"_Vbul"]**2 + sparc[gal_name+"_Vgas"]**2)
    vobsorig = np.sqrt(vbaryorig**2 + sparc[gal_name+"_VDM"]**2)
    
    vbarynew = np.sqrt(ML_disknew * sparc[gal_name+"_Vdisk"]**2 + ML_bulgenew * sparc[gal_name+"_Vbul"]**2 + sparc[gal_name+"_Vgas"]**2)
    vobsnew = np.sqrt(vbarynew**2 + sparc[gal_name+"_VDM"]**2)
    
    
    Vflat_orig = 0.5 * (vobsnew[-1] + vobsnew[-2])
    Vbar_sum = vobsnew[-1] + vobsnew[-2]
    
    if (np.abs(vobsnew[-3] - Vflat_orig)/Vflat_orig <= 0.05):
        return Vflat_recurse(3, Vbar_sum, vobsnew)
    else:
        return 0
        
    #print(name[x], VflatML_arr[x])
    
#     plt.plot(sparc[gal_name+"_rad"], vobsnew, label = 'new vobs')
#     plt.plot(sparc[gal_name+"_rad"], vobsorig, ls = ':', color = 'red', label = 'original vobs')
#     plt.title(gal_name)
#     plt.legend()
#     #plt.figtext(.1, .05, '{:.3f}, {}'.format(VflatML_arr[x], vobsnew))
#     plt.figtext(.1, .05, 'Vflat = {:.1f} km/s'.format(VflatML_arr[x]))
#     #plt.savefig('test galaxies/'+str(gal_name) + '.png', dpi=300)
#     plt.close()

def graphing_ML2(ML_disk, ML_bulge):
    for x in range(len(name)):
        VflatML_arr[x] = changing_ML2(ML_disk, ML_bulge, name[x])

    nonzero = np.where(VflatML_arr > 0)

   #  plt.figure(figsize=(10,8))
#     plt.scatter(Vflat[nonzero], VflatML_arr[nonzero] / Vflat[nonzero])
#     plt.xlabel("Fredrico's Vflat [km/s]")
#     plt.ylabel(r"$\mathrm{\frac{my \ Vflat}{Fredrico's \ \ Vflat}}$", rotation = 0, fontsize = 15, horizontalalignment = 'right')
#     plt.yticks(rotation=0)
#     plt.ylim(.8, 2)

#interact(graphing_ML2, ML_disk=widgets.FloatSlider(value=0.5, min=0.01, max=5,step=0.01, continuous_update=False),
      #  ML_bulge=widgets.FloatSlider(value=0.7, min=0.01, max=5,step=0.01, continuous_update=False)) 


# In[112]:


#for x in range(len(name)):
   # print(name[x], Vflat[x])

