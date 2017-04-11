"""
PlotRates_CRESST_Earth.py

Plot Mollweide projected Earth Map of CRESST-II Rate

Last Modified: 11/04/2017
Contact: Bradley Kavanagh (bradkav@gmail.com)

"""

#------------------------------------------
# Importing modules and setting up 
# matplotlib parameters
#------------------------------------------
import numpy as np
from scipy.interpolate import interp1d
from EarthAngles import calcGamma
import matplotlib.colors as colors

#Load the basemap module!
from mpl_toolkits.basemap import Basemap

import matplotlib as mpl
mpl.rc('text', usetex = True)
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
import pylab as pl
pl.rcParams.update(params)

font = { 'size'   : 15, 'family': 'serif'}
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

import os.path
import sys

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#------------------------------------------
# Options
#------------------------------------------

#Specify which operator
oper = 1

#Specify the value of alpha (angle between Earth's motion and z-axis)
alpha = (0.5*(36.3+49.3))

#Time in hours (t = 0 is some arbitrary time...)
time = 6


#------------------------------------------
# Main Calc/Plotting Code
#------------------------------------------

#Set up the map - molleweide projection looking down on 0N, 0W
#See http://matplotlib.org/basemap/users/examples.html for more examples
fig, ax = pl.subplots()

map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='c')
map.drawcoastlines(linewidth=0.5)
map.drawparallels(np.arange(-90.,120.,30.), labels=[1, 0, 0, 0],fontsize=14)
map.drawmeridians(np.arange(0.,420.,60.))

#NB: if you want to focus on a specific part of the Earth
# then the 'ortho' projection (orthographic) is useful.


#Load results from file and generate interpolation (as a function of angle gamma)
data = np.loadtxt("../results/CRESST-II/CRESSTRate_op=" + str(oper)+"_mx=0.50_ps=0.10.dat")
Rate = interp1d(np.pi*data[:,0], data[:,3], kind='linear')

#Calculate a regular lat/lon grid.
nlats = 145; nlons = 145; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (-np.pi + delta*np.indices((nlats,nlons))[1,:,:])

#Compute map projection coordinates of lat/lon grid.
x, y = map(lons*180./np.pi, lats*180./np.pi)


#Calculate the value of gamma at each point on the map
gamma = lats*0.0
for ilat in range(nlats):
    for ilon in range(nlons):
        gamma[ilat, ilon] = calcGamma(alpha*np.pi/180.0,lats[ilat,ilon],time + 24*lons[ilat,ilon]/(2.0*np.pi))


#Calculate "overhead" point and plot it as a cross
xa, ya = map(-time*(180.0/12.0), alpha)
pointa, = map.plot(xa,ya, 'k+', markersize=15.0, mew=2.0)

#Set up the colormap
maxlim = 1.2
minlim = 0.8
cm = pl.get_cmap("seismic")
new_cmap = truncate_colormap(cm, minval=0.05, 
        maxval=0.95, n=100)

#Plot the colormap of the rate
cshow = map.pcolor(x, y, Rate(gamma),
            cmap=new_cmap, vmin=minlim, vmax=maxlim, linewidths=0)

cb = map.colorbar(cshow,"right", extend='both')


#Labels for the plot
pl.title(r'$\mathrm{Operator}$ $\hat{\mathcal{O}}_{' + str(int(oper))+'}$ - $m_\chi = 0.5 \,\,\mathrm{GeV}$',
                y=1.03, fontsize=16)

timelabel = ax.text(0.82, 0.94, r't = ' + str(time) + ' hr', fontsize=14, 
            ha='left', va='center', transform=ax.transAxes)

ax.text(0.5, -0.07, r'Relative rate enhancement due to Earth-scattering \textit{(attenuation + deflection)}', fontsize=12, 
            ha='center', va='center', transform=ax.transAxes)

#pl.savefig('../plots/EarthMap_ortho_O' +str(oper)+'.pdf', bbox_inches='tight')
pl.show()

