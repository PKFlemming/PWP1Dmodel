import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import time_splitters
from datetime import datetime

datafile = r"D:\UKuni\3rdYr\MEP\cw\nc\AS\MLmethod\dsig\AS_dsig_reup_all_1.nc"
ds = nc.Dataset(datafile)
tempdata = ds["temp"]       # get temp profile data
zs = ds["depth"]            # depth data
jdays = ds["time"]          # get julian day data for all profiles
ms = ds["time2"]            # time2 tells us how many milliseconds we are into the day for each profile
jtimepairs = zip(jdays, ms) # zip jday and ms into a tuple for each profile (effectively [day, time])
ts = [time_splitters.JDayToGreg(jtimepair[0], jtimepair[1]) for jtimepair in jtimepairs] # convert to gregorian calendar


depthlimit = 80             # set depth limit to display
# make contoured graphic
fig = plt.figure(figsize=(24, 4))
heat_ax = plt.subplot2grid((1, 6), (0, 0), colspan=6)
cmap = plt.cm.rainbow

studyarea = "AS"

if studyarea == "COARE":
    tmin = 28                 # set max and min temperatures for colour bar (28.5,30.5 for COARE; 24,32 for AS)
    tmax = 30
    levels = np.linspace(tmin, tmax, 25) # COARE
    time_min = datetime(1994, 11, 5)
    myformat = '%.2f'

if studyarea == "AS":
    tmin = 24
    tmax = 32
    levels = np.linspace(tmin, tmax, 25)
    plt.axvline(x=datetime(1995, 2, 16), color="k", ls="--")
    plt.axvline(x=datetime(1995, 6, 16), color="k", ls="--")
    time_min = datetime(1994, 11, 5)
    myformat = '%.0f'

cs = heat_ax.contourf(ts[:], zs[:], tempdata[:,:], levels, cmap=cmap, extend="both", vmin=tmin, vmax=tmax)
cs.changed()
plt.colorbar(cs, format=myformat)
heat_ax.set_ylim(depthlimit, 0)

heat_ax.set_xlim(xmin=time_min)
plt.show()