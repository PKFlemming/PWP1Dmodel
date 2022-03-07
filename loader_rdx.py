import numpy as np
import matplotlib.pyplot as plt
from fetcher import getProfile
import csv
from os.path import join
from parameters import param_dict
from datetime import datetime, timedelta
from time_splitters import get_local_time, get_local_timediff
from matplotlib.patches import Rectangle
import gsw
import matplotlib.dates as mdates

"""
For Github version:
neaten this whole section
separate out function defs where possible
separate off dicts
centralise lat lon data, or get to import from param_dict
"""

study_area = "AS"
# study_area = "AS", "BF", "COARE"
# specify the location and names of the data files
COAREsource_file_dict = {"dir": r"D:\UKuni\3rdYr\MEP\cw\data\COARE",
                        "SFileName": "coaresal",
                        "TFileName": "coaretemp",
                        "velFileName": "coarevel"}
ASsource_file_dict = {"dir": r"D:\UKuni\3rdYr\MEP\cw\data\ArabianSea",
                        "SFileName": "arabianseasal",
                        "TFileName": "arabianseatemp",
                        "velFileName": "arabianseavel"}

# to get plots to display in separate window:
# settings > tools > python scientific > show plot tools in window = disabled
ASminmaxdict = {"tmin" : 27,        # moved up for Summer; for full yr should be 25
                "tmax" : 28.5,
                "salmin" : 35.5,
                "salmax" : 36.5,
                "currmin" : -1,
                "currmax" : 1,
                "densmin" : -2,
                "densmax" : 1}

COAREminmaxdict = {"tmin" : 28.5,
                   "tmax" : 30.5,
                   "salmin" : 34,
                   "salmax" : 35,
                   "currmin" : -2,
                   "currmax" : 2,
                   "densmin" : -4,
                   "densmax" : -2}

bf_minmaxdict = {"tmin" : -0.9,
                 "tmax" : 0.6,
                 "salmin" : 25,
                 "salmax" : 27,
                 "currmax" : 0.5,
                 "currmin" : -0.5,
                 "densmin" : 0.5,
                 "densmax" : 1}

if study_area == "COARE":
    source_file_dict = COAREsource_file_dict
    minmaxdict = COAREminmaxdict
    lon = 156
if study_area == "AS":
    source_file_dict = ASsource_file_dict
    minmaxdict = ASminmaxdict
    lon = 63.5
if study_area == "BF":
    minmaxdict = bf_minmaxdict
    lon = 145 # TODO: correct this


# dir  = source_file_dict["dir"]
# outfilestem = join(dir,"sourceprofiles",study_area)

tmin = minmaxdict["tmin"]                   # hardcoded plot boundaries for each dateset
tmax = minmaxdict["tmax"]
salmin = minmaxdict["salmin"]
salmax = minmaxdict["salmax"]
currmin = minmaxdict["currmin"]
currmax = minmaxdict["currmax"]
densmin = minmaxdict["densmin"]
densmax = minmaxdict["densmax"]
salnums = np.arange(salmin, salmax+0.001, 0.2)
currnums = np.arange(currmin, currmax+0.001, 1, dtype="i")
tnums = np.arange(tmin, tmax+0.001, 0.5, dtype="f")

fig = plt.figure(figsize=(18, 6))

# obs_profile_source is a dict specifying the location of the netcdf file holding the observed profile data
# it's set when this is called in main, if needed
def plotcolumn(col, timestep, forcings, interval, depthlimit, ts_source, obs_profile_source=None, show_initial=False, show_obs=False):
    set_of_forcings = forcings.getset(timestep)
    n_timesteps = col.n_timesteps
    # if we're getting the time stamp based on the initial time specified in param_dict:
    initial_time = datetime.strptime(param_dict["initial_time"], '%d/%m/%Y %H:%M')
    timenow = initial_time + timedelta(seconds = timestep*col.dt)
    timestamp = get_local_time(timenow, lon) # param_dict["lon"]) # TODO: set true local time
    # if we're getting the timestamps from the forcings csv:
    if ts_source == "csv":
        local_timediff = get_local_timediff(lon)
        timestamp_UTC = set_of_forcings["tsval"]
        timestamp_local = timestamp_UTC + local_timediff
    zs = col.zs[:depthlimit]

    heat_ax = plt.subplot2grid((1, 6), (0, 0), colspan=3)
    cmap = plt.cm.rainbow
    levels = np.linspace(tmin, tmax, 20)
    cs = heat_ax.contourf(col.trans_tprofiles, levels, cmap=cmap, extend="both", vmin=tmin, vmax=tmax)
    cs.changed()
    plt.colorbar(cs, format='%.1f')
    heat_ax.set_xlim(0, n_timesteps)
    heat_ax.set_ylim(depthlimit, 0)

    heat_flux_ax = heat_ax.twinx() # TODO: make these plot just once and stay up
    timearray = np.arange(n_timesteps)
    solq_series = forcings.getseries("solq")
    sltq_series = forcings.getseries("sltq")
    solq_thru_time, = heat_flux_ax.plot(timearray, solq_series, label="solq")
    sltq_thru_time, = heat_flux_ax.plot(timearray, sltq_series, label="sltq")
    netq = np.array(solq_series) + np.array(sltq_series)
    netq_thru_time, = heat_flux_ax.plot(timearray, netq, label="netq")
    fluxes = [solq_thru_time, sltq_thru_time, netq_thru_time]
    heat_flux_ax.legend(fluxes, [flux.get_label() for flux in fluxes], loc="lower right")
    ts_series = np.array(forcings.getseries("tsvals")) + local_timediff
    # this should definitely be a separate function. Finds the first midnight in the timeseries
    # and starts xaxis tick locations as near as possible to it
    firstmidnightidx = next(idx for idx, ts in enumerate(ts_series) if
                         ts.hour == 23 and ts.minute > 30 or
                         ts.hour == 0 and ts.minute < 31)
    # add 1 hr because time is within half hour of midnight. This ensures that the date shown is of the day beginning
    # not sure how it'll handle e.g. 0020 tho... will probably display Y m d 01, which would be unhelpful
    tsf_series = [(ts+timedelta(hours = 1)).strftime("%Y/%m/%d %Hh") for ts in ts_series]
    allticklocs = np.arange(firstmidnightidx, n_timesteps, 12)
    n_ticklocs = len(allticklocs)
    ticklocspacing = int(n_ticklocs/10)
    ticklocs = allticklocs[0:n_ticklocs:ticklocspacing]
    ticklabels = np.array(tsf_series)[ticklocs]
    # plt.sca(heat_flux_ax)
    plt.xticks(ticklocs, ticklabels)
    fig.autofmt_xdate() # rotate date labels on x axis

    sal_ax = plt.subplot2grid((1, 6), (0, 3), colspan=1)
    sal_ax.set_ylim((depthlimit, 0))
    sal_ax.set_xlim((salmin,salmax))
    salprof, = sal_ax.plot(col.sprofile[:depthlimit], zs, label="salinity", color="blue")  # plot s profile
    sal_ax.set_xticks(salnums)
    sal_ax.tick_params(axis="x", colors="blue")
    sal_ax.set_xlabel("Salinity")
    sal_ax.xaxis.label.set_color("blue")
    # sal_ax.set_title("Temperature and salinity")
    sal_ax.grid(True)

    temp_ax = sal_ax.twiny()  # plot temp on same subplot as sal, with same y axis
    temp_ax.set_xlim((tmin, tmax))
    temp_ax.set_xticks(tnums)
    temp_ax.tick_params(axis="x", colors="red")
    temp_ax.set_xlabel("Temperature, $\degree$C")
    temp_ax.xaxis.label.set_color("red")
    tprof, = temp_ax.plot(col.tprofile[:depthlimit], zs, label="temperature", color="red")  # plot t profile
    # add dummy blank rectangle in legend with label that tells us what the fw flux is at any given time
    # helps understand why the salinity profile changes the way it does
    # positive fw flux usually accompanied by decrease in surface salinity
    fw_flux_mm = round(set_of_forcings["fw"]*1000,1)
    fw_txt = f"fw_flux = {fw_flux_mm}mm/hr"
    extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label=fw_txt) # blank rectangle
    # add legend
    tsalprofs = [tprof, salprof, extra]
    temp_ax.legend(tsalprofs, [prof.get_label() for prof in tsalprofs], loc="lower right")

    dens_ax = plt.subplot2grid((1, 6), (0, 4), colspan=1)
    dens_ax.set_ylim((depthlimit, 0))
    dens_ax.set_xlim((densmin, densmax))
    dvals = col.dprofile[:depthlimit] - param_dict["reference_density"]
    dens_ax.plot(dvals, zs, label="density anomaly from 1020", color="k")  # plot d profile
    dens_ax.grid(True)
    dens_ax.set_xlabel(r"Density anomaly, kg m$^{-3}$")

    current_ax = plt.subplot2grid((1, 6), (0, 5), colspan=1)
    current_ax.plot(col.uprofile[:depthlimit], zs, label="zonal", color="green")   # plot currents
    current_ax.plot(col.vprofile[:depthlimit], zs, label="meridional", color="purple")
    current_ax.set_xticks(currnums)
    current_ax.set_ylim([depthlimit, 0])
    current_ax.set_xlim([currmin, currmax])
    current_ax.legend(loc="lower left")
    current_ax.grid(True)
    current_ax.set_xlabel("Current velocity, m s$^{-1}$")

    # add horizontal line at mixed layer base
    for ax in [temp_ax, dens_ax, current_ax]:
        ax.axhline(y=col.ml_base_idx*col.dz, color="darkgray", linestyle="--")

    if show_initial is True:                                    # plot initial profiles
        initial_svals = col.initial_sprofile[:depthlimit]
        sal_ax.plot(initial_svals, zs, label="initial salinity", color="lightblue", zorder=1)
        initial_dvals = col.initial_dprofile[:depthlimit] - param_dict["reference_density"]
        dens_ax.plot(initial_dvals, zs, label="initial density anomaly", color="darkgrey", zorder=1)
        initial_tvals = col.initial_tprofile[:depthlimit]
        temp_ax.plot(initial_tvals, zs, label="initial temperature", color="lightpink", zorder=1)

    if show_obs is True:
        source_file_dict = obs_profile_source
        dt_model = col.dt
        dt_observations = source_file_dict["dt_observations"]
        profile_idx = col.initial_profile_idx + timestep * dt_model / dt_observations
        print(f"observed profile index = {int(profile_idx)}")
        try:               # wrapped in simple try/except to handle missing profiles
            observed_svals = getProfile("sal", source_file_dict, record_number=profile_idx)[1]  # get observed salinity profile
            observed_tvals = getProfile("temp", source_file_dict, record_number=profile_idx)[1]  # get temp profile
            observed_uvals = getProfile("east", source_file_dict, record_number=profile_idx)[1]  # get u current profile
            observed_vvals = getProfile("north", source_file_dict, record_number=profile_idx)[1]  # get u current profile
            shortest_profile_depth = min(len(observed_svals), len(observed_tvals), len(observed_uvals))
            profile_limit = min(shortest_profile_depth, depthlimit) # in case shortest profile depth < depth limit
            for observed_val_set in [observed_svals, observed_tvals, observed_uvals, observed_vvals]:
                del observed_val_set[profile_limit:]
            ps = gsw.p_from_z(-1 * col.zs, col.lat)[:profile_limit]
            observed_dvals = gsw.density.rho_t_exact(observed_svals, observed_tvals, ps) - param_dict["reference_density"]

            sal_ax.plot(observed_svals, zs, label="observed salinity", color="blue", zorder=1, ls="--")
            temp_ax.plot(observed_tvals, zs, label="observed temperature", color="red", zorder=1, ls="--")
            dens_ax.plot(observed_dvals, zs, label="observed density", color="k", zorder=1, ls="--")

            obs_current_ax = current_ax.twiny()
            obs_current_ax.plot(observed_uvals, zs, label="observed zonal velocity", color="green", zorder=1, ls="--")
            obs_current_ax.plot(observed_vvals, zs, label="observed merid. velocity", color="purple", zorder=1, ls="--")
            obs_current_ax.set_xlim([-60, 60])
        except:
            print("some observed profile data missing")
    plt.suptitle(f"{study_area} Profiles for {timestamp_local} ({timestamp_UTC} UTC)")
    plt.show(block=False)
    plt.pause(interval)
    plt.cla()

def plotdailies():
    pass