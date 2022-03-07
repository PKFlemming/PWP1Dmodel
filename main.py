"""
This is a python instantiation of the Price-Weller-Pinkwell one-dimensional ocean surface model. It takes meteorological
data and models the upper layer of the ocean, including heat diffusion, mixing and so on. The run command, including
guidance on arguments, is at the bottom.
The mathematical model was designed by Jame F. Price, ROBERT A. Weller and Robert Pinkel, as described in their 1986
paper Diurnal Cycling: Observations and Models of the Upper Ocean Response to Diurnal Heating, Cooling, and Wind Mixing.
This paper can be found (behind paywall) at https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JC091iC07p08411
Free versions are also available.
"""

from column import Column
from forcings import Forcings
import gsw
from loader_rdx import plotcolumn
from day_log import Day_log
# import as_source_file_dict or coare_source_file_dict
from parameters import as_source_file_dict as source_file_dict
import nc_TG

# TODO: currently have to change AS/COARE in too many places- parameters, forcings csv in call to run, from parameters
# import, loader_rdx
# also have to set max/min T in loader_rdx minmaxdict for plotting. Best to set wide, observe a run, then adjust

def run_pwp(obs_profile_source, obs_format, forcings_csv, depthlimit, interval=None, plotstep=None, dates=None,
            recordingstyle="rate", wind_as_vel=False, ml_method="dd_by_dz", m_show="all", mode="run",
            show_initial=False, show_obs=False, nctag="", cc_depth=None):

    # step 0) read in profile data, set mode, read in forcing data
    # construct column with initial T, S, U, V profiles taken from netcdf
    col = Column(obs_profile_source, obs_format=obs_format, dates=dates)
    # density is calculated as part of the column object initialisation, and additional properties dt, dz and latitude
    # are set from the param dict

    col.mode = mode # set column mode ("run" or "test")

    # read in forcing data. wind_as_vels=False if wind data are stored as stresses.
    forcings = Forcings(forcings_csv, wind_as_vel=wind_as_vel, dates=dates, recordingstyle=recordingstyle)
    col.n_timesteps = len(forcings)

    if plotstep == None or plotstep == 0: # in these cases we only want to show plots at the end of the run
        plotstep = col.n_timesteps

    day_log = Day_log() # construct a day_log object to record average day SST and MLD

    for n in range(col.n_timesteps):
        set_of_forcings = forcings.getset(n) # get set of forcings for this timestep

        # step 1) apply heat fluxes
        col.apply_sltq(set_of_forcings["sltq"])     # apply sensible, latent and terrestrial heat fluxes
        col.apply_solq(set_of_forcings["solq"])     # apply solar heat flux

        # step 2) apply freshwater flux
        col.sprofile.apply_freshwater(set_of_forcings["fw"]) # apply freshwater forcing

        # step 3) recalculate density profile; run convective mixing
        col.calculate_dprofile_gsw(loop=n)
        if ml_method == "dsigma_from_surface": # make potential density profile if needed
            col.calculate_sigmaprofile()
        col.mix_convective(m_show, cc_depth=cc_depth)

        # step 4) rotate; apply momentum flux to mixed layer
        # first, find and log mixed layer base index
        col.set_ml_base_idx(ml_method = ml_method) # this is now accessible as col.ml_base_idx

        coriolis_param = gsw.f(col.lat)                                         # get rotation parameter
        half_rotation = -coriolis_param * col.dt / 2
        col.rotate_by_angle(half_rotation)                                      # do half of rotation
        col.apply_wind_stresses([set_of_forcings["u"], set_of_forcings["v"]])   # apply wind stress
        col.rotate_by_angle(half_rotation)                                      # do second half of rotation

        # step 5) bulk Richardson mixing
        col.set_bulk_richardson()                                             # calculate bulk richardson number
        col.bulk_r_mix(m_show, day_log, nudgepoint=28, nudgeval=0.5)          # check for and remove any instabilities
        # bulk mixing occasionally goes too far and mixes in e.g. 90 column levels in a row in a way that disagrees
        # with observations. These nudge values allow one to avoid that. See column.py for details.

        # step 6) gradient richardson mixing
        col.set_grprofile()                            # calculate gradient richardson number for every level in column
        col.gr_mix_col(m_show, day_log)                # if any are subcritical, mix till they aren't

        # step 7) readouts and plotting
        pc = round(100 * (n + 1) / col.n_timesteps, 1)                  # get completion %
        print(f"loop {n} ({pc}% complete). ML thickness = {int(col.ml_thickness)}m; bulk Richardson no. = "
              f"{round(col.R_b, 1)}")

        # log info on daily-mean sea surface temp and daily-mean mixed layer depth
        day_log.log_sst_mld(col.tprofile[0], col.ml_thickness, n*col.dt)
        # log the processed tprofile
        col.log_tprofile()

        if (n+1) % plotstep == 0: # if timestep is multiple of plotstep, call plotter function
            print("plotting...")
            plotcolumn(col, n, forcings, interval, depthlimit, ts_source="csv", obs_profile_source=source_file_dict,
                       show_initial=show_initial, show_obs=show_obs)

    # write to netcdf
    host_dir = r"D:\UKuni\3rdYr\S1\MEP\cw\nc\AS\noFW"
    nc_path = nc_TG.make_ncpath(host_dir, "AS", dates, nctag)
    # notes will be written to the netcdf file as an attribute of the dataset
    notes = "no freshwater, ml_method = dd"
    ts_series = forcings.getseries("tsvals")
    nc_TG.write_nc(nc_path, ts_series, col.zs, col.trans_tprofiles[:,1:], col.dt, day_log, notes)

# set wind_as_vel=True if wind forcing data is stored as velocities; False if it's stored as momentum fluxes
# ml_method determines the method used to find the mixed layer base index (and thus thickness, bulk Richardson
# number etc.). If ml_method == "dd_by_dz", these are calculated in terms of density gradient (ddensity/dz).
# If ml_method == "dd_from_surface", they are calculated by density-at-level minus density-at-surface.
# ml_method == "dt_from_surface" -> temperature criterion (0.2 or 0.5 deg C, set in parameters)
# ml_method == "dsigma_from_surface" -> potential density criterion (0.125 kg m^-3)
# m_show (mixings to show) is a list of some combination of "bulk", "convective" and "gradient", or "all", to determine
# what info is read out at each iteration while the model runs.
# Set mode = "test" to see a more detailed readout in the console of what step is being carried out at each moment
# dates can be None (default value) or a list of min,max dates formatted ["dd/mm/yyyy","dd/mm/yyyy"]. These will then be
# the min and max date values (inclusive) taken from the forcings csv.
# plotstep determines how frequently plots are presented. If plotstep == 1, profile plots will be shown at every loop.
# If plotstep == 10, plots will be shown every ten loops. Default is None, which -> plots only shown on conclusion of
# model run
dates = ["4/11/1994", "11/11/1994"]
m_show = [""]
forcings_csv = source_file_dict["forcings_csv"]
run_pwp(obs_profile_source=source_file_dict, obs_format="netcdf", forcings_csv=forcings_csv, depthlimit=100,
        interval=0.01, plotstep=4, dates=dates, wind_as_vel=False, ml_method="dt_from_surface",
        m_show=m_show, mode="run", recordingstyle="rate", show_initial=True, show_obs=True, nctag="DT02_CTL_chk4")