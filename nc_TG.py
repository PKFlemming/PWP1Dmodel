import netCDF4 as nc
from os.path import join
import numpy as np
import datetime
from time_splitters import datetime_to_jday, datetime_to_ms
from fetcher import getProfile
import mixers
from parameters import param_dict
from parameters import coare_source_file_dict as source_file_dict

# TODO: write_nc and write_nc2 are identical except for how they take in time data. Separate that off and rewrite as one

def write_nc(nc_path, timestamps, depths, temp_arr, dt, day_log, notes=""):

    print("writing data to netcdf")

    ncfile = nc.Dataset(nc_path, mode='w', format='NETCDF4_CLASSIC')
    ncfile.notes = notes # write notes, e.g. on mixings applied

    day_av_ssts = day_log.day_av_ssts
    day_av_mlds = day_log.day_av_mlds
    diurnal_mags = day_log.diurnal_cycle_mags
    n_nudges = day_log.nudges
    n_gradmaxes = day_log.gradmaxes

    ntimes = len(timestamps)
    print(f"timestamps = {timestamps}")
    ndepths = len(depths)
    ndays = len(day_av_ssts)

    print(f"ntimes = {ntimes}")

    # create dimensions
    depth_dim = ncfile.createDimension('depth', ndepths) # depth axis
    time_dim = ncfile.createDimension('time', ntimes) # time axis
    time2_dim = ncfile.createDimension('time2', ntimes)  # time axis
    day_av_sst_dim = ncfile.createDimension('day_av_sst', ndays)  # axis for daily-mean SST
    day_av_mld_dim = ncfile.createDimension('day_av_mld', ndays)  # axis for daily-mean MLD
    diurnal_mag_dim = ncfile.createDimension('diurnal_mag', ndays)  # axis for diurnal cycle magnitude
    nudge_dim = ncfile.createDimension('nudges', 1)  # to record how many times R_b was nudged
    gradmax_dim = ncfile.createDimension('gradmaxes', 1)  # to record how many times R_b was nudged
    day_idx_dim = ncfile.createDimension('day_idx', ndays)  # axis to link daily data with time_dim

    # create variables; assign units and descriptions
    depth = ncfile.createVariable('depth', np.float32, ('depth',))
    depth.units = 'metres'
    time = ncfile.createVariable('time', np.int, ('time',))
    time.units = 'Julian day'
    time2 = ncfile.createVariable('time2', np.int, ('time2',))
    time2.units = 'ms since start of day (UTC)'
    day_av_sst = ncfile.createVariable('day_av_sst', np.float32, ('day_av_sst',))
    day_av_mld = ncfile.createVariable('day_av_mld', np.float32, ('day_av_mld',))
    diurnal_mag = ncfile.createVariable('diurnal_mag', np.float32, ('diurnal_mag',))
    diurnal_mag.long_name = "diurnal_temperature_cycle_magnitude_for_specified_day"
    diurnal_mag.units = "deg C"
    nudges = ncfile.createVariable('nudges', np.int, ('nudges',))
    gradmaxes = ncfile.createVariable('gradmaxes', np.int, ('gradmaxes',))
    day_idx = ncfile.createVariable('day_idx', np.int, ('day_idx',))
    day_idx.description = "This links the daily data (SST, MLD, diurnal cycle magnitude) to the time data. If " \
                          "day_idx[i] == n, day_av_mld[i] is the average mixed layer depth for the day whose julian " \
                          "day value is stored at time[n], etc."

    temp = ncfile.createVariable('temp',np.float64,('depth','time')) # note: unlimited dimension is leftmost
    temp.units = 'deg C' # degrees Centigrade

    # this section matches day data (day average SST, say average MLD, etc.) to time, time2 indexes
    # the time2 index is not important, but this aims to assign an index that is near the middle of the day
    # This is not entirely robust- there might be problems if e.g. the time series starts part way through a day
    records_per_day = 24 * 3600 / dt
    first_day_idx = int(records_per_day/2)
    day_idxs = [first_day_idx + n*records_per_day for n in range(ndays)]

    # populate depth, time, time2 and day dimensions
    # time is julian day; time2 is milliseconds since start of day
    # see description above for day_idx
    depth[:] = depths
    times = [datetime_to_jday(timestamp) for timestamp in timestamps]
    time2s = [datetime_to_ms(timestamp) for timestamp in timestamps]
    time[:] = times
    time2[:] = time2s
    day_idx[:] = day_idxs

    day_av_sst[:] = day_av_ssts
    day_av_mld[:] = day_av_mlds
    diurnal_mag[:] = diurnal_mags

    temp[:,:] = temp_arr

    nudges[:] = [n_nudges]
    gradmaxes[:] = [n_gradmaxes]

    print("-- Wrote data, temp.shape is now ", temp.shape)
    print("-- Min/Max values:", np.round(temp[:,:].min(),2), np.round(temp[:,:].max(),2))

    ncfile.close()

def make_ncpath(hostdir, study_area, dates, suffix=""):
    date0, date1 = dates[0], dates[1]
    date0f = date0.replace("/","_")
    date1f = date1.replace("/","_")
    ncname = study_area + "_" + date0f + "__" + date1f + suffix + ".nc"
    ncpath = join(hostdir, ncname)
    return ncpath

# max_depth is max depth of profile you want returned
# can be max depth of original observed data, if you want the whole profile, or less
def get_checked_tprofile(source_file_dict, profile_idx, max_depth):
    # tnb is the truth value of a check conducted on the raw, uninterpolated, observed profile data
    # if there are enough values at the top and bottom of the profile for interpolation to be conducted safely, tnb
    # will return True and we can use the interpolated profile.
    # If it returns False we fill an array with error values and log that. We need to log a profile of some sort even
    # in the absence of reliable data to ensure consistent indexing.
    # time and time2 are the timestamp of the tprofile, where time = julian day and time2 = ms since start of day
    tprofile, tnb, time, time2 = getProfile("temp", source_file_dict, record_number=profile_idx, get_missing =\
        "topnbottom", get_ts = "as_jdms")[1:5]
    tprofile = tprofile[:max_depth]
    if tnb is False:        # if data is available but it is not safe to interpolate, overwrite with an error profile
        tprofile = [np.nan] * max_depth
    return(tprofile, time, time2)

# this is reused from column. Could generalise into separate function that works for both
# here it returns ml_depth, not ml_base_idx, bc tprofile here is with depths, not idxs
def get_ml_depth_nc(tprofile):
    tprofile = np.array(tprofile)
    dt_from_surface = abs(tprofile - tprofile[0])
    try:
        ml_depth = mixers.first_to_exceed(dt_from_surface, param_dict["ml_dt_from_surface_threshold"])
    except ValueError:
        ml_depth = np.nan
    return ml_depth

def make_tprofile_array(nmin, nmax, source_file_dict, max_depth):
    trans_tprofiles = []
    times = []
    time2s = []
    day_mlds = []
    day_av_mlds = []
    day_ssts = []
    day_av_ssts = []
    diurnal_cycle_mags = []
    n = nmin
    i = 0               # useful to have a 0 indexed counter for checking between loops
    while n <= nmax:
        tprofile, time, time2 = get_checked_tprofile(source_file_dict, n, max_depth)
        times.append(time)
        time2s.append(time2)
        mld = get_ml_depth_nc(tprofile)
        sst = tprofile[0]
        wrapped_tprofile = np.array([tprofile])         # make tprofile 2d
        trans_tprofile = np.transpose(wrapped_tprofile) # make vertical
        if n == nmin:                                   # if first loop, no array yet, so array = first profile
            trans_tprofiles = trans_tprofile
        else:                                           # if not, add profile to existing array
            trans_tprofiles = np.append(trans_tprofiles, trans_tprofile, axis=1)
            if time == times[i-1]:                      # if this loop belongs to the same day as the previous loop
                day_mlds.append(mld)                    # append mld to mlds for this day
                day_ssts.append(sst)                    # likewise for sst
            if time != times[i-1] or n == nmax:         # if it does not, this loop is the first of a new day
                                                        # if n == nmax, we're on the last loop
                day_av_mld = np.nanmean(day_mlds)       # in either case, we average the recorded mlds
                day_av_mlds.append(day_av_mld)          # append the average to the list of daily average mlds
                day_mlds = [mld]                        # and replace the list of mlds with a list containing only
                                                        # the first mld for today
                day_av_sst = np.nanmean(day_ssts)       # likewise for sst
                day_av_ssts.append(day_av_sst)          # TODO: this is all being reused from day_log. Streamline.
                max_sst = max(day_ssts)                 # as we now have a full day's data, get max and min SST
                min_sst = min(day_ssts)
                diurnal_cycle_mag = max_sst - min_sst   # max-min = diurnal SST cycle magnitude
                diurnal_cycle_mags.append(diurnal_cycle_mag)
                day_ssts = [sst]                        # clear day's SSTs and start afresh
        n+=1
        i+=1
        print(i)
    return trans_tprofiles, times, time2s, day_av_mlds, day_av_ssts, diurnal_cycle_mags

def write_nc2(nc_path, times, time2s, depths, temp_arr, day_av_ssts, day_av_mlds, diurnal_mags, dt):

    ncfile = nc.Dataset(nc_path, mode='w', format='NETCDF4_CLASSIC')
    ntimes = len(times)
    ndepths = len(depths)
    ndays = len(day_av_ssts)

    # create dimensions
    depth_dim = ncfile.createDimension('depth', ndepths) # depth axis
    time_dim = ncfile.createDimension('time', ntimes) # time axis
    time2_dim = ncfile.createDimension('time2', ntimes)  # time axis
    day_av_sst_dim = ncfile.createDimension('day_av_sst', ndays)  # axis for daily-mean SST
    day_av_mld_dim = ncfile.createDimension('day_av_mld', ndays)  # axis for daily-mean MLD
    diurnal_mag_dim = ncfile.createDimension('diurnal_mag', ndays)  # axis for diurnal cycle magnitude
    day_idx_dim = ncfile.createDimension('day_idx', ndays)  # axis to link daily data with time_dim

    # create variables; assign units and descriptions
    depth = ncfile.createVariable('depth', np.float32, ('depth',))
    time = ncfile.createVariable('time', np.int, ('time',))
    time2 = ncfile.createVariable('time2', np.int, ('time2',))
    day_av_sst = ncfile.createVariable('day_av_sst', np.float32, ('day_av_sst',))
    day_av_mld = ncfile.createVariable('day_av_mld', np.float32, ('day_av_mld',))
    diurnal_mag = ncfile.createVariable('diurnal_mag', np.float32, ('diurnal_mag',))
    day_idx = ncfile.createVariable('day_idx', np.int, ('day_idx',))
    temp = ncfile.createVariable('temp',np.float64,('depth','time')) # note: unlimited dimension is leftmost
    temp.units = 'deg C' # degrees Centigrade

    # this section matches day data (day average SST, say average MLD, etc.) to time, time2 indexes
    # the time2 index is not important, but this aims to assign an index that is near the middle of the day
    # This is not entirely robust- there might be problems if e.g. the time series starts part way through a day
    records_per_day = 24 * 3600 / dt
    first_day_idx = int(records_per_day/2)
    day_idxs = [first_day_idx + n*records_per_day for n in range(ndays)]

    # populate depth, time, time2 and day dimensions
    # time is julian day; time2 is milliseconds since start of day
    # see description above for day_idx
    depth[:] = depths
    time[:] = times
    time2[:] = time2s
    temp[:,:] = temp_arr
    day_idx[:] = day_idxs

    day_av_sst[:] = day_av_ssts
    day_av_mld[:] = day_av_mlds
    diurnal_mag[:] = diurnal_mags

    print("-- Wrote data, temp.shape is now ", temp.shape)

    ncfile.close()

# nmin = 0
# nmax = 12845
# max_depth = 124
# dt = param_dict["dt"]
#
# depths = range(max_depth)
# nc_path = r"D:\UKuni\3rdYr\MEP\cw\nc\remakenc_coare_02_ddfrs.nc"
#
# trans_tprofiles, times, time2s, day_av_mlds, day_av_ssts, diurnal_cycle_mags = \
#     make_tprofile_array(nmin, nmax, source_file_dict, max_depth)
# write_nc2(nc_path, times, time2s, depths, trans_tprofiles, day_av_ssts, day_av_mlds, diurnal_cycle_mags, dt)

