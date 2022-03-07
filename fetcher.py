from os.path import join
import netCDF4 as nc
import numpy as np
from math import ceil, floor
from scipy.interpolate import interp1d, PchipInterpolator
from time_splitters import JDayToGreg, getLocalTS

# function to check whether input is a number
def isnumber(input):
    truthValue = False # by default assume it isn't
    intype = type(input) # get input type
    if intype == float or intype == int or isinstance(input, np.floating):
        truthValue = True # if it's a numeric type, set to true
    return(truthValue)

# takes a pair of lists (e.g, depths and salinities at those depths)
# removes pair of values if either value is not a number
# This is an issue with the Arabian Sea data, which is quite patchy and has lots of nan values
def cleanPairedLists(list1, list2):
    n = 0
    while n < len(list1):
        val1 = list1[n] # iterate over values in lists
        val2 = list2[n]
        valXval = (val1*val2) # get product of values. If either is nan, -> nan.
        if not isnumber(valXval): # if either is not a number, delete entry
            del list1[n:n + 1]
            del list2[n:n + 1]
            n -= 1 # if a pair has been deleted, the lists are shorter, so we need to decrement n
        n += 1 # increment n

# takes a pair of lists
# first list is depths
# second list is values at those depths
# deletes entries where the depth is repeated from earlier in the list (leaves first entry)
# this is an issue with the TOGA COARE data, which has a couple of duplicate measurements. The values are identical,
# so it's not repeat measurements which might be worth averaging
def delRepeatedZs(zs, vals):
    n = 0
    while n < len(zs):
        zsSoFar = zs[:n]
        if zs[n] in zsSoFar:
            del zs[n:n + 1]
            del vals[n:n + 1]
            n -= 1
        n+= 1

# as above, but for unrealistically large z values
def delErrorZs(zs, vals):
    n = 0
    while n < len(zs):
        z = zs[n]
        val = vals[n]
        if z > 10000 or val > 10000:  # there are some error values/ unrealistic values in the COARE ds
            del zs[n:n + 1]
            del vals[n:n + 1]
            n -= 1  # if a pair has been deleted, the lists are shorter, so we need to decrement n
        n += 1  # increment n


# takes two lists
# optional: interval, technique specification
# first list is of depths
# second list is of vals (e.g. temperatures, salinities etc.) at those depths
# default technique is PCHIP (Piecewise Cubic Hermite Interpolating Polynomial), which avoids overshoots
# alternative is cubic
def interpolate(zs, vals, dz=1, technique="PCHIP"):
    try:
        if min(zs) > 10:
            print(f"warning! No data in top {floor(min(zs))} metres!")
    except:
        print("warning! no observed profile data")
    try:
        if technique == "PCHIP":
            interpolator = PchipInterpolator(zs, vals) # make interpolator function
            intdMin = dz/2 # set interpolated min to middle of surface level
            intdMax = ceil(max(zs))-dz/2
            intdRange = intdMax - intdMin
            nPoints = int(intdRange/dz) +1 # number of points = range/interval, +1 bc ends
            intdZs = np.linspace(intdMin, intdMax, nPoints) # array of n points between min and max depths
            intdVals = interpolator(intdZs, extrapolate=True) # use interpolator to model missing vals
        if technique == "cubic":
            interpolator = interp1d(zs, vals, kind='cubic')
            intdMin = ceil(min(zs)) # can't extrapolate with interp1d, so min = min. val (would have to be done separately)
            intdMax = floor(max(zs))
            intdRange = intdMax - intdMin
            nPoints = int(intdRange / dz)
            intdZs = np.linspace(intdMin, intdMax, nPoints + 1)
            intdVals = interpolator(intdZs)
    except ValueError:                              # if there's no observed data, the interpolator throws value error
        intdZs, intdVals = [np.nan]*len(zs), [np.nan]*len(zs)
    return(intdZs, intdVals)


# get profile data from source files, and interpolate
# if get_missing == "all", also returns total number of missing values
# if get_missing == "topnbottom", returns True if there are enough values at the top and bottom for interpolation to
# proceed with confidence and False if not
def getProfile(variable, source_file_dict, lon=0, record_number=0, interpolation="PCHIP", dz=1, get_missing = False,
               get_ts = False):
    # unpack info about location of data, coming from main
    sourcedir = source_file_dict["sourcedir"]                 # dir where datafiles are stored
    sal_filename = source_file_dict["sal_filename"]     # filenames for salinity...
    temp_filename = source_file_dict["temp_filename"]     # ...temperature...
    vel_filename = source_file_dict["vel_filename"] # ...and current velocity
    file_suffix = source_file_dict["file_suffix"]   # file suffix (e.g. .netcdf, .epic...)

    sal_datafile = join(sourcedir, sal_filename + file_suffix)    # data file
    sal_ds = nc.Dataset(sal_datafile)                 # dataset
    TDataFile = join(sourcedir, temp_filename + file_suffix)
    temp_ds = nc.Dataset(TDataFile)
    velDataFile = join(sourcedir, vel_filename + file_suffix)
    vel_ds = nc.Dataset(velDataFile)

    # pair variable name with relevant dataset. This allows getProfile to be called by variable name
    variableDict = {"sal": sal_ds,
                    "temp": temp_ds,
                    "east": vel_ds,
                    "north": vel_ds}

    ds = variableDict[variable]                             # get dataset
    zs = list(ds.variables["depth"][:])                     # get depths
    val_nests = list(ds.variables[variable][record_number])  # get variable values
    vals = [val_nest[0][0] for val_nest in val_nests]       # unpack variable values (the ds store them as [[var]])

    originalLen = len(zs)                                   # TODO: doesn't a float wrapper also unpack these vals?
                                                            # or casting to array? (see get_times, below)

    if get_missing == "topnbottom":                # check original, uncleaned profile for values missing at top
        missing_from_top = check_topnbottom(vals)

    cleanPairedLists(zs, vals)                          # delete non-data entries from lists
    delRepeatedZs(zs, vals)                             # delete duplicate depth entries
    delErrorZs(zs, vals)
    nMissingVals = originalLen - len(zs)                # note number of non-data entries
    if interpolation != None:
        zs, vals = interpolate(zs, vals, dz, interpolation)
    infoToReturn = [list(zs), list(vals)]                           # return depths and values

    if get_missing == "all":
        infoToReturn.append(nMissingVals)               # if asked, also return number of missing values
    elif get_missing == "topnbottom":
        infoToReturn.append(missing_from_top)


    if get_ts == True:
        JDayArray = ds["time"]
        JDay = int(JDayArray[record_number])        # Julian day of record
        msArray = ds["time2"]
        ms = int(msArray[record_number])            # milliseconds since start of day
        try:
            timestamp = JDayToGreg(JDay, ms)        # function to convert JDay and ms to gregorian timestamp
            timestamp = getLocalTS(timestamp, lon)
        except:
            timestamp = "NoTimestamp"
        infoToReturn.append(timestamp)

    if get_ts == "as_jdms":
        jday_array = ds["time"]
        jday = int(jday_array[record_number])        # Julian day of record
        ms_array = ds["time2"]
        ms = int(ms_array[record_number])            # milliseconds since start of day
        try:
            infoToReturn.append(jday)
            infoToReturn.append(ms)
        except:
            pass

    return infoToReturn

def find_matching_jday(source_file_dict, jday_to_match):
    sourcedir = source_file_dict["sourcedir"]           # dir where datafiles are stored
    sal_filename = source_file_dict["sal_filename"]     # filename for sal (files are assumed to have the same times)
    file_suffix = source_file_dict["file_suffix"]       # file suffix (e.g. .netcdf, .epic...)
    sal_datafile = join(sourcedir, sal_filename + file_suffix)  # data file
    sal_ds = nc.Dataset(sal_datafile)                   # dataset
    jday_array = np.array(sal_ds["time"])               # make array of times
    instances = np.where(jday_array == jday_to_match)[0]# make array of values that match target day. This comes back
                                                        # nested; first member of nested array is array of instances
    first_instance = instances[0]                       # take first instance
    return first_instance                               # this is index of first profile on target day

def first_present_value(val_array):
    idx = min(idx for idx, val in enumerate(val_array) if val > 0)
    return idx

# this function checks the top and bottom of the source profile.
# It returns False if either the top or the bottom of the profile fails to meet the checks
# Because PCHIP interpolation is used, if the top and bottom are fixed all values in between will at least be reasonable
def check_topnbottom(val_array):
    truth_value = False # False by default
    # switch to True if we have the top value or the second *and* third from top
    if isnumber(val_array[0]) or (isnumber(val_array[1]) and isnumber(val_array[2])):
        truth_value = True
        if isnumber(val_array[-1]) is False: # if we don't have the bottom value, switch back to False
            truth_value = False

    return truth_value