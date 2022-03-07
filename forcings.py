import csv
import numpy as np
from parameters import param_dict
from datetime import datetime, timedelta
class Forcings:
    # set wind_as_vel to false if wind data stored as stresses
    # recordingstyle defines whether forcings are rates (heat fluxes as W m^-2, evap and pptn as m s^-1) or amounts
    # accumulated over the recording period (J m^-2, m). The latter is the case for ERA5 data.
    # dates is an optional list of two dates- a min and a max- to get a subsection of longer series.
    def __init__(self, csvfile, wind_as_vel=True, dates=None, recordingstyle="rate"):
        with open(csvfile) as f:  # read in data
            data = list(csv.reader(f))      # all values are positive downward
        headers = data[0]                   # first row of csv is headers
        tsidx = headers.index("timeStamp")  # timestamp
        sqidx = headers.index("sensq")      # sensible heat, J m^-2 (not W- this is 1hr accumulated energy)
        lqidx = headers.index("latq")       # latent (all energy fluxes are same units)
        solqidx = headers.index("solq")     # solar radiation
        terrqidx = headers.index("terrq")   # terrestrial radiation
        uidx = headers.index("u")           # zonal wind velocities, m s^-1
        vidx = headers.index("v")           # meridional
        fwidx = headers.index("fw")         # freshwater flux (evap + precip, already calculated externally)

        # we store the forcing data in two ways:
        # 1) as series, each of which is a *full* timeseries for *one* forcing. These are used for graphing
        # Currently this is only done for the heat fluxes
        self.series_of_forcings = {"solq" : [],
                                   "sltq" : [],
                                   "ts" : [],
                                   "tsvals" : [],}
        # 2) As sets. Each set is a set of *all* forcings (heat fluxes, water fluxes etc.) for *one* timestep
        # These are used for iterating the model: at each timestep, we want one set of forcings. This makes each
        # forcing easily accessible
        self.sets_of_forcings = []

        if dates != None:
            dates = [datetime.strptime(date, '%d/%m/%Y') for date in dates]
            datemin, datemax = min(dates), max(dates)
            datemax += timedelta(hours = 24)

        for row in data[1:]:
            ts = row[tsidx]                                 # get timestamp
            tsval = datetime.strptime(ts, '%d/%m/%Y %H:%M') # convert to datetime object
                                                            # this could work on the jday and ms columns
            if dates != None and tsval > datemax:    # if dates have been set and the datetime object exceeds datemax...
                break                                # ...stop looking
            if dates == None or datemin < tsval < datemax:    # if date limits haven't been set, or we're between them
                sltq = float(row[sqidx])+ float(row[lqidx]) + float(row[terrqidx]) # sensible and latent combined
                solq = float(row[solqidx])
                u = float(row[uidx])            # if wind forcing data is stored as stresses they can be read in directly
                v = float(row[vidx])
                if wind_as_vel == True:         # if stored as velocities, need to be converted to stresses
                    u = vel_to_stress(u)        # my data are stored as velocities, so this is the default
                    v = vel_to_stress(v)
                fw = float(row[fwidx])
                set_of_forcings = {"ts" : ts,           # string
                                   "tsval" : tsval,     # datetime object
                                   "sltq" : sltq,
                                   "solq" : solq,
                                   "u" : u,
                                   "v" : v,
                                   "fw" : fw
                                   }
                for key in ["sltq", "solq", "u", "v", "fw"]:             # catch error values
                    if set_of_forcings[key] > 10**6:
                        set_of_forcings[key] = 0
                        ts = set_of_forcings["ts"]
                        print(f"warning: no {key} forcing data for timestep {ts}. {key} set to 0 for this timestep")

                # if fluxes are rates (W rather than j, ms^-1 rather than m accumulated in time period), multiply by dt
                if recordingstyle == "rate":
                    for key in ["sltq", "solq", "fw"]:
                        set_of_forcings[key] = set_of_forcings[key]*param_dict["dt"]
                # add set of forcing to dict of sets
                self.sets_of_forcings.append(set_of_forcings)

                # add forcings to series of forcings, in the case of those that will be used for longitudinal graphing
                self.series_of_forcings["solq"].append(solq)
                self.series_of_forcings["sltq"].append(sltq)
                self.series_of_forcings["ts"].append(ts)
                self.series_of_forcings["tsvals"].append(tsval)

    def getset(self, position):
        return self.sets_of_forcings[position]
    def getseries(self, key):
        return self.series_of_forcings[key]
    def __len__(self):
        return len(self.sets_of_forcings)

# d_air is density of air in kg m^-3; c_drag is dimensionless drag coefficient
def vel_to_stress(vel, d_air=param_dict["d_air"], c_drag=param_dict["c_drag"]):
    stress = d_air * c_drag * abs(vel) * vel
    return stress