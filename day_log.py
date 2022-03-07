import numpy as np
class Day_log:
    def __init__(self):
        self.day_ssts = []                      # sea surface temperatures for this day
        self.day_mlds = []                      # mixed layer depths for this day
        self.day_av_ssts = []                   # daily-mean ssts for model run
        self.day_av_mlds = []                   # daily-mean mlds for model run
        self.diurnal_cycle_mags = []            # diurnal cycle magnitudes for model run
        # TODO: record nudges and gradmaxes in separate error log
        # also, more generally, get it to record all settings in the netcdf:
        #       mixings applied
        #       MLD method
        #       name of forcings csv
        self.nudges = 0                         # number of times bulk Richardson mixing has been nudged
        self.gradmaxes = 0                      # number of times gradient Richardson mixing has maxed out
    def log_sst_mld(self, sst, mld, time_since_beginning):
        self.day_ssts.append(sst)  # log surface temp to day record
        self.day_mlds.append(mld)  # log mixed layer depth to day record
        if time_since_beginning % (24 * 3600) == 0 and time_since_beginning > 0:  # if hrs since start is multiple of 24
            # process sea surface temperatures
            # find and log day average SST
            day_av_sst = np.mean(self.day_ssts)  # take average of recorded SSTs
            self.day_av_ssts.append(day_av_sst)  # log to record of day-average SSTs
            # find and log magnitude of diurnal cycle
            max_sst = max(self.day_ssts)
            min_sst = min(self.day_ssts)
            diurnal_cycle_mag = max_sst - min_sst
            self.diurnal_cycle_mags.append(diurnal_cycle_mag)
            # clear day record of SSTs so we can start again for the next day
            self.day_ssts = []
            # find and log day average mixed layer depth
            day_av_mld = np.mean(self.day_mlds)  # take average of recorded MLDs
            self.day_av_mlds.append(day_av_mld)  # log to record of day-average MLDs
            self.day_mlds = []  # clear day record