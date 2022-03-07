import numpy as np
import csv
from math import exp, cos, sin
import gsw  # https://github.com/TEOS-10/GSW-Python
from profile import Profile
from tprofile import TProfile
from sprofile import SProfile
from dprofile import DProfile
from currentprofile import CurrentProfile
import mixers
from parameters import param_dict
import absorption
from fetcher import getProfile, find_matching_jday
from jdcal import gcal2jd
from datetime import datetime, timedelta
import netCDF4 as nc

# TODO: there's a problem with self.zs if observed profile depth changes. Make self.zs dynamic?

class Column:
    def __init__(self, obs_profile_source, dates, obs_format): # TODO: rewrite here to take sourcefile, obs_format. ft = csv or netcdf
        self.mode="run"
        if obs_format == "csv":
            print("running csv section of if statement")
            with open(obs_profile_source) as f:  # read in data
                data = list(csv.reader(f))
            headers = data[0]           # first row of csv is headers Z, S, T, U, V for depth, sal, temp, currents
            data = data[1:]             # remaining rows are data
            zidx = headers.index("Z")
            sidx = headers.index("S")
            tidx = headers.index("T")
            uidx = headers.index("U")
            vidx = headers.index("V")
            self.nlevels = len(data)
            self.make_zero_profiles()          # make profiles with zeros as placeholders
            for level in range(self.nlevels):  # populate profiles with values from source file
                row = data[level]
                self.sprofile[level] = float(row[sidx])  # all wrapped in floats in case data stored as str
                self.tprofile[level] = float(row[tidx])
                self.uprofile[level] = float(row[uidx])
                self.vprofile[level] = float(row[vidx])
                self.zs[level] = float(row[zidx])

        if obs_format == "netcdf":
            source_file_dict = obs_profile_source
            profiledate_greg = dates[0]                                         # gregorian date as string dd/mm/yyyy
            day_greg, month_greg, year_greg = profiledate_greg.split("/")       # get day, month, year values
            profiledate_jd = sum(gcal2jd(year_greg, month_greg, day_greg))+0.5  # convert to Julian Day
            profile_idx = find_matching_jday(source_file_dict, profiledate_jd)  # find first match in jdays of profiles
            self.initial_profile_idx = profile_idx                              # log for later use
            print(self.initial_profile_idx)
            zs, SVals = getProfile("sal", source_file_dict, record_number=profile_idx)[:2]  # get s profile, depth points
            TVals = getProfile("temp", source_file_dict, record_number=profile_idx)[1]      # get temp profile
            UVals = getProfile("east", source_file_dict, record_number=profile_idx)[1]      # get current profiles
            VVals = getProfile("north", source_file_dict, record_number=profile_idx)[1]
            self.nlevels = min(len(SVals), len(TVals), len(UVals))              # find the shortest profile
            for p in (zs, SVals, TVals, UVals, VVals):                          # truncate all profiles to this depth
                del p[self.nlevels:]
            self.zs = np.array([z - 0.5 for z in zs])  # profile fetcher was originally configured to return midpoints
            # mixers.smooth_inversions(SVals)
            # mixers.smooth_inversions(TVals)
            self.sprofile = SProfile(SVals)            # populate s and t profiles
            self.tprofile = TProfile(TVals)
            n_zeros = np.zeros(self.nlevels)
            self.uprofile = CurrentProfile(n_zeros) # populate current and density profiles with 0s
            self.vprofile = CurrentProfile(n_zeros)
            self.dprofile = DProfile(n_zeros)
            # self.uprofile = CurrentProfile(UVals)  # populate current and density profiles with 0s
            # self.vprofile = CurrentProfile(VVals)


        self.tprofiles=[self.tprofile]                  # make a list where we store all the tprofiles for graphing
                                                        # another tprofile will be produced at each model iteration
        self.trans_tprofiles = np.transpose(self.tprofiles)

        self.lat = param_dict["lat"]                    # lat attribute is set here bc needed for density calculations
        self.calculate_dprofile_gsw()

        self.initial_tprofile = self.tprofile.data.copy()       # log initial profiles for graphing
        self.initial_sprofile = self.sprofile.data.copy()
        self.initial_dprofile = self.dprofile.data.copy()

        self.dz = self.zs[1] - self.zs[0]               # extract dz from zs
        for profile in self:                            # set dz for all profiles (needed for e.g. fw flux on sprofile)
            profile.dz = self.dz
        self.dprofile.dz = self.dz                      # density done separately (see comment on __getitem__, below
        self.dt = param_dict["dt"]

    def make_zero_profiles(self):
        n_zeros = np.zeros(self.nlevels)                 # construct empty profiles using relevant class constructors
        self.sprofile = SProfile(n_zeros)
        self.tprofile = TProfile(n_zeros)
        self.uprofile = CurrentProfile(n_zeros)
        self.vprofile = CurrentProfile(n_zeros)
        self.dprofile = DProfile(n_zeros)
        self.sigmaprofile = DProfile(n_zeros)
        self.zs = n_zeros

    # getter method to allow iteration over profiles in column. Doesn't include density because on most occasions when
    # we want to iterate over the profiles we don't want density- for mixing, e.g., we mix u, v, s and t and then
    # calculate density *from* s and t
    def __getitem__(self, position):
        profiles = [self.sprofile, self.tprofile, self.vprofile, self.uprofile]
        return profiles[position]

    # method to calculate densities using the gsw python package (Gibbs Seawater; TEOS 10)
    def calculate_dprofile_gsw(self, loop=0):
        if self.mode == "test":
            print("calculating density profile")
        ps = gsw.p_from_z(-1* self.zs, self.lat)
        sprofile = self.sprofile
        tprofile = self.tprofile
        dprofile = gsw.density.rho_t_exact(sprofile, tprofile, ps)
        self.dprofile[:] = dprofile

    def calculate_sigmaprofile(self):
        ps = gsw.p_from_z(-1 * self.zs, self.lat)
        p_ref = [1000]*len(ps)
        self.sigmaprofile = gsw.pot_rho_t_exact(self.sprofile, self.tprofile, ps, p_ref)

    # apply sensible, latent and terrestrial heat fluxes to top level of column
    def apply_sltq(self, heat_in, cp_sw = param_dict["cp_sw"]):   # specific heat capacity of seawater J kg^-1 K^-1
        if self.mode == "test":
            print("applying sltq")
        self.tprofile[0] += heat_in /(self.dz * self.dprofile[0] * cp_sw)

    def apply_solq(self, incident_solq):
        # I_lw(/sw) is fraction of solar radiation presumed to be longwave(/shortwave)
        if self.mode == "test":
            print("applying solq")
        n_levels = len(self.tprofile)
        for level in range(n_levels):
            dz = self.dz
            z = level*dz
            fraction_absorbed = absorption.fraction_absorbed(z, dz)
            q_in_at_level = incident_solq * fraction_absorbed
            self.tprofile[level] += (q_in_at_level /(dz * self.dprofile[0] * param_dict["cp_sw"]))

    # abortive attempts at damping functions
    def apply_damper(self, damping_factor_in_per_days):
        dt_in_days = self.dt / (24*3600)
        damper = damping_factor_in_per_days * dt_in_days
        self.tprofile.data *= (1-damper)

    def apply_sin_damper(self, n, damping_factor_in_per_days, damping_proportion=0):
        dt_in_days = self.dt / (24*3600)
        damper = sin((n * damping_factor_in_per_days * dt_in_days)%1 * 2 * np.pi)
        self.tprofile.data *= (1 + damper*damping_proportion)

    def apply_linest_damper(self, damping_factor_in_deg_per_day):
        self.tprofile.data -= (damping_factor_in_deg_per_day * self.dt / (24 * 3600))

    # apply freshwater flux
    # sign convention: *all* fluxes are positive in the down (towards earth) direction
    # evaporation is positive down, i.e. usually negative
    # precipitation is positive down, i.e. usually positive
    # fw flux is ... = evap + precip (calculated prior to running model and input as a single forcing)
    # fw > 0 -> *gain* of freshwater -> *decrease* in salinity
    # salinity[0] -> salinity[0] * (1-fw/dz) (no dt term bc w is accumulated over 1 hr, like other fluxes)
    def apply_freshwater(self, fw):
        if self.mode == "test":
            print("applying fw")
        self.sprofile[0] *= (1-fw/self.dz)

    def mix_convective(self, mixings_to_show, cc_depth):
        if self.mode == "test":
            print("checking for convective instability")
        nloops = 0 # counter for convenience during dev; remove later
        dprofile = self.dprofile
        stuvprofiles = [self.sprofile, self.tprofile, self.uprofile, self.vprofile] # could update. Column getter now excludes density
        density_deltasigns = mixers.get_deltasigns(dprofile)
        # mixers.mix_to_n(stuvprofiles, cc_depth) # used in sensitivity testing to conv. mix to constant depth
        while any(density_deltasign < 0 for density_deltasign in density_deltasigns):
            if self.mode == "test":
                print("found instability")
            last_neg_ddelt = mixers.get_last_neg(density_deltasigns)
            if any(val in mixings_to_show for val in ["convective", "all"]):
                print(f"convective mixing to {last_neg_ddelt + 1}")
            mixers.mix_to_n(stuvprofiles, last_neg_ddelt+1)
            self.calculate_dprofile_gsw()
            density_deltasigns = mixers.get_deltasigns(dprofile)
        #     nloops += 1
            # print(nloops)
    # if ml_method == dd_by_dz, ML base index is defined as the first level at which the density gradient dd/dz
    # (d density/d depth) > some threshold value
    # returns index of bottom level of mixed layer, before the density discontinuity
    # need to see how consistently this picks out the ML... could pick some other measure of heterogeneity... 2nd diffs?
    # This gives a (very plausibly) thick ML, and therefore doesn't require much bulk richardson mixing
    # It therefore remains somewhat heterogeneous
    # if ml_method == dd_from_surface, ml base idx is the first level at which (density value minus the surface density
    # value) > threshold. This tends to set idx = 1 (with this threshold value)
    def set_ml_base_idx(self, ml_method):
        if self.mode == "test":
            print("setting ml base index")
        dprofile = self.dprofile
        dds = np.diff(dprofile)  # make array of first differences of density
        ddbydzs = dds / self.dz  # divide by depth step to get change of density with depth between level
                                 # and next level
        if ml_method == "dd_by_dz":
            try:
                ml_base_idx = mixers.first_to_exceed(ddbydzs, param_dict["ml_ddbydz_threshold"])
                # ML base index = level at which density gradient first exceeds threshold...
            except:
                ml_base_idx = np.argmax(ddbydzs)
                # ...unless the gradient never exceeds that threshold. Then idx = level of maximum density gradient.
        elif ml_method == "dd_from_surface":
            dd_from_surface = abs(dprofile - dprofile[0])
            try:
                ml_base_idx = mixers.first_to_exceed(dd_from_surface, param_dict["ml_dd_from_surface_threshold"])
            except ValueError:
                ml_base_idx = np.argmax(dd_from_surface)
        elif ml_method == "dt_from_surface":
            tprofile = self.tprofile
            dt_from_surface = abs(tprofile - tprofile[0])
            ml_base_idx = mixers.first_to_exceed(dt_from_surface, param_dict["ml_dt_from_surface_threshold"])
        elif ml_method == "dsigma_from_surface":
            sigmaprofile = self.sigmaprofile
            dsigma_from_surface = abs(sigmaprofile - sigmaprofile[0])
            ml_base_idx = mixers.first_to_exceed(dsigma_from_surface, param_dict["ml_dsigma_from_surface_threshold"])
        # print(ml_base_idx, ml_base_idx_sigma, np.round(abs(ml_base_idx-ml_base_idx_sigma)/ml_base_idx),2)
        if ml_base_idx < 1:
            ml_base_idx = 1
        self.ml_base_idx = ml_base_idx
        self.ml_thickness = ml_base_idx * self.dz

    def rotate_by_angle(self, angle):
        if self.mode == "test":
            print("rotating")
        u = self.uprofile[:]
        v = self.vprofile[:]
        u2 = u * cos(angle) - v * sin(angle)
        v2 = u * sin(angle) + v * cos(angle)
        self.uprofile[:] = u2
        self.vprofile[:] = v2

    # pseudo code:
    # acceleration = windstress/(ml_thickness * average density)
    # can think of this as (force per unit area)/(mass per unit area); F=ma -> F/m = a
    # change in velocity dvel = acceleration*dt (dt in seconds)
    # vel += dvel for array down to ml_idx
    # NOTE: this gives tiny values for dv (OOM 10**-3), but this appears to be realistic
    # wind stresses in ASea met are same OOM as those calculated from the ERA5 wind velocity data, tho they
    # are *not* the same values (this should be addressed)

    def apply_wind_stresses(self, wind_stresses_u_v):
        if self.mode == "test":
            print("applying winds")
        ml_base_idx = self.ml_base_idx
        ml_thickness = self.ml_thickness
        d_average_in_ml = np.mean(self.dprofile[:ml_base_idx])
        accelerations_u_v = [wind_stress/(ml_thickness * d_average_in_ml) for wind_stress in wind_stresses_u_v]
        [dvel_u, dvel_v] = [acceleration * param_dict["dt"] for acceleration in accelerations_u_v]
        self.uprofile[:ml_base_idx] += dvel_u
        self.vprofile[:ml_base_idx] += dvel_v

    def get_ml_delta(self, profile):
        ml_base_idx = self.ml_base_idx
        value_in_ml = profile[0]  # this is approach taken in fortran
                                  # PROBLEM: this is difference from top to bottom of ML, not difference from base
                                  # to next layer. It's the same iff the ML is homogenous
        value_in_next_level = profile[ml_base_idx]      # was using ml_base_idx+1... try this for testing purps.
        ml_delta = value_in_next_level - value_in_ml
        return ml_delta

    def set_bulk_richardson(self):
        # pseudo code:
        # R_b = (g * delta(density) * mixed_layer_thickness)/(reference density * delta(velocity)^2)
        # delta(var) = var(level below mixed layer) - var(mixed layer)
        ### PROBLEM: assumes homogeneous ML ###
        # fortran uses var(surface)... will try that. Maybe density smoothing will make ml homogenous...?
        if self.mode == "test":
            print("setting bulk Richardson number")
        ml_base_idx = self.ml_base_idx
        ml_thickness = self.ml_thickness
        ddensity = self.get_ml_delta(self.dprofile)
        du = self.get_ml_delta(self.uprofile)
        dv = self.get_ml_delta(self.vprofile)

        reference_density = param_dict["reference_density"]
        # reference_density = self.dprofile[0] # changed this as an experiment... deep bulk mixing is causing probs
        g = param_dict["g"]
        dvel_sq = du**2 + dv**2
        if ml_base_idx > self.nlevels - 2:
            R_b = 0.65
            # when spinning up, high currents can lead to bulk mixing down to base of column. In that case we set
            # R_b to 0.65 and continue. Otherwise it tries to mix off the bottom of the column.
        elif dvel_sq > 0:
            R_b = (g * ddensity * ml_thickness) / (reference_density * dvel_sq)
        else:
            R_b = np.inf
        self.R_b = R_b

    def bulk_r_mix(self, mixings_to_show, day_log, nudgepoint=None, nudgeval=0):
        if self.mode == "test":
            print("bulk richardson mixing")
        R_b_critical = 0.64 # 0.65 canonically; changed to avoid a specific over-mixing event
        i=0
        while self.R_b < R_b_critical:                                   # while(/if) it's less than the critical value...
            if any(val in mixings_to_show for val in ["bulk", "all"]):
                print(f"R_b = {round(self.R_b,3)}; bulk mixing to {int(self.ml_base_idx * self.dz)}m")
            for profile in self:
                mixers.mix_to_n(profile, self.ml_base_idx)   # ... mix stuv profiles down to the base level...
            self.calculate_dprofile_gsw()                    # ...recalculate density...
            self.ml_base_idx += 1                            # ...add another level to the mixed layer to be mixed in next time..
            self.set_bulk_richardson()                       # and recalculate bulk richardson number
            if i == nudgepoint:
                day_log.nudges += 1                         # this logs how many times R_b has been nudged
                R_b_critical = nudgeval
                print("nudging R_b")
                # this is a fudge. The bulk stability criterion occasionally entrains >80 layers in one loop
                # and cools the SST hugely in a way that doesn't square at all with observations
                # In my runs of ~one year of hourly data, it does this 0-1 times per run
                # This fudge relaxes stringency of the bulk stability criterion once we've mixed {nudgepoint} layers.
                # R_b_critical is reset when bulk_r_mix is next called
                # On my data, with nudgepoint = 10, this function was called on about 1 loop in 20, and raised MLD
                # by average 5m. This was a problem, so the nudgepoint was raised to 30
            i += 1

    # function to get gradient richardson number at a given level
    def get_gradient_richardson(self, j):
        # pseudo code:
        # R_g = (g * dd * dz) / (d_ref * (du^2 + dv^2)) where d{var} = var_(j+1) - var_j
        g = param_dict["g"]
        d_ref = param_dict["reference_density"]
        dd = self.dprofile[j+1] - self.dprofile[j]
        dz = self.dz
        du = self.uprofile[j+1] - self.uprofile[j]
        dv = self.vprofile[j+1] - self.vprofile[j]
        if du*dv == 0:
            R_g = np.inf
        else:
            R_g = (g * abs(dd) * dz) / (d_ref * (du**2 + dv**2))
            # abs val of dd has been taken here as a patch. When processing 13/06/95, AS data, some dd came back
            # negative, -> R_g < 0, and it couldn't be resolved
        return R_g

    # set the gradient richardson numbers for each level in the column (except the bottom one)
    def set_grprofile(self):
        if self.mode == "test":
            print("calculating gradient richardson profile")
        nlevels = self.nlevels
        self.grprofile = np.empty(nlevels-1)
        for j in range(nlevels - 1):
            self.grprofile[j] = self.get_gradient_richardson(j)

    # scan gr profile
    # if min value is below critical level, apply gradient richardson mixing
    def gr_mix_col(self, mixings_to_show, day_log):
        if self.mode == "test":
            print("gradient richardson mixing")
        min_gr = min(self.grprofile)
        i = 0
        while min_gr < 0.25 and i < 1000:                       # could put this in param dict for flexibility
            min_gr_index = np.argmin(self.grprofile)            # this level and the one below are mixed
            if any(val in mixings_to_show for val in ["gradient", "all"]):
                print(f"min. R_g = {round(min_gr,3)}; gradient mixing at {int(min_gr_index*self.dz)}m")
            for profile in self:
                profile.gr_mix_profile(min_gr_index, min_gr)
            self.calculate_dprofile_gsw()
            gr_mix_min_idx = max(min_gr_index-1,0)              # these four layers are recalculated
            gr_mix_max_idx = min(min_gr_index+2, self.nlevels-1)
            for j in range(gr_mix_min_idx, gr_mix_max_idx):
                self.grprofile[j] = self.get_gradient_richardson(j)
            min_gr = min(self.grprofile)
            i += 1            # to prevent endless loops.
            if i == 1000:
                day_log.gradmaxes += 1                          # add 1 to record of times gr mixing has maxed out
                print("gradient Richardson mixing maxed out")

    # we log both tprofiles and transposed tprofiles. A tprofile has temps for each depth at a given time:
    # [T_(n=n, z=0), T_(n=n, z=1)... T_(n=n, z=max.depth)]
    # a transposed t profile has all temps generated so far for a particular depth:
    # [T_(n=0,z=z), T_(n=1, z=z)... T_(n=latest timestep, z=z)]
    # transposed tprofiles are used for heatmapping and recording to csv.
    # TODO: do tprofiles need to persist?
    def log_tprofile(self):
        if self.mode == "test":
            print("logging t profiles")
        tprofile_to_append = (self.tprofile.data).copy()                    # make a copy of the profile
        self.tprofiles.append(tprofile_to_append)                           # append to profiles
        trans_tprofile_to_append = np.transpose([tprofile_to_append])       # transpose
        self.trans_tprofiles = np.append(self.trans_tprofiles, trans_tprofile_to_append, axis=1)    # append to array


    # TODO: make profile writer include lat
