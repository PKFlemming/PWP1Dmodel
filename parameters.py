param_dict = {"cp_sw" : 4005,        # from Sharqawy et al. 2010, setting T and S ~= mean values from Arabian Sea dataset
              "lat" : 15.5,          # 15.5 for Arabian Sea, 74 for bf, -1.75 for COARE TODO: put all site deets in one dict- source file info and lat, lon
              "lon" : 61.5,          # 61.5 for Arabian Sea, 156 for COARE
              "d_air" : 1.2,         # density of air in kg m^-3. d_air and c_drag are used to convert the 10m wind vels
                                     # to wind stresses
              "c_drag" : 10**-3,      # drag coefficient; dimensionless
              "ml_ddbydz_threshold" : 0.1, # in fortran this is 10**-4... sigma units tho
                                          # have chosen this value to target the visible "elbow" on the density profiles
                                          # mention in parameterisation report?
              "ml_dd_from_surface_threshold" : 0.1, # find actual value for this
              "ml_dt_from_surface_threshold" : 0.2, # https://www.nodc.noaa.gov/OC5/WOA94/mix.html
              "ml_dsigma_from_surface_threshold" : 0.125, # get auth. value for this?
              "dt" : 1 * 3600,              # 3 for bf; switch back to 1 for AS
              "initial_time" : "19/10/1994 18:00",
              "reference_density" : 1024, # from fortran
              "g" : 9.81,
              "rG" : 0.3
              }

# this dict sets the file sources and parameters for the profiles, if these are to be built from netcdf
# it also includes the path to the forcings csv
as_source_file_dict = {"sourcedir": r"D:\UKuni\3rdYr\S1\MEP\cw\data\ArabianSea",
                       "sal_filename" : "arabianseasal",
                       "temp_filename" : "arabianseatemp",
                       "vel_filename" : "arabianseavel",
                       "file_suffix" : ".epic",
                       "dt_observations" : 900000/1000,    # timestep between observations, in ms, /1000 (to get in s)
                       "forcings_csv" :
                           r"D:\UKuni\3rdYr\S1\MEP\cw\data\ArabianSea\forcings\AS_forcings_94_10_19_1800on_EPIC.csv",
                           # whole series; epic data. Rates, wind_as_vel=False
                       }

coare_source_file_dict = {"sourcedir": r"D:\UKuni\3rdYr\S1\MEP\cw\data\COARE",
                        "sal_filename": "coaresal",
                        "temp_filename": "coaretemp",
                        "vel_filename": "coarevel",
                        "file_suffix" : ".epic",
                        "dt_observations" : 900000/1000,    # timestep between observations, in ms, /1000 (to get in s)
                        "coare_forcings_csv" :
                              r"D:\UKuni\3rdYr\S1\MEP\cw\data\COARE\forcings\COARE_forcings_EPIC_i1.csv",
                              # whole series; epic data. Rates, wind_as_vel=False
                          }


settings_dict = {"winds_as_vel" : False,
                 "show_initial" : True,
                 "show_obs" : True}
