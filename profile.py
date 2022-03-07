import numpy as np
from parameters import param_dict
class Profile:
    def __init__(self, data):
        self.data = np.array(data)
    def __getitem__(self, position):
        return self.data[position]
    def __setitem__(self, position, value):
        self.data[position] = value

    def __repr__(self):
        return str(self.data) # print(profile) -> print data
    def __len__(self):
        return len(self.data)
    def gr_mix_profile(self, j, r_min, r_critical = 0.25):
        rG = param_dict["rG"]
        # pseudo code:
        # mu{var} = (var[j] - var[j+1])/2
        # d{var} = (1- rg / rG) * mu{var}
        # var[j] = var[j] -d{var}
        # var[j+1] = var[j+1] + d{var}
        # mu_var = (self[j] - self[j+1])/2
        # d_var = (1 - r_min/rG) * mu_var
        # self[j] = self[j] - d_var
        # self[j+1] = self[j+1] + d_var

        rcon = 0.02 + (r_critical - r_min) / 2  # rc = critical richardson value (0.25); r = current (minimum) value
        rnew = r_critical + rcon / 5.
        f = 1 - r_min / rnew  # d_var = f * mu_var
        mu_var = (self[j] - self[j + 1]) / 2 # if top one (_j) is bigger, this is positive      smaller -> negative
        d_var = f * mu_var                   # ...so this is positive                                   -> negative
        self[j] = self[j] - d_var            # ...so the top one gets *smaller*                 top gets bigger
        self[j + 1] = self[j + 1] + d_var    # ...and the bottom one gets *bigger*              bottom gets smaller
                                             # ...and the difference between them gets *smaller* (in both cases)