from profile import Profile
class SProfile(Profile):
    # apply freshwater flux
    # sign convention: *all* fluxes are positive in the down (towards earth) direction
    # evaporation is positive down, i.e. usually negative
    # precipitation is positive down, i.e. usually positive
    # fw flux is ... = evap + precip (calculated prior to running model and input as a single forcing)
    # fw > 0 -> *gain* of freshwater -> *decrease* in salinity
    # salinity[0] -> salinity[0] * (1-fw/dz) (no dt term bc w is accumulated over 1 hr, like other fluxes)
    def apply_freshwater(self, fw):
        dz = self.dz
        self[0] *= (1-fw/dz)