from os.path import join
import csv
import numpy as np
from statistics import mean
from parameters import param_dict
import matplotlib.pyplot as plt

# takes an array and an index (0-indexed)
# sets all members of the array up to *and including* array[n] to the mean of said members
def mix_to_n_single(array_in, n):
    array_in[:n+1] = mean(array_in[:n+1])

# takes an array or a list of arrays (*not* a list of vals)                             used
# functions as above, but can be applied to one array *or* to several at once
def mix_to_n(array_or_arrays_in, n):
    if isinstance(array_or_arrays_in, list):
        arrays_in = array_or_arrays_in
        for array_in in arrays_in:
            mix_to_n_single(array_in, n)
    else:
        array_in = array_or_arrays_in
        mix_to_n_single(array_in, n)

# takes an array                                                                          used
# takes the first differences of the members of that array
# takes the signs of those differences. Differences are forward-looking, i.e. difference[0] = array[1]-array[0]
def get_deltasigns(array_in):
    deltas = np.diff(array_in)
    deltasigns = np.sign(deltas)
    return deltasigns

# this gets the last negative number in an array. If there is no negative number, it throws an error
# need to check that there is a neg. number *before* using this function
def get_last_neg(array_in):                                                             # used
    if min(array_in) < 0:
        last_neg_idx = max(idx for idx, val_in in enumerate(array_in) if val_in < 0)
    return last_neg_idx

def get_last_neg_delt(array_in):                                                          # not used
    deltasigns = get_deltasigns(array_in)
    if any(deltasign < 0 for deltasign in deltasigns):
        last_neg_delt_idx = get_last_neg(deltasigns)
    return last_neg_delt_idx

# takes an array                                                                            not used
# finds the last pair that have a negative first difference, i.e. array[n+1] < array[n]
# mixes (averages) array members up to *and including* the second member of that pair
def mix_to_last_neg_single(array_in):
    last_neg_delt_idx = get_last_neg_delt(array_in)
    array_in = mix_to_n(array_in, last_neg_delt_idx+1)

# takes an array or a list of arrays                                                       # not used
# functions as above, but can be applied to one array *or* to several at once
def mix_to_last_neg(array_or_arrays_in):
    if isinstance(array_or_arrays_in, list):
        arrays_in = array_or_arrays_in
        for array_in in arrays_in:
            mix_to_last_neg_single(array_in)
    else:
        array_in = array_or_arrays_in
        mix_to_last_neg_single(array_in)

# takes an array                                                                           ## not used
# finds the last negative delta; mixes up to and including
# re-checks for negative deltas; keeps mixing until array is monotonic increasing
def mix_to_monotonicity(array_in):
    deltasigns = get_deltasigns(array_in)
    while any(deltasign < 0 for deltasign in deltasigns):
        mix_to_last_neg(array_in)
        deltasigns = get_deltasigns(array_in)

# takes an array and a threshold value                                  # used
# finds first value in array to exceed threshold value
# this is to find the base of the mixed layer, defined as the point at which ddensity/dz > threshold
def first_to_exceed(val_array, threshold):
    idx = min(idx for idx, val in enumerate(val_array) if val > threshold)
    return idx

# these next two aren't used yet
# their function is to remove any inversions by mixing adjacent levels- may be helpful for prep.ing initial profiles
# rather than mixing down from surface, which smooths all the way down to an inversion, this applies local smoothing
# at the site of the inversion. This produces a monotonic profile. In practice this pair of functions would need to be
# split further bc the smoothing is applied to the t,sprofiles and the inversion-checking to the dprofile. Just logging
# the idea here for possible future revisiting.
def get_neg_idxs(array_in):
    if min(array_in) < 0:
        neg_idxs = [idx for idx, val_in in enumerate(array_in) if val_in < 0]
    return neg_idxs

def smooth_inversions(array_in, direction):
    delts = np.diff(array_in)
    zs = np.arange(len(array_in))
    plt.plot(array_in, zs)
    plt.ylim(max(zs), 0)
    plt.show()
    plt.pause(1)
    while any(delt < 0 for delt in delts):
        print(min(delts))
        neg_idxs = get_neg_idxs(delts)
        for neg_idx in neg_idxs:
            array_in[neg_idx] = array_in[neg_idx+1] = np.mean([array_in[neg_idx], array_in[neg_idx+1]])
        delts = np.diff(array_in)
        plt.plot(array_in, zs)
        plt.ylim(max(zs),0)
        plt.show()
        plt.pause(0.1)

def gr_mix_profile(profile, j, rg):
    rG = param_dict["rG"]
    # pseudo code:
    # mu{var} = (var[j] - var[j+1])/2
    # d{var} = (1- rg / rG) * mu{var}
    # var[j] = var[j] -d{var}
    # var[j+1] = var[j+1] + d{var}
    mu_var = (profile[j] - profile[j+1])/2
    d_var = (1 - rg/rG) * mu_var
    profile[j] = profile[j] - d_var
    profile[j+1] = profile[j+1] + d_var
    return profile
# tvals = np.array([28.009063947300664, 28.016413890100434, 28.013720377314854, 28.013818740844727, 28.019866943359375, 28.006534933659818, 28.010179578953693, 28.016091116233373, 28.022831889873064, 28.02896424424696, 28.033334823296805, 28.036869599901067, 28.040096341921135, 28.043331609995185, 28.04689196476141, 28.050995897083855, 28.05533667762138, 28.060043681605237, 28.065340479448075, 28.071450641562542, 28.079365937256, 28.09242757423483, 28.10781556739358, 28.121776004279354, 28.13055497243925, 28.13094187185051, 28.124697961352023, 28.11437617696024, 28.102236730784927, 28.09053983493586, 28.081133747705096, 28.072746780318962, 28.064517847349933, 28.055880019635175, 28.046266368011857, 28.03563613580514, 28.025855170657994, 28.01503183402758, 28.00066784175834, 27.980264909694718, 27.88195658200228, 27.373217599205137, 26.602424534281997, 25.796897915751895, 25.18395827213388, 24.92679166354327, 24.776096308670752, 24.66041865803841, 24.562237908415806, 24.464033256572513, 24.353160568983927, 24.24612542412919, 24.14246439514079, 24.03736698138703, 23.92602268223619, 23.80247999399955, 23.660982929020967, 23.510659119201094, 23.362432814088436, 23.22722826323149, 23.113992918810858, 23.01570522766404, 22.925762039129005, 22.83886088674457, 22.749699304049525, 22.65073710791473, 22.530047362470825, 22.397472340565397, 22.265847790669525, 22.148009461254272, 22.056793100790685, 22.005034457749844, 21.999317964538932, 22.009097622707486, 22.025946913287044, 22.047037983313203, 22.069542979821563, 22.090634049847722, 22.10748334042728, 22.117262998595834, 22.112722068224794, 22.068493897682252, 21.98754929387351, 21.877726591347795, 21.746864124654344, 21.60280022834238, 21.45337323696114, 21.30642148505985, 21.169783307187743, 21.05129703789405, 20.955468912635425, 20.86698360544255, 20.782269072088614, 20.70073413788792, 20.621787628154777, 20.544838368203486, 20.46929518334835, 20.394566898903676, 20.32006234018377, 20.245190332502936, 20.169300301097127, 20.09176609901602, 20.012758579131695, 19.932543708639653, 19.851387454735377, 19.769555784614358, 19.687314665472087, 19.604930064504057, 19.52266794890575, 19.44079428587266, 19.35957504260028, 19.27927618628409, 19.20016368411959, 19.12250350330227, 19.04656161102761, 18.972603974491104, 18.900896560888242, 18.831705337414515, 18.765296271265413, 18.701935329636424, 18.64188847972304, 18.585421688720746, 18.532800923825036, 18.484292152231397, 18.440161341135322, 18.400046230960232, 18.360397948823287, 18.320738009897266, 18.28121110285525, 18.24196191637031, 18.20313513911553, 18.164875459763987, 18.127327566988757, 18.09063614946292, 18.054945895859554, 18.020401494851725, 17.987147635112528, 17.955329005315033, 17.925090294132318, 17.896576190237457, 17.869931382303534, 17.84530055900363, 17.822828409010807, 17.80265962099816, 17.784938883638755, 17.769810885605676, 17.757420315572, 17.747911862210806, 17.741430214195166, 17.73812006019817, 17.73788826196289, 17.73935237536621, 17.74216326904297, 17.74619850842285, 17.751335658935545, 17.757452286010743, 17.764425955078124, 17.772134231567385, 17.780454680908203, 17.789264868530275, 17.79844235986328, 17.807864720336912, 17.81740951538086, 17.826954310424806, 17.836376670898435, 17.845554162231444, 17.854364349853515, 17.862684799194334, 17.870393075683594, 17.877366744750976, 17.883483371826173, 17.888620522338865, 17.89265576171875, 17.89546665539551, 17.89693076879883, 17.895177627319335, 17.880114821899415, 17.851196334838868, 17.80968177355957, 17.7568307454834, 17.693902858032228, 17.62215771862793, 17.542854934692382, 17.45725411364746, 17.366614862915036, 17.272196789916993, 17.175259502075193, 17.077062606811523, 16.97886571154785, 16.881928423706054, 16.787510350708008, 16.696871099975585, 16.611270278930665, 16.531967494995115, 16.460222355590822, 16.39729446813965, 16.344443440063475, 16.30292887878418, 16.274010391723632, 16.258947586303712, 16.257152959838866, 16.258295284301756, 16.26048838806152, 16.263636745727542, 16.26764483190918, 16.27241712121582, 16.277858088256835, 16.2838722076416, 16.290363953979494, 16.297237801879884, 16.30439822595215, 16.311749700805663, 16.319196701049805, 16.326643701293943, 16.33399517614746, 16.341155600219725, 16.34802944812012, 16.35452119445801, 16.360535313842774, 16.36597628088379, 16.370748570190433, 16.37475665637207, 16.377905014038085, 16.38009811779785, 16.381240442260744, 16.38070726045019, 16.375396909371297, 16.365116734709936, 16.350222069505854, 16.33106824679879, 16.308010599628496, 16.281404461034718, 16.2516051640572, 16.218968041735685, 16.183848427109925, 16.146601653219662, 16.10758305310464, 16.067147959804615, 16.02565170635932, 15.983449625808504, 15.940897051191918, 15.898349315549305, 15.85616175192041, 15.81468969334498, 15.77428847286276, 15.735313423513498, 15.698119878336936, 15.663063170372823, 15.630498632660903, 15.600781598240923])
# tvals2 = np.array([1, 2, 3, 4, 5, 6, 4, 8, 3, 5, 2, 11, 6, 4, 20, 9, 10, 9, 11, 12])
# smooth_inversions(-tvals, direction = "decreasing")
