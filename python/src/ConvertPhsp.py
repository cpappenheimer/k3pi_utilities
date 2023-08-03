import ROOT
import sys

PY_UTILS_PKG_DIR = "/data/home/pappenheimer/develop_k3pistudies/K3PiStudies/extern/k3pi_utilities/python/src"

sys.path.append(PY_UTILS_PKG_DIR)   
import py_k3pi_utilities.utils

def main():
    ##### CONFIG
    UTILS_BUILD_DIR = "/data/home/pappenheimer/develop_k3pistudies/build/extern/k3pi_utilities/"
    UTILS_INC_DIR = "/data/home/pappenheimer/develop_k3pistudies/K3PiStudies/extern/k3pi_utilities/include/"

    D0_M_MEV = 1864.84
    #############

    py_k3pi_utilities.utils.loadK3PiCUtils(UTILS_BUILD_DIR, UTILS_INC_DIR)

    # px, py, pz, pE
    kAmpGen_GeV = [-0.226055, 0.370587, -0.0468854, 0.659053]
    osPi1AmpGen_GeV = [0.0753974, 0.244695, 0.209527, 0.359085]
    osPi2AmpGen_GeV = [0.0735886, -0.242084, -0.301654, 0.417726]
    ssPiAmpGen_GeV = [0.0770686, -0.373198, 0.139013, 0.428976]

    # convert 4 vecs from AmpGen to MeV, TLorentzVector
    k = py_k3pi_utilities.utils.ampGentoTLorentzVector(kAmpGen_GeV)
    osPi1 = py_k3pi_utilities.utils.ampGentoTLorentzVector(osPi1AmpGen_GeV)
    osPi2 = py_k3pi_utilities.utils.ampGentoTLorentzVector(osPi2AmpGen_GeV)
    ssPi = py_k3pi_utilities.utils.ampGentoTLorentzVector(ssPiAmpGen_GeV)
    d0 = ROOT.K3PiStudies.K3PiStudiesUtils.toTLorentzVector(D0_M_MEV, 0.0, 0.0, 0.0)

    # convert to m12, m34, c12, c34, phi
    phsp = ROOT.K3PiStudies.K3PiStudiesUtils.calc_phsp(d0, k, osPi1, ssPi, osPi2)
    m12_MeV = phsp[0]
    m34_MeV = phsp[1]
    c12 = phsp[2]
    c34 = phsp[3]
    phi_rad_0_to_2pi = phsp[4]
    phi_rad_neg_pi_to_pi = ROOT.K3PiStudies.K3PiStudiesUtils.changeAngleRange_neg_pi_to_pi(phi_rad_0_to_2pi)
    print("PhspPoint4Body testPoint = {{ {}, {}, {}, {}, {} }};".format(m12_MeV, m34_MeV, c12, c34, phi_rad_neg_pi_to_pi))
# end main

if __name__ == "__main__":
    main()