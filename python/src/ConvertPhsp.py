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
    kAmpGen_GeV = [-0.22605460233259722, 0.37058687639201848, -0.046885439376411875, 0.65905276036464722]
    osPi1AmpGen_GeV = [0.075397408921232992, 0.24469544143911467, 0.20952672690121868, 0.35908482669738223]
    osPi2AmpGen_GeV = [0.07358860140319394, -0.24208436188963289, -0.30165403210059527, 0.41772611931236503]
    ssPiAmpGen_GeV = [0.077068592008170317, -0.37319795594150029, 0.13901274457578858, 0.42897629362560541]

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
    print("PhspPoint4Body testPoint [GeV] [rad -pi to pi]= {{ {}, {}, {}, {}, {} }};".format(m12_MeV/1000.0, m34_MeV/1000.0, c12, c34, phi_rad_neg_pi_to_pi))
    print("phi, 0 to 2pi = {}".format(phi_rad_0_to_2pi))

    print("PhspPoint4Body testPoint M12 = {}", (k + osPi1).M() / 1000.0)
    print("PhspPoint4Body testPoint M13 = {}", (k + osPi2).M() / 1000.0)
    print("PhspPoint4Body testPoint M14 = {}", (k + ssPi).M() / 1000.0)
    print("PhspPoint4Body testPoint M23 = {}",  (osPi1 + osPi2).M() / 1000.0)
    print("PhspPoint4Body testPoint M24 = {}",  (osPi1 + ssPi).M() / 1000.0)
    print("PhspPoint4Body testPoint M34 = {}",  (osPi2 + ssPi).M() / 1000.0)

    # # cross check boosts in GooFit

    # # boost k, os pi2 to kpi2 rest frame
    # vec_kPi2 = k + osPi2
    # kPi2Boost = vec_kPi2.BoostVector()
    # k.Boost(kPi2Boost)
    # osPi2.Boost(kPi2Boost)
    # print("In K Pi2 rest frame:")
    # print("px, py, pz, E [GeV]")
    # print("K = {}, {}, {}, {}".format(
    #     k.Px()/1000.0, 
    #     k.Py()/1000.0, 
    #     k.Pz()/1000.0, 
    #     k.E()/1000.0))
    # print("m(K) = {}".format(k.M()))
    # print("OS Pi 2 = {}, {}, {}, {}".format(
    #     osPi2.Px()/1000.0, 
    #     osPi2.Py()/1000.0, 
    #     osPi2.Pz()/1000.0, 
    #     osPi2.E()/1000.0))
    # print("m(OS Pi2) = {}".format(osPi2.M()))

    # # boost os pi1, ss pi to pi pi rest frame
    # vec_PiPi = osPi1 + ssPi
    # piPiBoost = vec_PiPi.BoostVector()
    # osPi1.Boost(piPiBoost)
    # ssPi.Boost(piPiBoost)
    # print("In Pi Pi rest frame:")
    # print("px, py, pz, E [GeV]")
    # print("OS Pi 1 = {}, {}, {}, {}".format(
    #     osPi1.Px()/1000.0, 
    #     osPi1.Py()/1000.0, 
    #     osPi1.Pz()/1000.0, 
    #     osPi1.E()/1000.0))
    # print("m(OS Pi 1) = {}".format(osPi1.M()))
    # print("SS Pi = {}, {}, {}, {}".format(
    #     ssPi.Px()/1000.0, 
    #     ssPi.Py()/1000.0, 
    #     ssPi.Pz()/1000.0, 
    #     ssPi.E()/1000.0))
    # print("m(SS Pi 1) = {}".format(ssPi.M()))

    # test if output vecs from GooFit get4Vecs gives same invariants
    testK = [0.389060, -0.140074, -0.140162, 0.659053]
    testOSPi1 = [0.264945, 0.140074, 0.140162, 0.359085]
    testOSPi2 = [-0.319720, -0.229770, 0.000000, 0.417726]
    testSSPi = [-0.334285, 0.229770, 0.000000, 0.428976]
    testKVec = py_k3pi_utilities.utils.ampGentoTLorentzVector(testK)
    testOSPi1Vec = py_k3pi_utilities.utils.ampGentoTLorentzVector(testOSPi1)
    testOSPi2Vec = py_k3pi_utilities.utils.ampGentoTLorentzVector(testOSPi2)
    testSSPiVec = py_k3pi_utilities.utils.ampGentoTLorentzVector(testSSPi)
    print("From GooFit M12 = {}", (testKVec + testOSPi1Vec).M() / 1000.0)
    print("From GooFit M13 = {}", (testKVec + testOSPi2Vec).M() / 1000.0)
    print("From GooFit M14 = {}", (testKVec + testSSPiVec).M() / 1000.0)
    print("From GooFit M23 = {}",  (testOSPi1Vec + testOSPi2Vec).M() / 1000.0)
    print("From GooFit M24 = {}",  (testOSPi1Vec + testSSPiVec).M() / 1000.0)
    print("From GooFit M34 = {}",  (testOSPi2Vec + testSSPiVec).M() / 1000.0)
    phspTest = ROOT.K3PiStudies.K3PiStudiesUtils.calc_phsp(d0, testKVec, testOSPi1Vec, testSSPiVec, testOSPi2Vec)
    m12_MeVTest = phspTest[0]
    m34_MeVTest = phspTest[1]
    c12Test = phspTest[2]
    c34Test = phspTest[3]
    phi_rad_0_to_2piTest = phspTest[4]
    phi_rad_neg_pi_to_piTest = ROOT.K3PiStudies.K3PiStudiesUtils.changeAngleRange_neg_pi_to_pi(phi_rad_0_to_2piTest)
    print("GooFit computed testPoint [GeV] [rad -pi to pi]= {{ {}, {}, {}, {}, {} }}".format(m12_MeVTest/1000.0, m34_MeVTest/1000.0, c12Test, c34Test, phi_rad_neg_pi_to_piTest))
# end main

if __name__ == "__main__":
    main()