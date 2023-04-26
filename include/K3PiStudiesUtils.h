#pragma once

#include <string>
#include <vector>
#include <utility>

#include <TMath.h>
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>
#include <TH1.h>

namespace K3PiStudies
{

	struct InvalidDecayError : std::exception
	{
		const std::string _msg;

		InvalidDecayError(const std::string &msg) : _msg(msg)
		{
		}

		const char *what() const noexcept override
		{
			return _msg.c_str();
		}
	}; // end InvalidDecayError

	struct ComputationError : std::exception
	{
		const std::string _msg;

		ComputationError(const std::string &msg) : _msg(msg)
		{
		}

		const char *what() const noexcept override
		{
			return _msg.c_str();
		}
	}; // end ComputationError

	// names for the particles used in the "*D0Fit*" vars in the ntuple
	enum D0Fit_PNames
	{
		D0Fit_D0_Kplus,
		D0Fit_D0_piplus_0,
		D0Fit_D0_piplus_1,
		D0Fit_D0_piplus
	};

	// names for the particles used in the "*ReFit*" vars in the ntuple
	enum ReFit_PNames
	{
		ReFit_D0_Kplus,
		ReFit_D0_piplus_0,
		ReFit_D0_piplus_1,
		ReFit_D0_piplus
	};

	class K3PiStudiesUtils final
	{
	public:
		static inline const std::string _ALL_REGION_FLAG = "ALL";
		static inline const std::string _SIG_REGION_FLAG = "SIGNAL";
		static constexpr double _DELTAM_PDG_MEV = 145.4258;
		static constexpr double _SIG_REGION_LOW_MD0_BOUND_MEV = 1850;
		static constexpr double _SIG_REGION_HIGH_MD0_BOUND_MEV = 1882;
		static constexpr double _SIG_REGION_LOW_DELTAM_BOUND_MEV = _DELTAM_PDG_MEV - 0.7;
		static constexpr double _SIG_REGION_HIGH_DELTAM_BOUND_MEV = _DELTAM_PDG_MEV + 0.7;
		static constexpr double _ALL_REGS_D0_MASS_AXIS_MIN_MEV = 1800.0;
		static constexpr double _ALL_REGS_D0_MASS_AXIS_MAX_MEV = 1920.0;
		static constexpr double _ALL_REGS_DELTAM_AXIS_MIN_MEV = 138.0;
		static constexpr double _ALL_REGS_DELTAM_AXIS_MAX_MEV = 158.0;
		static const double _C_M_PER_SEC;
		static constexpr double _MM_TO_M = 1.0 / 1000.0;
		static const double _SEC_TO_NS;
		static const unsigned int _KAON_ID = 321;
		static const unsigned int _PION_ID = 211;
		static constexpr double _PI = TMath::Pi();
		static constexpr double _GEV_TO_MEV = 1000.0;
		static constexpr double _D0_LIFETIME_PS = 0.4103;
		static constexpr double _NS_TO_PS = 1000.0;
		static const std::string _RS_FLAG;
		static const std::string _WS_FLAG;
		static const std::string _BOTH_FLAG;
		static constexpr double _KAON_MASS = 493.677;
		static constexpr double _PION_MASS = 139.57061;
		static const std::string _REFIT_FLAG;
		static const std::string _D0_FIT_FLAG;
		static const std::string _P_FLAG;
		static constexpr double _COMPARE_EPS = std::numeric_limits<double>::epsilon();

		K3PiStudiesUtils() = default;
		K3PiStudiesUtils(const K3PiStudiesUtils &copyMe) = default;
		K3PiStudiesUtils(K3PiStudiesUtils &&moveMe) = default;
		~K3PiStudiesUtils() = default;
		K3PiStudiesUtils &operator=(const K3PiStudiesUtils &copyMe) = default;
		K3PiStudiesUtils &operator=(K3PiStudiesUtils &&moveMe) = default;

		static void makeTLegendBkgTransparent(TLegend& leg);

		static void makeTPaveTextBkgTransparent(TPaveText& pt);

		static std::pair<double,double> getRegionAxisBoundsDeltaMMeV(const std::string& regionName);

		static std::pair<double,double> getRegionAxisBoundsMD0MeV(const std::string& regionName);

		static bool isInDeltaMRegion(const std::string& regionName, double deltaMMeV);

		static bool isInD0MassRegion(const std::string& regionName, double d0MassMeV);

		static std::string printRegionBoundsDeltaM(const std::string& regionName);

		static std::string printRegionBoundsMD0(const std::string& regionName);

		static void makeNormalizedComparisonPlot(
			TH1 *const h1,
			TH1 *const h2,
			const TString &legLine1,
			const TString &legLine2,
			bool addNumEntries,
			const TString &saveName);

		static double angleBetweenDecayPlanesKutschke(
			const TVector3 &d4_motherRestFrame,
			const TVector3 &d5_motherRestFrame,
			const TVector3 &d6_motherRestFrame,
			const TVector3 &d7_motherRestFrame);

		static TString makeTitleStr(
			const TString &title,
			const TString &xLabel,
			const TString &yLabel);

		static TString makeYAxisLabel(
			int numBins,
			double axisMin,
			double axisMax,
			const TString &unit,
			bool normalizedPlot);

		static void changeToRainbowPalette();

		static double verifyAngle(
			const TVector3 &v1,
			const TVector3 &v2,
			double v1v2Angle,
			bool v1v2AngleIsNegPiToPi,
			const std::string &angleName,
			bool printDiff);

		static float helicity_angle_func(
			float d0_px,
			float d0_py,
			float d0_pz,
			float d0_m,
			float pis_px,
			float pis_py,
			float pis_pz,
			float pis_m);

		static float helicity_angle_func(
			float d0_px,
			float d0_py,
			float d0_pz,
			float d0_m,
			const ROOT::RVec<float> &pis_px,
			const ROOT::RVec<float> &pis_py,
			const ROOT::RVec<float> &pis_pz,
			float pis_m);

		static double compute_delta_angle(
			double extra_px,
			double extra_py,
			double extra_pz,
			double extra_m,
			double d_px,
			double d_py,
			double d_pz,
			double d_m);

		static double compute_delta_angle(
			double extra_px,
			double extra_py,
			double extra_pz,
			double d_px,
			double d_py,
			double d_pz);

		static double getPhi(
			double px,
			double py,
			double pz,
			double pE);

		static double getEta(
			double px,
			double py,
			double pz,
			double pE);

		static double getPT(
			double px,
			double py,
			double pz,
			double pE);

		static bool isWithinDecayTimeBin(
			double dtime,
			const std::pair<double, double> &decayTimeLimits);

		static unsigned int determineQuadrant(double sin2ThetaA, double sin2ThetaC);

		static bool isKPi1LowerMassPair(
			const TLorentzVector &kminus_4vec,
			const TLorentzVector &piplus1_4vec,
			const TLorentzVector &piplus2_4vec);

		static bool isKaonNeg(
			int kaonInd,
			int D0_P0_ID,
			int D0_P1_ID,
			int D0_P2_ID,
			int D0_P3_ID);

		static int findSSPion(
			bool kaonIsNeg,
			int D0_P0_ID,
			int D0_P1_ID,
			int D0_P2_ID,
			int D0_P3_ID);

		static std::vector<int> findOSPions(
			bool kaonIsNeg,
			int D0_P0_ID,
			int D0_P1_ID,
			int D0_P2_ID,
			int D0_P3_ID);

		static int findKaon(
			int D0_P0_ID,
			int D0_P1_ID,
			int D0_P2_ID,
			int D0_P3_ID);

		static bool isD0(int dStarPiID);

		static bool isRS(bool isD0, bool isKaonNeg);

		static std::vector<std::string> buildListFromCommaSepStr(const std::string &filesString);

		static double cTauMMToTauNS(double cTauMM);

		static double tauNSToTauPS(double tauNS);

		static std::string d0TimeBinToString(
			const std::pair<double, double> &decayTimeLimits,
			const std::string &unit);

		static std::vector<std::pair<double, double>> makeTimeBins(const std::vector<double> &upperBinEdges);

		static double getD0Part_M(
			int ind,
			double D0_P0_M,
			double D0_P1_M,
			double D0_P2_M,
			double D0_P3_M);

		static double getD0Part_PX(
			int ind,
			double D0_P0_PX,
			double D0_P1_PX,
			double D0_P2_PX,
			double D0_P3_PX);

		static double getD0Part_PY(
			int ind,
			double D0_P0_PY,
			double D0_P1_PY,
			double D0_P2_PY,
			double D0_P3_PY);

		static double getD0Part_PZ(
			int ind,
			double D0_P0_PZ,
			double D0_P1_PZ,
			double D0_P2_PZ,
			double D0_P3_PZ);

		static double getD0Fit_PE(
			D0Fit_PNames pName,
			double Dst_D0Fit_D0_Kplus_PE,
			double Dst_D0Fit_D0_piplus_0_PE,
			double Dst_D0Fit_D0_piplus_1_PE,
			double Dst_D0Fit_D0_piplus_PE);

		static double getD0Fit_PX(
			D0Fit_PNames pName,
			double Dst_D0Fit_D0_Kplus_PX,
			double Dst_D0Fit_D0_piplus_0_PX,
			double Dst_D0Fit_D0_piplus_1_PX,
			double Dst_D0Fit_D0_piplus_PX);

		static double getD0Fit_PY(
			D0Fit_PNames pName,
			double Dst_D0Fit_D0_Kplus_PY,
			double Dst_D0Fit_D0_piplus_0_PY,
			double Dst_D0Fit_D0_piplus_1_PY,
			double Dst_D0Fit_D0_piplus_PY);

		static double getD0Fit_PZ(
			D0Fit_PNames pName,
			double Dst_D0Fit_D0_Kplus_PZ,
			double Dst_D0Fit_D0_piplus_0_PZ,
			double Dst_D0Fit_D0_piplus_1_PZ,
			double Dst_D0Fit_D0_piplus_PZ);

		static D0Fit_PNames findD0FitKaon(
			int Dst_D0Fit_D0_Kplus_ID,
			int Dst_D0Fit_D0_piplus_0_ID,
			int Dst_D0Fit_D0_piplus_1_ID,
			int Dst_D0Fit_D0_piplus_ID);

		static D0Fit_PNames indexToD0Fit_PName(int index);

		static std::vector<D0Fit_PNames> findD0FitOSPions(
			bool kaonIsNeg,
			int Dst_D0Fit_D0_Kplus_ID,
			int Dst_D0Fit_D0_piplus_0_ID,
			int Dst_D0Fit_D0_piplus_1_ID,
			int Dst_D0Fit_D0_piplus_ID);

		static D0Fit_PNames findD0FitSSPion(
			bool kaonIsNeg,
			int Dst_D0Fit_D0_Kplus_ID,
			int Dst_D0Fit_D0_piplus_0_ID,
			int Dst_D0Fit_D0_piplus_1_ID,
			int Dst_D0Fit_D0_piplus_ID);

		static bool isD0FitKaonNeg(
			D0Fit_PNames kaonName,
			int Dst_D0Fit_D0_Kplus_ID,
			int Dst_D0Fit_D0_piplus_0_ID,
			int Dst_D0Fit_D0_piplus_1_ID,
			int Dst_D0Fit_D0_piplus_ID);

		static double getProbNNx(
			double D0_P0_ProbNNx,
			double D0_P1_ProbNNx,
			double D0_P2_ProbNNx,
			double D0_P3_ProbNNx,
			int ind);

		static ReFit_PNames indexToReFit_PName(int index);

		static ReFit_PNames findReFitKaon(
			int Dst_ReFit_D0_Kplus_ID,
			int Dst_ReFit_D0_piplus_0_ID,
			int Dst_ReFit_D0_piplus_1_ID,
			int Dst_ReFit_D0_piplus_ID);

		static std::vector<ReFit_PNames> findReFitOSPions(
			bool kaonIsNeg,
			int Dst_ReFit_D0_Kplus_ID,
			int Dst_ReFit_D0_piplus_0_ID,
			int Dst_ReFit_D0_piplus_1_ID,
			int Dst_ReFit_D0_piplus_ID);

		static ReFit_PNames findReFitSSPion(
			bool kaonIsNeg,
			int Dst_ReFit_D0_Kplus_ID,
			int Dst_ReFit_D0_piplus_0_ID,
			int Dst_ReFit_D0_piplus_1_ID,
			int Dst_ReFit_D0_piplus_ID);

		static bool isReFitKaonNeg(
			ReFit_PNames kaonName,
			int Dst_ReFit_D0_Kplus_ID,
			int Dst_ReFit_D0_piplus_0_ID,
			int Dst_ReFit_D0_piplus_1_ID,
			int Dst_ReFit_D0_piplus_ID);

		static double getReFit_PE(
			ReFit_PNames pName,
			double Dst_ReFit_D0_Kplus_PE,
			double Dst_ReFit_D0_piplus_0_PE,
			double Dst_ReFit_D0_piplus_1_PE,
			double Dst_ReFit_D0_piplus_PE);

		static double getReFit_PX(
			ReFit_PNames pName,
			double Dst_ReFit_D0_Kplus_PX,
			double Dst_ReFit_D0_piplus_0_PX,
			double Dst_ReFit_D0_piplus_1_PX,
			double Dst_ReFit_D0_piplus_PX);

		static double getReFit_PY(
			ReFit_PNames pName,
			double Dst_ReFit_D0_Kplus_PY,
			double Dst_ReFit_D0_piplus_0_PY,
			double Dst_ReFit_D0_piplus_1_PY,
			double Dst_ReFit_D0_piplus_PY);

		static double getReFit_PZ(
			ReFit_PNames pName,
			double Dst_ReFit_D0_Kplus_PZ,
			double Dst_ReFit_D0_piplus_0_PZ,
			double Dst_ReFit_D0_piplus_1_PZ,
			double Dst_ReFit_D0_piplus_PZ);

		static std::vector<double> calc_phsp(
			double K_D0Fit_PT,
			double K_D0Fit_ETA,
			double K_D0Fit_PHI,
			double Pi_SS_D0Fit_PT,
			double Pi_SS_D0Fit_ETA,
			double Pi_SS_D0Fit_PHI,
			double Pi_OS1_D0Fit_PT,
			double Pi_OS1_D0Fit_ETA,
			double Pi_OS1_D0Fit_PHI,
			double Pi_OS2_D0Fit_PT,
			double Pi_OS2_D0Fit_ETA,
			double Pi_OS2_D0Fit_PHI,
			bool pi1GoesWithK,
			bool verifyAngles,
			bool printDiff);

		static double getD0Part_PE(
			int ind,
			double D0_P0_PE,
			double D0_P1_PE,
			double D0_P2_PE,
			double D0_P3_PE);

		static bool areDoublesEqual(
			std::function<bool(double, double)> isEqualFunc,
			double d1,
			double d2,
			const std::string &varName,
			bool printDiff);

		static double changeAngleRange_neg_pi_to_pi(double angle_0_to_2pi);

		static double radToDeg(double angleRad);

		static double changeAngleRange_0_to_2pi(double angle_neg_pi_to_pi);

		static void silenceROOTHistSaveMsgs();

		static bool combinedToleranceCompare(double x, double y);

		static bool isExactlyEqual(double d1, double d2);

	}; // end K3PiStudiesUtils class

	struct Phsp4Body
	{
		const double _m12_MeV;
		const double _m34_MeV;
		const double _cos12;
		const double _cos34;
		const double _phi_rad;

		/** difference between phi values we computed when using 2 different methods (ours and TVector3's .Angle()) to calculate (use to check for calculation errors) */
		const double _phi_diff;
		const double _phiA_diff;

		// alternative sets of variables
		const double _m14_MeV;
		const double _m32_MeV;
		const double _phiA_rad;

		/** difference between phi values we computed when using 2 different methods (ours and Eq. 42 in Kutschke's An Angular Distribution Cookbook) to calculate (use to check for calculation errors) */
		const double _phi_diff2;
		const double _phiA_diff2;

		Phsp4Body(
			double m12_MeV,
			double m34_MeV,
			double cos12,
			double cos34,
			double phi_rad,
			double phi_diff,
			double m14_MeV,
			double m32_MeV,
			double phiA_rad,
			double phiA_diff,
			double phi_diff2,
			double phiA_diff2) : _m12_MeV(m12_MeV), _m34_MeV(m34_MeV), _cos12(cos12), _cos34(cos34), _phi_rad(phi_rad), _phi_diff(phi_diff), _m14_MeV(m14_MeV), _m32_MeV(m32_MeV), _phiA_rad(phiA_rad), _phiA_diff(phiA_diff), _phi_diff2(phi_diff2), _phiA_diff2(phi_diff2)
		{
		}

		int compare5(const Phsp4Body &other, std::function<bool(double, double)> isEqualFunc, int eventNum, bool printSanityChecks) const
		{
			std::string evt = std::to_string(eventNum);
			bool isEqual[5];
			isEqual[0] = K3PiStudiesUtils::areDoublesEqual(isEqualFunc, this->_m12_MeV, other._m12_MeV, "Event " + evt + " m12", printSanityChecks);
			isEqual[1] = K3PiStudiesUtils::areDoublesEqual(isEqualFunc, this->_m34_MeV, other._m34_MeV, "Event " + evt + " m34", printSanityChecks);
			isEqual[2] = K3PiStudiesUtils::areDoublesEqual(isEqualFunc, this->_cos12, other._cos12, "Event " + evt + " cos12", printSanityChecks);
			isEqual[3] = K3PiStudiesUtils::areDoublesEqual(isEqualFunc, this->_cos34, other._cos34, "Event " + evt + " cos34", printSanityChecks);
			isEqual[4] = K3PiStudiesUtils::areDoublesEqual(isEqualFunc, this->_phi_rad, other._phi_rad, "Event " + evt + " phi", printSanityChecks);

			int numDiffs = 0;
			for (int i = 0; i < 5; i++)
			{
				if (!isEqual[i])
				{
					numDiffs++;
				}
			}

			return numDiffs;
		}
	};

} // end namespace K3PiStudies