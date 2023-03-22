#include <sstream>

#include <TLorentzVector.h>
#include <Math/Vector4D.h>

#include "K3PiStudiesUtils.h"

namespace K3PiStudies
{

	const std::string K3PiStudiesUtils::_RS_FLAG = "RS";
	const std::string K3PiStudiesUtils::_WS_FLAG = "WS";
	const std::string K3PiStudiesUtils::_BOTH_FLAG = "BOTH";
	const double K3PiStudiesUtils::_C_M_PER_SEC = 3.0 * TMath::Power(10, 8);
	const double K3PiStudiesUtils::_SEC_TO_NS = TMath::Power(10, 9);
	const std::string K3PiStudiesUtils::_REFIT_FLAG = "REFIT";
	const std::string K3PiStudiesUtils::_D0_FIT_FLAG = "D0_FIT";
	const std::string K3PiStudiesUtils::_P_FLAG = "P";


	bool K3PiStudiesUtils::areDoublesEqual(double d1, double d2, const std::string& varName, bool printDiff)
	{
		if (d1 != d2)
		{
			if (printDiff)
			{
				std::cout << "Found difference for " << varName << ": " << d1 << ", " << d2 << std::endl;
			}
			
			return false;
		}
		else
		{
			return true;
		}
	}

	/**
	 * For the D0_P0_*, D0_P1_*, D0_P2_*, D0_P3_* vars
	 */
	double K3PiStudiesUtils::getD0Part_M(
		int ind,
		double D0_P0_M,
		double D0_P1_M,
		double D0_P2_M,
		double D0_P3_M)
	{
		double m;
		switch (ind)
		{
		case 0:
			m = D0_P0_M;
			break;
		case 1:
			m = D0_P1_M;
			break;
		case 2:
			m = D0_P2_M;
			break;
		case 3:
			m = D0_P3_M;
			break;
		default:
			throw InvalidDecayError("getD0Part_M: Cannot find daughter with index " + std::to_string(ind) + " in daughters.");
		}

		return m;
	}


	/**
	 * For the D0_P0_*, D0_P1_*, D0_P2_*, D0_P3_* vars
	 */
	double K3PiStudiesUtils::getD0Part_PE(
		int ind,
		double D0_P0_PE,
		double D0_P1_PE,
		double D0_P2_PE,
		double D0_P3_PE)
	{
		double pE;
		switch (ind)
		{
		case 0:
			pE = D0_P0_PE;
			break;
		case 1:
			pE = D0_P1_PE;
			break;
		case 2:
			pE = D0_P2_PE;
			break;
		case 3:
			pE = D0_P3_PE;
			break;
		default:
			throw InvalidDecayError("getD0Part_PE: Cannot find daughter with index " + std::to_string(ind) + " in daughters.");
		}

		return pE;
	}

	/**
	 * For the D0_P0_*, D0_P1_*, D0_P2_*, D0_P3_* vars
	 */
	double K3PiStudiesUtils::getD0Part_PZ(
		int ind,
		double D0_P0_PZ,
		double D0_P1_PZ,
		double D0_P2_PZ,
		double D0_P3_PZ)
	{
		double pz;
		switch (ind)
		{
		case 0:
			pz = D0_P0_PZ;
			break;
		case 1:
			pz = D0_P1_PZ;
			break;
		case 2:
			pz = D0_P2_PZ;
			break;
		case 3:
			pz = D0_P3_PZ;
			break;
		default:
			throw InvalidDecayError("getD0Part_PZ: Cannot find daughter with index " + std::to_string(ind) + " in daughters.");
		}

		return pz;
	}

	/**
	 * For the D0_P0_*, D0_P1_*, D0_P2_*, D0_P3_* vars
	 */
	double K3PiStudiesUtils::getD0Part_PY(
		int ind,
		double D0_P0_PY,
		double D0_P1_PY,
		double D0_P2_PY,
		double D0_P3_PY)
	{
		double py;
		switch (ind)
		{
		case 0:
			py = D0_P0_PY;
			break;
		case 1:
			py = D0_P1_PY;
			break;
		case 2:
			py = D0_P2_PY;
			break;
		case 3:
			py = D0_P3_PY;
			break;
		default:
			throw InvalidDecayError("getD0Part_PY: Cannot find daughter with index " + std::to_string(ind) + " in daughters.");
		}

		return py;
	}

	/**
	 * For the D0_P0_*, D0_P1_*, D0_P2_*, D0_P3_* vars
	 */
	double K3PiStudiesUtils::getD0Part_PX(
		int ind,
		double D0_P0_PX,
		double D0_P1_PX,
		double D0_P2_PX,
		double D0_P3_PX)
	{
		double px;
		switch (ind)
		{
		case 0:
			px = D0_P0_PX;
			break;
		case 1:
			px = D0_P1_PX;
			break;
		case 2:
			px = D0_P2_PX;
			break;
		case 3:
			px = D0_P3_PX;
			break;
		default:
			throw InvalidDecayError("getD0Part_PX: Cannot find daughter with index " + std::to_string(ind) + " in daughters.");
		}

		return px;
	}

	double K3PiStudiesUtils::getPhi(
		double px,
		double py,
		double pz,
		double pE)
	{
		ROOT::Math::PxPyPzEVector v(px, py, pz, pE);
		return v.Phi();
	}

	double K3PiStudiesUtils::getEta(
		double px,
		double py,
		double pz,
		double pE)
	{
		ROOT::Math::PxPyPzEVector v(px, py, pz, pE);
		return v.Eta();
	}

	double K3PiStudiesUtils::getPT(
		double px,
		double py,
		double pz,
		double pE)
	{
		ROOT::Math::PxPyPzEVector v(px, py, pz, pE);
		return v.Pt();
	}

	/*
	 * Function to calculate phase space from John's apply_full_selection.py code
	 */
	std::vector<double> K3PiStudiesUtils::calc_phsp(
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
		bool ordered)
	{
		TLorentzVector piGoesWithPi, ssPiVec, kVec, piGoesWithK;
		ssPiVec.SetPtEtaPhiM(Pi_SS_D0Fit_PT, Pi_SS_D0Fit_ETA, Pi_SS_D0Fit_PHI, _PION_MASS);
		kVec.SetPtEtaPhiM(K_D0Fit_PT, K_D0Fit_ETA, K_D0Fit_PHI, _KAON_MASS);

		if (ordered)
		{
			TLorentzVector pi1, pi2;
			pi1.SetPtEtaPhiM(Pi_OS1_D0Fit_PT, Pi_OS1_D0Fit_ETA, Pi_OS1_D0Fit_PHI, _PION_MASS);
			pi2.SetPtEtaPhiM(Pi_OS2_D0Fit_PT, Pi_OS2_D0Fit_ETA, Pi_OS2_D0Fit_PHI, _PION_MASS);

			if (isKPi1LowerMassPair(kVec, pi1, pi2))
			{
				piGoesWithPi.SetPtEtaPhiM(pi2.Pt(), pi2.Eta(), pi2.Phi(), _PION_MASS);
				piGoesWithK.SetPtEtaPhiM(pi1.Pt(), pi1.Eta(), pi1.Phi(), _PION_MASS);
			}
			else
			{
				piGoesWithPi.SetPtEtaPhiM(pi1.Pt(), pi1.Eta(), pi1.Phi(), _PION_MASS);
				piGoesWithK.SetPtEtaPhiM(pi2.Pt(), pi2.Eta(), pi2.Phi(), _PION_MASS);
			}
		}
		else
		{
			// Randomise the assignment of piGoesWithPi and piGoesWithK, which should be the two OS pions.
			if (((double)std::rand() / (RAND_MAX)) < 0.5) // FIXME use seed
			{
				piGoesWithPi.SetPtEtaPhiM(Pi_OS1_D0Fit_PT, Pi_OS1_D0Fit_ETA, Pi_OS1_D0Fit_PHI, _PION_MASS);
				piGoesWithK.SetPtEtaPhiM(Pi_OS2_D0Fit_PT, Pi_OS2_D0Fit_ETA, Pi_OS2_D0Fit_PHI, _PION_MASS);
			}
			else
			{
				piGoesWithK.SetPtEtaPhiM(Pi_OS1_D0Fit_PT, Pi_OS1_D0Fit_ETA, Pi_OS1_D0Fit_PHI, _PION_MASS);
				piGoesWithPi.SetPtEtaPhiM(Pi_OS2_D0Fit_PT, Pi_OS2_D0Fit_ETA, Pi_OS2_D0Fit_PHI, _PION_MASS);
			}
		}

		// Boost everything to D0 restframe
		auto mum = piGoesWithPi + ssPiVec + kVec + piGoesWithK;
		double m12 = (kVec + piGoesWithK).M();
		double m34 = (piGoesWithPi + ssPiVec).M();
		double m13 = (piGoesWithPi + kVec).M();
		piGoesWithPi.Boost(-mum.BoostVector());
		ssPiVec.Boost(-mum.BoostVector());
		kVec.Boost(-mum.BoostVector());
		piGoesWithK.Boost(-mum.BoostVector());

		TLorentzVector piPiVec, kPiVec;
		piPiVec = piGoesWithPi + ssPiVec;
		kPiVec = kVec + piGoesWithK;

		TVector3 piGoesWithPin = piGoesWithPi.Vect().Unit();
		TVector3 ssPiVecn = ssPiVec.Vect().Unit();
		TVector3 kVecn = kVec.Vect().Unit();
		TVector3 piGoesWithKn = piGoesWithK.Vect().Unit();
		TVector3 piPiVecn = piPiVec.Vect().Unit();
		TVector3 kPiVecn = kPiVec.Vect().Unit();

		TVector3 n1 = piGoesWithPin.Cross(ssPiVecn);
		TVector3 n2 = kVecn.Cross(piGoesWithKn);
		TVector3 n3 = n1.Unit().Cross(n2.Unit());

		// Calculation of the angle Phi between the planes.
		double cosp = n1.Unit().Dot(n2.Unit());
		double sinp = n3.Dot(kPiVecn);
		double phi = acos(cosp);
		if (sinp < 0)
			phi *= -1;

		// Vectors in rest fram of their resonance.
		TLorentzVector piGoesWithPir = piGoesWithPi;
		TLorentzVector kVecr = kVec;
		piGoesWithPir.Boost(-piPiVec.BoostVector());
		kVecr.Boost(-kPiVec.BoostVector());
		TVector3 piGoesWithPirn = piGoesWithPir.Vect().Unit();
		TVector3 kVecrn = kVecr.Vect().Unit();

		// helicity angle for piPiVec and kPiVec frame
		double cos1 = kPiVecn.Dot(kVecrn);
		double cos2 = piPiVecn.Dot(piGoesWithPirn);
	
		// *12 = k, os pi that makes lower mass
		// *34 = ss pi, other os pi 
		std::vector<double> vars = {m12, m34, cos1, cos2, phi, m13};
		return vars;
	}

	bool K3PiStudiesUtils::isReFitKaonNeg(
		ReFit_PNames kaonName,
		int Dst_ReFit_D0_Kplus_ID,
		int Dst_ReFit_D0_piplus_0_ID,
		int Dst_ReFit_D0_piplus_1_ID,
		int Dst_ReFit_D0_piplus_ID)
	{
		bool isKaonNeg;
		switch (kaonName)
		{
		case ReFit_PNames::ReFit_D0_Kplus:
			isKaonNeg = Dst_ReFit_D0_Kplus_ID < 0;
			break;
		case ReFit_PNames::ReFit_D0_piplus_0:
			isKaonNeg = Dst_ReFit_D0_piplus_0_ID < 0;
			break;
		case ReFit_PNames::ReFit_D0_piplus_1:
			isKaonNeg = Dst_ReFit_D0_piplus_1_ID < 0;
			break;
		case ReFit_PNames::ReFit_D0_piplus:
			isKaonNeg = Dst_ReFit_D0_piplus_ID < 0;
			break;
		default:
			throw InvalidDecayError("isReFitKaonNeg: Cannot find daughter with name " + std::to_string(kaonName) + " in daughters.");
		}

		return isKaonNeg;
	}

	bool K3PiStudiesUtils::isD0FitKaonNeg(
		D0Fit_PNames kaonName,
		int Dst_D0Fit_D0_Kplus_ID,
		int Dst_D0Fit_D0_piplus_0_ID,
		int Dst_D0Fit_D0_piplus_1_ID,
		int Dst_D0Fit_D0_piplus_ID)
	{
		bool isKaonNeg;
		switch (kaonName)
		{
		case D0Fit_PNames::D0Fit_D0_Kplus:
			isKaonNeg = Dst_D0Fit_D0_Kplus_ID < 0;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_0:
			isKaonNeg = Dst_D0Fit_D0_piplus_0_ID < 0;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_1:
			isKaonNeg = Dst_D0Fit_D0_piplus_1_ID < 0;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus:
			isKaonNeg = Dst_D0Fit_D0_piplus_ID < 0;
			break;
		default:
			throw InvalidDecayError("isD0FitKaonNeg: Cannot find daughter with name " + std::to_string(kaonName) + " in daughters.");
		}

		return isKaonNeg;
	}

	double K3PiStudiesUtils::getD0Fit_PE(
		D0Fit_PNames pName,
		double Dst_D0Fit_D0_Kplus_PE,
		double Dst_D0Fit_D0_piplus_0_PE,
		double Dst_D0Fit_D0_piplus_1_PE,
		double Dst_D0Fit_D0_piplus_PE)
	{
		double pE;
		switch (pName)
		{
		case D0Fit_PNames::D0Fit_D0_Kplus:
			pE = Dst_D0Fit_D0_Kplus_PE;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_0:
			pE = Dst_D0Fit_D0_piplus_0_PE;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_1:
			pE = Dst_D0Fit_D0_piplus_1_PE;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus:
			pE = Dst_D0Fit_D0_piplus_PE;
			break;
		default:
			throw InvalidDecayError("getD0Fit_PE: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return pE;
	}

	double K3PiStudiesUtils::getD0Fit_PX(
		D0Fit_PNames pName,
		double Dst_D0Fit_D0_Kplus_PX,
		double Dst_D0Fit_D0_piplus_0_PX,
		double Dst_D0Fit_D0_piplus_1_PX,
		double Dst_D0Fit_D0_piplus_PX)
	{
		double px;
		switch (pName)
		{
		case D0Fit_PNames::D0Fit_D0_Kplus:
			px = Dst_D0Fit_D0_Kplus_PX;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_0:
			px = Dst_D0Fit_D0_piplus_0_PX;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_1:
			px = Dst_D0Fit_D0_piplus_1_PX;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus:
			px = Dst_D0Fit_D0_piplus_PX;
			break;
		default:
			throw InvalidDecayError("getD0Fit_PX: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return px;
	}

	double K3PiStudiesUtils::getD0Fit_PY(
		D0Fit_PNames pName,
		double Dst_D0Fit_D0_Kplus_PY,
		double Dst_D0Fit_D0_piplus_0_PY,
		double Dst_D0Fit_D0_piplus_1_PY,
		double Dst_D0Fit_D0_piplus_PY)
	{
		double py;
		switch (pName)
		{
		case D0Fit_PNames::D0Fit_D0_Kplus:
			py = Dst_D0Fit_D0_Kplus_PY;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_0:
			py = Dst_D0Fit_D0_piplus_0_PY;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_1:
			py = Dst_D0Fit_D0_piplus_1_PY;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus:
			py = Dst_D0Fit_D0_piplus_PY;
			break;
		default:
			throw InvalidDecayError("getD0Fit_PY: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return py;
	}

	double K3PiStudiesUtils::getD0Fit_PZ(
		D0Fit_PNames pName,
		double Dst_D0Fit_D0_Kplus_PZ,
		double Dst_D0Fit_D0_piplus_0_PZ,
		double Dst_D0Fit_D0_piplus_1_PZ,
		double Dst_D0Fit_D0_piplus_PZ)
	{
		double pz;
		switch (pName)
		{
		case D0Fit_PNames::D0Fit_D0_Kplus:
			pz = Dst_D0Fit_D0_Kplus_PZ;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_0:
			pz = Dst_D0Fit_D0_piplus_0_PZ;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus_1:
			pz = Dst_D0Fit_D0_piplus_1_PZ;
			break;
		case D0Fit_PNames::D0Fit_D0_piplus:
			pz = Dst_D0Fit_D0_piplus_PZ;
			break;
		default:
			throw InvalidDecayError("getD0Fit_PZ: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return pz;
	}

	double K3PiStudiesUtils::cTauMMToTauNS(double cTauMM)
	{
		double cTauM = cTauMM * _MM_TO_M;
		double tauSec = cTauM / _C_M_PER_SEC;
		double tauNS = tauSec * _SEC_TO_NS;

		return tauNS;
	}

	double K3PiStudiesUtils::tauNSToTauPS(double tauNS)
	{
		return tauNS * _NS_TO_PS;
	}

	bool K3PiStudiesUtils::isWithinDecayTimeBin(
		double dtime,
		const std::pair<double, double> &decayTimeLimits)
	{
		double lowerBinEdge = decayTimeLimits.first;
		double upperBinEdge = decayTimeLimits.second;

		return dtime >= lowerBinEdge && dtime < upperBinEdge;
	}

	unsigned int K3PiStudiesUtils::determineQuadrant(double sin2ThetaA, double sin2ThetaC)
	{
		unsigned int quadrant = 0;

		if (sin2ThetaA < 0 && sin2ThetaC < 0)
		{
			quadrant = 1;
		}
		else if (sin2ThetaA < 0 && sin2ThetaC > 0)
		{
			quadrant = 2;
		}
		else if (sin2ThetaA > 0 && sin2ThetaC < 0)
		{
			quadrant = 3;
		}
		else if (sin2ThetaA > 0 && sin2ThetaC > 0)
		{
			quadrant = 4;
		}

		return quadrant;
	}

	bool K3PiStudiesUtils::isKPi1LowerMassPair(
		const TLorentzVector &kminus_4vec,
		const TLorentzVector &piplus1_4vec,
		const TLorentzVector &piplus2_4vec)
	{
		TLorentzVector kpiplus1_4vec = kminus_4vec + piplus1_4vec;
		TLorentzVector kpiplus2_4vec = kminus_4vec + piplus2_4vec;

		double mkpi1 = kpiplus1_4vec.M();
		double mkpi2 = kpiplus2_4vec.M();

		return mkpi1 < mkpi2;
	}

	bool K3PiStudiesUtils::isKaonNeg(
		int kaonInd,
		int D0_P0_ID,
		int D0_P1_ID,
		int D0_P2_ID,
		int D0_P3_ID)
	{
		bool isKaonNeg;
		switch (kaonInd)
		{
		case 0:
			isKaonNeg = D0_P0_ID < 0;
			break;
		case 1:
			isKaonNeg = D0_P1_ID < 0;
			break;
		case 2:
			isKaonNeg = D0_P2_ID < 0;
			break;
		case 3:
			isKaonNeg = D0_P3_ID < 0;
			break;
		default:
			throw InvalidDecayError("isKaonNeg: Cannot find kaon with index " + std::to_string(kaonInd) + " in daughters.");
		}

		return isKaonNeg;
	}

	int K3PiStudiesUtils::findSSPion(
		bool kaonIsNeg,
		int D0_P0_ID,
		int D0_P1_ID,
		int D0_P2_ID,
		int D0_P3_ID)
	{
		std::vector<int> ssSignPionIndices;
		ssSignPionIndices.reserve(1);

		int ssPionID = kaonIsNeg ? -1 * _PION_ID : _PION_ID;
		if (D0_P0_ID == ssPionID)
		{
			ssSignPionIndices.push_back(0);
		}
		if (D0_P1_ID == ssPionID)
		{
			ssSignPionIndices.push_back(1);
		}
		if (D0_P2_ID == ssPionID)
		{
			ssSignPionIndices.push_back(2);
		}
		if (D0_P3_ID == ssPionID)
		{
			ssSignPionIndices.push_back(3);
		}

		if (ssSignPionIndices.size() != 1)
		{
			throw InvalidDecayError("findSSPion: Did not find same sign pion in daughters.");
		}

		return ssSignPionIndices[0];
	}

	std::vector<int> K3PiStudiesUtils::findOSPions(
		bool kaonIsNeg,
		int D0_P0_ID,
		int D0_P1_ID,
		int D0_P2_ID,
		int D0_P3_ID)
	{
		std::vector<int> osPionIndices;
		osPionIndices.reserve(2);

		int osPionID = kaonIsNeg ? _PION_ID : -1 * _PION_ID;
		if (D0_P0_ID == osPionID)
		{
			osPionIndices.push_back(0);
		}
		if (D0_P1_ID == osPionID)
		{
			osPionIndices.push_back(1);
		}
		if (D0_P2_ID == osPionID)
		{
			osPionIndices.push_back(2);
		}
		if (D0_P3_ID == osPionID)
		{
			osPionIndices.push_back(3);
		}

		if (osPionIndices.size() != 2)
		{
			throw InvalidDecayError("findOSPions: Did not find the two opposite sign pions in daughters.");
		}

		return osPionIndices;
	}

	D0Fit_PNames K3PiStudiesUtils::findD0FitSSPion(
		bool kaonIsNeg,
		int Dst_D0Fit_D0_Kplus_ID,
		int Dst_D0Fit_D0_piplus_0_ID,
		int Dst_D0Fit_D0_piplus_1_ID,
		int Dst_D0Fit_D0_piplus_ID)
	{
		int index = findSSPion(kaonIsNeg,
							   Dst_D0Fit_D0_Kplus_ID,
							   Dst_D0Fit_D0_piplus_0_ID,
							   Dst_D0Fit_D0_piplus_1_ID,
							   Dst_D0Fit_D0_piplus_ID);

		return indexToD0Fit_PName(index);
	}

	ReFit_PNames K3PiStudiesUtils::findReFitSSPion(
		bool kaonIsNeg,
		int Dst_ReFit_D0_Kplus_ID,
		int Dst_ReFit_D0_piplus_0_ID,
		int Dst_ReFit_D0_piplus_1_ID,
		int Dst_ReFit_D0_piplus_ID)
	{
		int index = findSSPion(kaonIsNeg,
							   Dst_ReFit_D0_Kplus_ID,
							   Dst_ReFit_D0_piplus_0_ID,
							   Dst_ReFit_D0_piplus_1_ID,
							   Dst_ReFit_D0_piplus_ID);

		return indexToReFit_PName(index);
	}

	std::vector<D0Fit_PNames> K3PiStudiesUtils::findD0FitOSPions(
		bool kaonIsNeg,
		int Dst_D0Fit_D0_Kplus_ID,
		int Dst_D0Fit_D0_piplus_0_ID,
		int Dst_D0Fit_D0_piplus_1_ID,
		int Dst_D0Fit_D0_piplus_ID)
	{
		std::vector<int> indices = findOSPions(
			kaonIsNeg,
			Dst_D0Fit_D0_Kplus_ID,
			Dst_D0Fit_D0_piplus_0_ID,
			Dst_D0Fit_D0_piplus_1_ID,
			Dst_D0Fit_D0_piplus_ID);

		std::vector<D0Fit_PNames> osPionNames;
		osPionNames.reserve(2);
		for (int i = 0; i < 2; i++)
		{
			osPionNames[i] = indexToD0Fit_PName(indices[i]);
		}

		return osPionNames;
	}

	std::vector<ReFit_PNames> K3PiStudiesUtils::findReFitOSPions(
		bool kaonIsNeg,
		int Dst_ReFit_D0_Kplus_ID,
		int Dst_ReFit_D0_piplus_0_ID,
		int Dst_ReFit_D0_piplus_1_ID,
		int Dst_ReFit_D0_piplus_ID)
	{
		std::vector<int> indices = findOSPions(
			kaonIsNeg,
			Dst_ReFit_D0_Kplus_ID,
			Dst_ReFit_D0_piplus_0_ID,
			Dst_ReFit_D0_piplus_1_ID,
			Dst_ReFit_D0_piplus_ID);

		std::vector<ReFit_PNames> osPionNames;
		osPionNames.reserve(2);
		for (int i = 0; i < 2; i++)
		{
			osPionNames[i] = indexToReFit_PName(indices[i]);
		}

		return osPionNames;
	}

	D0Fit_PNames K3PiStudiesUtils::indexToD0Fit_PName(int index)
	{
		D0Fit_PNames pName;
		switch (index)
		{
		case 0:
			pName = D0Fit_PNames::D0Fit_D0_Kplus;
			break;
		case 1:
			pName = D0Fit_PNames::D0Fit_D0_piplus_0;
			break;
		case 2:
			pName = D0Fit_PNames::D0Fit_D0_piplus_1;
			break;
		case 3:
			pName = D0Fit_PNames::D0Fit_D0_piplus;
			break;
		default:
			throw InvalidDecayError("indexToD0Fit_PName: Cannot find particle with index " + std::to_string(index) + " in daughters.");
		}

		return pName;
	}

	D0Fit_PNames K3PiStudiesUtils::findD0FitKaon(
		int Dst_D0Fit_D0_Kplus_ID,
		int Dst_D0Fit_D0_piplus_0_ID,
		int Dst_D0Fit_D0_piplus_1_ID,
		int Dst_D0Fit_D0_piplus_ID)
	{
		int index = findKaon(
			Dst_D0Fit_D0_Kplus_ID,
			Dst_D0Fit_D0_piplus_0_ID,
			Dst_D0Fit_D0_piplus_1_ID,
			Dst_D0Fit_D0_piplus_ID);

		return indexToD0Fit_PName(index);
	}

	ReFit_PNames K3PiStudiesUtils::indexToReFit_PName(int index)
	{
		ReFit_PNames pName;
		switch (index)
		{
		case 0:
			pName = ReFit_PNames::ReFit_D0_Kplus;
			break;
		case 1:
			pName = ReFit_PNames::ReFit_D0_piplus_0;
			break;
		case 2:
			pName = ReFit_PNames::ReFit_D0_piplus_1;
			break;
		case 3:
			pName = ReFit_PNames::ReFit_D0_piplus;
			break;
		default:
			throw InvalidDecayError("indexToReFit_PName: Cannot find particle with index " + std::to_string(index) + " in daughters.");
		}

		return pName;
	}

	ReFit_PNames K3PiStudiesUtils::findReFitKaon(
		int Dst_ReFit_D0_Kplus_ID,
		int Dst_ReFit_D0_piplus_0_ID,
		int Dst_ReFit_D0_piplus_1_ID,
		int Dst_ReFit_D0_piplus_ID)
	{
		int index = findKaon(
			Dst_ReFit_D0_Kplus_ID,
			Dst_ReFit_D0_piplus_0_ID,
			Dst_ReFit_D0_piplus_1_ID,
			Dst_ReFit_D0_piplus_ID);

		return indexToReFit_PName(index);
	}

	// returns 0 if P0 is kaon, 1 if P1 is, 2 if P2 is, 3 if P3 is
	int K3PiStudiesUtils::findKaon(
		int D0_P0_ID,
		int D0_P1_ID,
		int D0_P2_ID,
		int D0_P3_ID)
	{
		std::vector<int> kaonIndices;
		kaonIndices.reserve(1);

		if (std::abs(D0_P0_ID) == _KAON_ID)
		{
			kaonIndices.push_back(0);
		}
		if (std::abs(D0_P1_ID) == _KAON_ID)
		{
			kaonIndices.push_back(1);
		}
		if (std::abs(D0_P2_ID) == _KAON_ID)
		{
			kaonIndices.push_back(2);
		}
		if (std::abs(D0_P3_ID) == _KAON_ID)
		{
			kaonIndices.push_back(3);
		}

		if (kaonIndices.size() != 1)
		{
			throw InvalidDecayError("findKaon: Did not find kaon in daughters.");
		}

		return kaonIndices[0];
	}

	bool K3PiStudiesUtils::isD0(int dStarPiID)
	{
		return dStarPiID > 0;
	}

	bool K3PiStudiesUtils::isRS(bool isD0, bool isKaonNeg)
	{
		bool isRS;
		if (isD0)
		{
			isRS = isKaonNeg;
		}
		else
		{
			isRS = !isKaonNeg;
		}
		return isRS;
	}

	// param filesString: comma separated list of input root files
	// returns a vector of the input root files
	std::vector<std::string> K3PiStudiesUtils::buildListFromCommaSepStr(
		const std::string &filesString)
	{
		// make a copy so we don't modify the input string
		std::string filesStringCopy = filesString;

		// first remove any whitespace
		filesStringCopy.erase(
			std::remove_if(filesStringCopy.begin(), filesStringCopy.end(), ::isspace),
			filesStringCopy.end());

		// now separate on commas and build result vector
		std::stringstream ss(filesStringCopy);
		std::vector<std::string> result;
		while (ss.good())
		{
			std::string substr;
			getline(ss, substr, ',');
			result.push_back(substr);
		}

		return result;
	}

	std::string K3PiStudiesUtils::d0TimeBinToString(
		const std::pair<double, double> &decayTimeLimits,
		const std::string &unit)
	{
		double lowerBinEdge = decayTimeLimits.first;
		double upperBinEdge = decayTimeLimits.second;

		return std::to_string(lowerBinEdge) + " <= D0 decay t < " + std::to_string(upperBinEdge) + " [" + unit + "]";
	}

	std::vector<std::pair<double, double>> K3PiStudiesUtils::makeTimeBins(const std::vector<double> &upperBinEdges)
	{
		int numBins = upperBinEdges.size();
		std::vector<std::pair<double, double>> decayTimeLimits;
		decayTimeLimits.reserve(numBins + 1);

		for (int b = 0; b < numBins + 1; b++)
		{
			double lowerBinEdge = (b == 0) ? -std::numeric_limits<double>::infinity() : upperBinEdges[b - 1];
			double upperBinEdge = (b == numBins) ? std::numeric_limits<double>::infinity() : upperBinEdges[b]; // last bin is overflow
			decayTimeLimits.push_back(std::make_pair(lowerBinEdge, upperBinEdge));
		}

		return decayTimeLimits;
	}

	double K3PiStudiesUtils::getProbNNx(
		double D0_P0_ProbNNx,
		double D0_P1_ProbNNx,
		double D0_P2_ProbNNx,
		double D0_P3_ProbNNx,
		int ind)
	{
		double probNNx;
		switch (ind)
		{
		case 0:
			probNNx = D0_P0_ProbNNx;
			break;
		case 1:
			probNNx = D0_P1_ProbNNx;
			break;
		case 2:
			probNNx = D0_P2_ProbNNx;
			break;
		case 3:
			probNNx = D0_P3_ProbNNx;
			break;
		default:
			throw InvalidDecayError("getProbNNx: Cannot find particle with index " + std::to_string(ind) + " in daughters.");
		}
		return probNNx;
	}

	/**
	 * From John's apply_full_selection.py code
	 */
	double K3PiStudiesUtils::compute_delta_angle(
		double extra_px,
		double extra_py,
		double extra_pz,
		double d_px,
		double d_py,
		double d_pz)
	{
		TLorentzVector extra, d;
		extra.SetXYZM(extra_px, extra_py, extra_pz, _PION_MASS);
		d.SetXYZM(d_px, d_py, d_pz, _PION_MASS);
		return d.Angle(extra.Vect());
	}

	/**
	 * From John's apply_full_selection.py code
	 */
	double K3PiStudiesUtils::compute_delta_angle(
		double extra_px,
		double extra_py,
		double extra_pz,
		double extra_m,
		double d_px,
		double d_py,
		double d_pz,
		double d_m)
	{
		TLorentzVector extra, d;
		extra.SetXYZM(extra_px, extra_py, extra_pz, extra_m);
		d.SetXYZM(d_px, d_py, d_pz, d_m);
		return d.Angle(extra.Vect());
	}

	/**
	 * function to calculate helicity angle of soft pion
	 * From John's apply_full_selection.py code
	 */
	float K3PiStudiesUtils::helicity_angle_func(
		float d0_px,
		float d0_py,
		float d0_pz,
		float d0_m,
		const ROOT::RVec<float> &pis_px,
		const ROOT::RVec<float> &pis_py,
		const ROOT::RVec<float> &pis_pz,
		float pis_m)
	{
		// create the D* rest frame vector via a boost from lab to rest frame
		TLorentzVector d0_vec, pis_vec;
		d0_vec.SetXYZM(d0_px, d0_py, d0_pz, d0_m);
		pis_vec.SetXYZM(pis_px[0], pis_py[0], pis_pz[0], pis_m);
		auto dstar_lab_vec = (d0_vec + pis_vec);

		auto lab_n = dstar_lab_vec.Vect().Unit();
		auto pis_vec_boost = pis_vec;

		// boost soft pion to D* rest frame
		pis_vec_boost.Boost(-dstar_lab_vec.BoostVector()); // minus sign for lab -> CM
		auto pis_vec_boost_n = pis_vec_boost.Vect().Unit();
		// Helicity angle is angle of pi soft relative to D* lab frame momentum
		double helicity_angle = pis_vec_boost_n.Angle(lab_n);

		return helicity_angle;
	}

	/**
	 * function to calculate helicity angle of soft pion
	 * From John's apply_full_selection.py code
	 */
	float K3PiStudiesUtils::helicity_angle_func(
		float d0_px,
		float d0_py,
		float d0_pz,
		float d0_m,
		float pis_px,
		float pis_py,
		float pis_pz,
		float pis_m)
	{
		// create the D* rest frame vector via a boost from lab to rest frame
		TLorentzVector d0_vec, pis_vec;
		d0_vec.SetXYZM(d0_px, d0_py, d0_pz, d0_m);
		pis_vec.SetXYZM(pis_px, pis_py, pis_pz, pis_m);
		auto dstar_lab_vec = (d0_vec + pis_vec);

		auto lab_n = dstar_lab_vec.Vect().Unit();
		auto pis_vec_boost = pis_vec;

		// boost soft pion to D* rest frame
		pis_vec_boost.Boost(-dstar_lab_vec.BoostVector()); // minus sign for lab -> CM
		auto pis_vec_boost_n = pis_vec_boost.Vect().Unit();
		// Helicity angle is angle of pi soft relative to D* lab frame momentum
		double helicity_angle = pis_vec_boost_n.Angle(lab_n);

		return helicity_angle;
	}

	double K3PiStudiesUtils::getReFit_PE(
		ReFit_PNames pName,
		double Dst_ReFit_D0_Kplus_PE,
		double Dst_ReFit_D0_piplus_0_PE,
		double Dst_ReFit_D0_piplus_1_PE,
		double Dst_ReFit_D0_piplus_PE)
	{
		double pE;
		switch (pName)
		{
		case ReFit_PNames::ReFit_D0_Kplus:
			pE = Dst_ReFit_D0_Kplus_PE;
			break;
		case ReFit_PNames::ReFit_D0_piplus_0:
			pE = Dst_ReFit_D0_piplus_0_PE;
			break;
		case ReFit_PNames::ReFit_D0_piplus_1:
			pE = Dst_ReFit_D0_piplus_1_PE;
			break;
		case ReFit_PNames::ReFit_D0_piplus:
			pE = Dst_ReFit_D0_piplus_PE;
			break;
		default:
			throw InvalidDecayError("getReFit_PE: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return pE;
	}

	double K3PiStudiesUtils::getReFit_PX(
		ReFit_PNames pName,
		double Dst_ReFit_D0_Kplus_PX,
		double Dst_ReFit_D0_piplus_0_PX,
		double Dst_ReFit_D0_piplus_1_PX,
		double Dst_ReFit_D0_piplus_PX)
	{
		double px;
		switch (pName)
		{
		case ReFit_PNames::ReFit_D0_Kplus:
			px = Dst_ReFit_D0_Kplus_PX;
			break;
		case ReFit_PNames::ReFit_D0_piplus_0:
			px = Dst_ReFit_D0_piplus_0_PX;
			break;
		case ReFit_PNames::ReFit_D0_piplus_1:
			px = Dst_ReFit_D0_piplus_1_PX;
			break;
		case ReFit_PNames::ReFit_D0_piplus:
			px = Dst_ReFit_D0_piplus_PX;
			break;
		default:
			throw InvalidDecayError("getReFit_PX: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return px;
	}

	double K3PiStudiesUtils::getReFit_PY(
		ReFit_PNames pName,
		double Dst_ReFit_D0_Kplus_PY,
		double Dst_ReFit_D0_piplus_0_PY,
		double Dst_ReFit_D0_piplus_1_PY,
		double Dst_ReFit_D0_piplus_PY)
	{
		double py;
		switch (pName)
		{
		case ReFit_PNames::ReFit_D0_Kplus:
			py = Dst_ReFit_D0_Kplus_PY;
			break;
		case ReFit_PNames::ReFit_D0_piplus_0:
			py = Dst_ReFit_D0_piplus_0_PY;
			break;
		case ReFit_PNames::ReFit_D0_piplus_1:
			py = Dst_ReFit_D0_piplus_1_PY;
			break;
		case ReFit_PNames::ReFit_D0_piplus:
			py = Dst_ReFit_D0_piplus_PY;
			break;
		default:
			throw InvalidDecayError("getReFit_PY: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return py;
	}

	double K3PiStudiesUtils::getReFit_PZ(
		ReFit_PNames pName,
		double Dst_ReFit_D0_Kplus_PZ,
		double Dst_ReFit_D0_piplus_0_PZ,
		double Dst_ReFit_D0_piplus_1_PZ,
		double Dst_ReFit_D0_piplus_PZ)
	{
		double pz;
		switch (pName)
		{
		case ReFit_PNames::ReFit_D0_Kplus:
			pz = Dst_ReFit_D0_Kplus_PZ;
			break;
		case ReFit_PNames::ReFit_D0_piplus_0:
			pz = Dst_ReFit_D0_piplus_0_PZ;
			break;
		case ReFit_PNames::ReFit_D0_piplus_1:
			pz = Dst_ReFit_D0_piplus_1_PZ;
			break;
		case ReFit_PNames::ReFit_D0_piplus:
			pz = Dst_ReFit_D0_piplus_PZ;
			break;
		default:
			throw InvalidDecayError("getReFit_PZ: Cannot find daughter with name " + std::to_string(pName) + " in daughters.");
		}

		return pz;
	}

} // end namespace K3PiStudies