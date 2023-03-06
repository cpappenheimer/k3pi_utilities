#include <sstream>

#include <TLorentzVector.h>

#include "K3PiStudiesUtils.h"


namespace K3PiStudies {

	const std::string K3PiStudiesUtils::_RS_FLAG = "RS";
    const std::string K3PiStudiesUtils::_WS_FLAG = "WS";
    const std::string K3PiStudiesUtils::_BOTH_FLAG = "BOTH";
	const double K3PiStudiesUtils::_C_M_PER_SEC = 3.0 * TMath::Power(10, 8);
	const double K3PiStudiesUtils::_SEC_TO_NS = TMath::Power(10, 9);

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
				  const std::pair<double, double>& decayTimeLimits)
	{
		double lowerBinEdge = decayTimeLimits.first;
		double upperBinEdge = decayTimeLimits.second;
	
		return dtime >= lowerBinEdge && dtime < upperBinEdge;
	}


    unsigned int K3PiStudiesUtils::determineQuadrant(double sin2ThetaA, double sin2ThetaC)
	{
		unsigned int quadrant = 0;
		
		if (sin2ThetaA < 0 && sin2ThetaC < 0) {
			quadrant = 1;
		} else if (sin2ThetaA < 0 && sin2ThetaC > 0) {
			quadrant = 2;
		} else if (sin2ThetaA > 0 && sin2ThetaC < 0) {
			quadrant = 3;
		} else if (sin2ThetaA > 0 && sin2ThetaC > 0) {
			quadrant = 4;
		}
		
		return quadrant;
	}


    bool K3PiStudiesUtils::doesPi1GoWithK(
		const TLorentzVector& kminus_4vec, 
		const TLorentzVector& piplus1_4vec, 
		const TLorentzVector& piplus2_4vec)
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
		switch(kaonInd) 
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


    int K3PiStudiesUtils::findOppSignPion(
		bool kaonIsNeg,
		int D0_P0_ID, 
		int D0_P1_ID, 
		int D0_P2_ID,
		int D0_P3_ID)
	{
		std::vector<int> oppSignPionIndices;
		oppSignPionIndices.reserve(1);

		int oppPionID = kaonIsNeg ? -1*_PION_ID : _PION_ID;
		if (D0_P0_ID == oppPionID)
		{
			oppSignPionIndices.push_back(0);
		}
		if (D0_P1_ID == oppPionID)
		{
			oppSignPionIndices.push_back(1);
		}
		if (D0_P2_ID == oppPionID)
		{
			oppSignPionIndices.push_back(2);
		}
		if (D0_P3_ID == oppPionID)
		{
			oppSignPionIndices.push_back(3);
		}

		if (oppSignPionIndices.size() != 1)
		{
			throw InvalidDecayError("findOppSignPion: Did not find opposite sign pion in daughters.");
		}

		return oppSignPionIndices[0];
	}


    std::vector<int> K3PiStudiesUtils::findSSPions(
		bool kaonIsNeg,
		int D0_P0_ID, 
		int D0_P1_ID, 
		int D0_P2_ID,
		int D0_P3_ID)
	{
		std::vector<int> sameSignPionIndices;
		sameSignPionIndices.reserve(2);
		
		int ssPionID = kaonIsNeg ? _PION_ID : -1*_PION_ID;
		if (D0_P0_ID == ssPionID)
		{
			sameSignPionIndices.push_back(0);
		}
		if (D0_P1_ID == ssPionID)
		{
			sameSignPionIndices.push_back(1);
		}
		if (D0_P2_ID == ssPionID)
		{
			sameSignPionIndices.push_back(2);
		}
		if (D0_P3_ID == ssPionID)
		{
			sameSignPionIndices.push_back(3);
		}

		if (sameSignPionIndices.size() != 2)
		{
			throw InvalidDecayError("findSSPions: Did not find two same sign pions in daughters.");
		}

		return sameSignPionIndices;
	}


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
 		while (ss.good()) {
    		std::string substr;
    		getline(ss, substr, ',');
   			result.push_back(substr);
  		}

  		return result;
	}


	std::string K3PiStudiesUtils::d0TimeBinToString(
		const std::pair<double, double>& decayTimeLimits,
		const std::string& unit)
	{
		double lowerBinEdge = decayTimeLimits.first;
		double upperBinEdge = decayTimeLimits.second;

		return std::to_string(lowerBinEdge) + " <= D0 decay t < " 
			+ std::to_string(upperBinEdge) + " [" + unit + "]";
	}


	std::vector<std::pair<double, double>> K3PiStudiesUtils::makeTimeBins(const std::vector<double>& upperBinEdges)
	{
		int numBins = upperBinEdges.size();
		std::vector<std::pair<double, double>> decayTimeLimits;
		decayTimeLimits.reserve(numBins+1);

		for (int b=0; b<numBins+1; b++)
		{
			double lowerBinEdge = (b==0) ? -std::numeric_limits<double>::infinity() : upperBinEdges[b-1];
			double upperBinEdge = (b==numBins) ? std::numeric_limits<double>::infinity() : upperBinEdges[b]; // last bin is overflow
			decayTimeLimits.push_back(std::make_pair(lowerBinEdge, upperBinEdge));
		}

		return decayTimeLimits;
	}

} // end namespace K3PiStudies