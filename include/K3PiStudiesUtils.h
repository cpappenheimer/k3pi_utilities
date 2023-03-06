#pragma once

#include <string>
#include <vector>
#include <utility>

#include <TMath.h>
#include <TLorentzVector.h>

namespace K3PiStudies {

struct InvalidDecayError : std::exception 
{
   const std::string _msg;

   InvalidDecayError(const std::string& msg) : _msg(msg)
   {
   }

   const char* what() const noexcept override
   {
      return _msg.c_str();
   }
}; // end InvalidDecayError

class K3PiStudiesUtils final {
   public:
      static constexpr double _C_M_PER_SEC = 3.0 * pow(10, 8);
      static constexpr double _MM_TO_M = 1.0 / 1000.0;
      static constexpr double _SEC_TO_NS = pow(10, 9);
      static const unsigned int _KAON_ID = 321;
      static const unsigned int _PION_ID = 211;
      static constexpr double _PI = TMath::Pi();
      static constexpr double _GEV_TO_MEV = 1000.0;
      static constexpr double _D0_LIFETIME_PS = 0.410;
      static constexpr double _NS_TO_PS = 1000.0;
      static const std::string _RS_FLAG;
      static const std::string _WS_FLAG;
      static const std::string _BOTH_FLAG;

      K3PiStudiesUtils() = default;
      K3PiStudiesUtils(const K3PiStudiesUtils& copyMe) = default;
      K3PiStudiesUtils(K3PiStudiesUtils&& moveMe) = default;
      ~K3PiStudiesUtils() = default;
      K3PiStudiesUtils& operator=(const K3PiStudiesUtils& copyMe) = default;
      K3PiStudiesUtils& operator=(K3PiStudiesUtils&& moveMe) = default;

      static bool isWithinDecayTimeBin(
		   double dtime,
			const std::pair<double, double>& decayTimeLimits);

      static unsigned int determineQuadrant(double sin2ThetaA, double sin2ThetaC);

      static bool doesPi1GoWithK(
         const TLorentzVector& kminus_4vec, 
         const TLorentzVector& piplus1_4vec, 
         const TLorentzVector& piplus2_4vec);

      static bool isKaonNeg(
         int kaonInd,
         int D0_P0_ID, 
		   int D0_P1_ID, 
		   int D0_P2_ID,
		   int D0_P3_ID);

      static int findOppSignPion(
		   bool kaonIsNeg,
		   int D0_P0_ID, 
		   int D0_P1_ID, 
		   int D0_P2_ID,
		   int D0_P3_ID);

      static std::vector<int> findSSPions(
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
		   const std::pair<double, double>& decayTimeLimits,
		   const std::string& unit);

      static std::vector<std::pair<double, double>> makeTimeBins(const std::vector<double>& upperBinEdges);
}; // end K3PiStudiesUtils class

} // end namespace K3PiStudies