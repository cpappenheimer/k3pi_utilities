#include <sstream>

#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TCanvas.h>

#include <boost/algorithm/string.hpp>

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

	TLorentzVector K3PiStudiesUtils::toTLorentzVector(
		double pE,
		double px,
		double py,
		double pz)
	{
		TLorentzVector v;
		v.SetPxPyPzE(px, py, pz, pE);
		return v;
	}

	/**
	 * @see https://en.wikipedia.org/wiki/Inverse-variance_weighting
	 *
	 * @return pair where ans.first = weighted mean, ans.second = error on weighted mean
	 */
	std::pair<double, double> K3PiStudiesUtils::invVarWeightedAvg(
		const std::vector<double> &vals,
		const std::vector<double> &errs)
	{
		const unsigned int N = vals.size();
		if (errs.size() != N)
		{
			std::cout << "Error calculating weighted mean. Inputs have different sizes." << std::endl;
			return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
		}

		std::vector<double> weightingFactors;
		weightingFactors.reserve(N);
		for (auto const &err : errs)
		{
			weightingFactors.push_back(1.0 / (err * err));
		}

		double sumInvWeightingFactors = 0.0;
		for (auto const &w : weightingFactors)
		{
			sumInvWeightingFactors += w;
		}
		// std::cout << "sumInvWeightingFactors: " << sumInvWeightingFactors << std::endl;

		double weightedMean = 0.0;
		for (int i = 0; i < N; i++)
		{
			weightedMean += vals[i] * weightingFactors[i];
		}
		// std::cout << "weightedMean: " << weightedMean << std::endl;
		weightedMean /= sumInvWeightingFactors;
		// std::cout << "weightedMean after division: " << weightedMean << std::endl;

		return std::make_pair(weightedMean, sqrt(1.0 / sumInvWeightingFactors));
	}

	/**
	 * Eqs. 1 and 2 in Mike's angular distributions ANA note
	 *
	 * @return pair where asym.first = asymmetry, asym.second = error
	 */
	std::pair<double, double> K3PiStudiesUtils::calcAsymmetry(double nAbove, double nBelow)
	{
		double asym = (nAbove - nBelow) / (nAbove + nBelow);
		double asymErr = sqrt((1.0 - asym * asym) / (nAbove + nBelow));

		return std::make_pair(asym, asymErr);
	}

	/**
	 * @param countPositiveEntries true if want to count # positive (>= 0) entries, false if want to count # negative entries
	 * @return # of entries that are either positive or negative, depending on the value of `countPositiveEntries`, after the applyToEntry function is applied to each entry
	 */
	int K3PiStudiesUtils::countFuncResult(
		const TH1 *const h,
		std::function<double(double)> applyToEntry,
		bool countPositiveEntries)
	{
		unsigned int nBins = h->GetNbinsX();
		// std::cout << "Num bins: " << nBins << std::endl;

		/**
		 * From ROOT doc:
		 * bin = 0;       underflow bin
		 * bin = 1;       first bin with low-edge xlow INCLUDED
		 * bin = nbins;   last bin with upper-edge xup EXCLUDED
		 * bin = nbins+1; overflow bin
		 */
		unsigned int numUnderflow = h->GetBinContent(0);
		unsigned int numOverflow = h->GetBinContent(nBins + 1);
		if (numUnderflow != 0 || numOverflow != 0)
		{
			std::cout << "Underflow/overflow bins not empty. Cannot calculate number positive/negative entries accurately for " << h->GetName() << "." << std::endl;
			return -1;
		}

		unsigned int numNeg = 0;
		unsigned int numPos = 0;
		for (int b = 1; b <= nBins; b++)
		{
			double binCenter = h->GetXaxis()->GetBinCenter(b);
			double valToTest = applyToEntry(binCenter);
			if (valToTest >= 0.0)
			{
				numPos += h->GetBinContent(b);
			}
			else
			{
				numNeg += h->GetBinContent(b);
			}
		}

		unsigned int totEntries = h->GetEntries();
		if (numPos + numNeg != totEntries)
		{
			std::cout << "Error calculating # positive/# negative entries for " << h->GetName() << "." << std::endl;
			return -1;
		}

		if (countPositiveEntries)
		{
			return numPos;
		}
		else
		{
			return numNeg;
		}
	}

	void K3PiStudiesUtils::makeTLegendBkgTransparent(TLegend &leg)
	{
		leg.SetBorderSize(0);
		leg.SetFillColorAlpha(kWhite, 0.0);
	}

	void K3PiStudiesUtils::makeTPaveTextBkgTransparent(TPaveText &pt)
	{
		pt.SetFillColorAlpha(kWhite, 0.0);
	}

	std::string K3PiStudiesUtils::printRegionBoundsDeltaM(const std::string &regionName)
	{
		double upperBound = 0.0;
		double lowerBound = 0.0;
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			lowerBound = -1.0 * std::numeric_limits<double>::infinity();
			upperBound = std::numeric_limits<double>::infinity();
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			lowerBound = _SIG_REGION_LOW_DELTAM_BOUND_MEV;
			upperBound = _SIG_REGION_HIGH_DELTAM_BOUND_MEV;
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			lowerBound = std::numeric_limits<double>::quiet_NaN();
			upperBound = std::numeric_limits<double>::quiet_NaN();
		}

		return std::to_string(lowerBound) + " <= delta M <= " + std::to_string(upperBound) + " [MeV]";
	}

	std::string K3PiStudiesUtils::printRegionBoundsMD0(const std::string &regionName)
	{
		double upperBound = 0.0;
		double lowerBound = 0.0;
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			lowerBound = -1.0 * std::numeric_limits<double>::infinity();
			upperBound = std::numeric_limits<double>::infinity();
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			lowerBound = _SIG_REGION_LOW_MD0_BOUND_MEV;
			upperBound = _SIG_REGION_HIGH_MD0_BOUND_MEV;
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			lowerBound = std::numeric_limits<double>::quiet_NaN();
			upperBound = std::numeric_limits<double>::quiet_NaN();
		}

		return std::to_string(lowerBound) + " <= m(D0) <= " + std::to_string(upperBound) + " [MeV]";
	}

	/**
	 * @return a pair where pair.first = the lower limit to use on a delta m axis, pair.second = the upper limit to use on a delta m axis
	 */
	std::pair<double, double> K3PiStudiesUtils::getRegionAxisBoundsDeltaMMeV(const std::string &regionName)
	{
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			return std::make_pair(_ALL_REGS_DELTAM_AXIS_MIN_MEV, _ALL_REGS_DELTAM_AXIS_MAX_MEV);
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			return std::make_pair(_SIG_REGION_LOW_DELTAM_BOUND_MEV, _SIG_REGION_HIGH_DELTAM_BOUND_MEV);
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
		}
	}

	/**
	 * @return a pair where pair.first = the lower limit to use on a D0 mass axis, pair.second = the upper limit to use on a D0 mass axis
	 */
	std::pair<double, double> K3PiStudiesUtils::getRegionAxisBoundsMD0MeV(const std::string &regionName)
	{
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			return std::make_pair(_ALL_REGS_D0_MASS_AXIS_MIN_MEV, _ALL_REGS_D0_MASS_AXIS_MAX_MEV);
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			return std::make_pair(_SIG_REGION_LOW_MD0_BOUND_MEV, _SIG_REGION_HIGH_MD0_BOUND_MEV);
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			return std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
		}
	}

	bool K3PiStudiesUtils::isInDeltaMRegion(
		const std::string &regionName,
		double deltaMMeV)
	{
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			return true;
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			return deltaMMeV >= _SIG_REGION_LOW_DELTAM_BOUND_MEV && deltaMMeV <= _SIG_REGION_HIGH_DELTAM_BOUND_MEV;
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			return false;
		}
	}

	bool K3PiStudiesUtils::isInD0MassRegion(
		const std::string &regionName,
		double d0MassMeV)
	{
		if (boost::iequals(regionName, _ALL_REGION_FLAG))
		{
			return true;
		}
		else if (boost::iequals(regionName, _SIG_REGION_FLAG))
		{
			return d0MassMeV >= _SIG_REGION_LOW_MD0_BOUND_MEV && d0MassMeV <= _SIG_REGION_HIGH_MD0_BOUND_MEV;
		}
		else
		{
			std::cout << "Unknown region " << regionName << "!" << std::endl;
			return false;
		}
	}

	void K3PiStudiesUtils::adjustYAxisForCompare(
		TH1 *const h1,
		TH1 *const h2)
	{
		const double yMax1 = h1->GetMaximum();
		const double yMax2 = h2->GetMaximum();
		const double overallMax = yMax1 > yMax2 ? yMax1 : yMax2;
		const double yMax = 1.1 * overallMax;
		h1->GetYaxis()->SetRangeUser(0.0, yMax);
		h2->GetYaxis()->SetRangeUser(0.0, yMax);
	}

	void K3PiStudiesUtils::makeNormalizedComparisonPlot(
		TH1 *const h1,
		TH1 *const h2,
		const TString &legLine1,
		const TString &legLine2,
		bool addNumEntries,
		const TString &saveName,
		const TString& unit,
		bool updateYLabel)
	{
		const unsigned int n1 = h1->GetEntries();
		const unsigned int n2 = h2->GetEntries();
		if (n1 == 0 || n2 == 0)
		{
			std::cout << "WARNING: For " << saveName << ", h1 or h2 has 0 entries. Cannot make comparison histogram!" << std::endl;
			return;
		}

		TCanvas c1;

		if (updateYLabel)
		{
			double xMin = h1->GetXaxis()->GetXmin(); 
    		double xMax = h1->GetXaxis()->GetXmax();
    		unsigned int numBins = h1->GetNbinsX();
			TString updatedYLabel = makeYAxisLabel(numBins, xMin, xMax, unit, true);
			h1->SetYTitle(updatedYLabel);
			h2->SetYTitle(updatedYLabel);
		}

		h1->Scale(1.0 / h1->Integral());
		h2->Scale(1.0 / h2->Integral());

		adjustYAxisForCompare(h1, h2);

		h1->SetLineColor(kBlue);
		h1->SetLineWidth(2);
		h1->Draw("HIST");

		h2->SetLineColor(kRed + 1);
		h2->SetLineWidth(2);
		h2->Draw("HIST SAME");

		TLegend leg(0.12, 0.76, 0.32, 0.89);
		makeTLegendBkgTransparent(leg);
		leg.AddEntry(h1, legLine1, "L");
		leg.AddEntry(h2, legLine2, "L");
		leg.Draw("SAME");

		if (addNumEntries)
		{
			TString e1 = std::to_string(n1);
			TString e2 = std::to_string(n2);

			TPaveText pt1(0.60, 0.8, 0.9, 0.9, "NDC"); // NDC sets coords
			makeTPaveTextBkgTransparent(pt1);
			pt1.AddText("n(" + legLine1 + ") = " + e1);
			pt1.AddText("n(" + legLine2 + ") = " + e2);
			pt1.Draw("SAME");

			c1.SaveAs(saveName);
		}
		else
		{
			c1.SaveAs(saveName);
		}
	}

	/**
	 * See Eq. 42 in Kutschke's An Angular Distribution Cookbook
	 * @return angle between the (4,5) decay plane and the (6,7) decay plane in mother rest frame, ranging from -pi to pi
	 */
	double K3PiStudiesUtils::angleBetweenDecayPlanesKutschke(
		const TVector3 &d4_motherRestFrame,
		const TVector3 &d5_motherRestFrame,
		const TVector3 &d6_motherRestFrame,
		const TVector3 &d7_motherRestFrame)
	{
		TVector3 nPrime = (d4_motherRestFrame.Unit()).Cross(d5_motherRestFrame.Unit());
		TVector3 nHatPrime = nPrime.Unit();

		TVector3 nDoublePrime = (d6_motherRestFrame.Unit()).Cross(d7_motherRestFrame.Unit());
		TVector3 nHatDoublePrime = nDoublePrime.Unit();

		TVector3 p2Hat = (d4_motherRestFrame + d5_motherRestFrame).Unit();

		double cosPhi = nHatDoublePrime.Dot(nHatPrime);
		double sinPhi = (nHatDoublePrime.Cross(nHatPrime)).Dot(p2Hat);

		return TMath::ATan2(sinPhi, cosPhi);
	}

	TString K3PiStudiesUtils::makeTitleStr(
		const TString &title,
		const TString &xLabel,
		const TString &yLabel)
	{
		return title + ";" + xLabel + ";" + yLabel;
	}

	TString K3PiStudiesUtils::makeYAxisLabel(
		int numBins,
		double axisMin,
		double axisMax,
		const TString &unit,
		bool normalizedPlot)
	{
		double axisLength = axisMax - axisMin;
		double binSize = axisLength / numBins;

		TString yType = (normalizedPlot) ? "Fraction" : "Events";

		return yType + " / " + std::to_string(binSize) + " " + unit;
	}

	void K3PiStudiesUtils::changeToRainbowPalette()
	{
		gStyle->SetPalette(kRainBow);
	}

	/**
	 * @param v1v2AngleIsNegPiToPi true if `v1v2Angle` ranges from -pi to pi, false if it ranges from 0 to 2 pi
	 * @return diff between .Angle method and our angle calculation (`v1v2Angle`)
	 */
	double K3PiStudiesUtils::verifyAngle(
		const TVector3 &v1,
		const TVector3 &v2,
		double v1v2Angle,
		bool v1v2AngleIsNegPiToPi,
		const std::string &angleName,
		bool printDiff)
	{
		double angleDiff = 0.0;

		double angleToCompare = v1.Angle(v2);
		// .Angle uses acos, which ranges from 0 to pi. Need to get correct quadrant and put it in range -pi to pi
		if (sin(v1v2Angle) < 0.0)
		{
			angleToCompare *= -1.0;
		}

		if (!v1v2AngleIsNegPiToPi)
		{
			angleToCompare = changeAngleRange_0_to_2pi(angleToCompare);
		}

		// std::cout << angleName << ": " << v1v2Angle << std::endl;
		// std::cout << "From .Angle(): " << angleToCompare << std::endl;
		// std::cout << angleName << " (deg) : " << radToDeg(v1v2Angle) << std::endl;
		// std::cout << "From .Angle() (deg) : " << radToDeg(angleToCompare) << std::endl;

		areDoublesEqual(combinedToleranceCompare, v1v2Angle, angleToCompare, angleName + " / .Angle()", printDiff);
		angleDiff = std::fabs(v1v2Angle - angleToCompare);

		return angleDiff;
	}

	bool K3PiStudiesUtils::isExactlyEqual(double d1, double d2)
	{
		return d1 == d2;
	}

	/**
	 * @see https://stackoverflow.com/a/15012792
	 */
	bool K3PiStudiesUtils::combinedToleranceCompare(double x, double y)
	{
		double maxXYOne = std::max({1.0, std::fabs(x), std::fabs(y)});

		// std::cout << "maxXYOne: " << maxXYOne << std::endl;
		// std::cout << "EPS: " << std::numeric_limits<double>::epsilon() * maxXYOne << std::endl;
		// std::cout << "diff: " << std::fabs(x - y) << std::endl;

		return std::fabs(x - y) <= _COMPARE_EPS * maxXYOne;
	}

	void K3PiStudiesUtils::silenceROOTHistSaveMsgs()
	{
		gErrorIgnoreLevel = kWarning;
	}

	double K3PiStudiesUtils::radToDeg(double angleRad)
	{
		return TMath::RadToDeg() * angleRad;
	}

	/**
	 * @param angle_0_to_2pi angle in range 0 to 2pi
	 * @return angle in range -pi to pi
	 */
	double K3PiStudiesUtils::changeAngleRange_neg_pi_to_pi(double angle_0_to_2pi)
	{
		if (angle_0_to_2pi > _PI)
		{
			return angle_0_to_2pi - 2.0 * _PI;
		}
		else
		{
			return angle_0_to_2pi;
		}
	}

	/**
	 * @param angle_neg_pi_to_pi angle that ranges from -pi to pi
	 * @return angle that ranges from 0 to 2pi
	 */
	double K3PiStudiesUtils::changeAngleRange_0_to_2pi(double angle_neg_pi_to_pi)
	{
		if (angle_neg_pi_to_pi < 0.0)
		{
			// std::cout << "Neg phi" << std::endl;
			return angle_neg_pi_to_pi + 2.0 * K3PiStudiesUtils::_PI;
		}
		else
		{
			return angle_neg_pi_to_pi;
		}
	}

	/**
	 * @param isEqualFunc returns true if the d1, d2 are equal; false otherwise
	 */
	bool K3PiStudiesUtils::areDoublesEqual(
		std::function<bool(double, double)> isEqualFunc,
		double d1,
		double d2,
		const std::string &varName,
		bool printDiff)
	{
		if (!isEqualFunc(d1, d2))
		{
			if (printDiff)
			{
				std::cout << "Found difference for " << varName << ": " << d1 << ", " << d2 << "; diff = " << d1 - d2 << std::endl;
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

	/**
	 * @param pA_IN_D0CM K
	 * @param pB_IN_D0CM OS pi 1
	 * @param pC_IN_D0CM SS pi
	 * @param pD_IN_D0CM OS pi 2
	 * @returns vector with entries m12, m34, cos12, cos34, phi
	 */
	std::vector<double> K3PiStudiesUtils::calc_phsp(
		const TLorentzVector &pD0_IN_D0CM,
		const TLorentzVector &pA_IN_D0CM, // K-
		const TLorentzVector &pB_IN_D0CM, // OS pi 1
		const TLorentzVector &pC_IN_D0CM, // SS pi
		const TLorentzVector &pD_IN_D0CM) // OS pi 2
	{
		//  note that _pA_IN_D0CM_MEV, _pB_IN_D0CM_MEV, etc., are in the D0 CM.
		//  we are going to define zhat as the pA_3vec+pB_3vec direction.
		//  to consider the helicity angles of the AB and CD pairs in their
		//  respective CMs, we should make Lorentz transformations along
		//  the zhat (or -zhat) directions. Note that the CD system is moving
		//  along the -zhat direction to start.
		const TLorentzVector pAB_4vec = pA_IN_D0CM + pB_IN_D0CM;
		const double mAB = pAB_4vec.M(); // m12

		const TLorentzVector pCD_4vec = pC_IN_D0CM + pD_IN_D0CM;
		const double mCD = pCD_4vec.M(); // m34

		const TVector3 pA_3vec = pA_IN_D0CM.Vect();
		const TVector3 pB_3vec = pB_IN_D0CM.Vect();
		const TVector3 pC_3vec = pC_IN_D0CM.Vect();
		const TVector3 pD_3vec = pD_IN_D0CM.Vect();
		const TVector3 pAB_3vec = pAB_4vec.Vect();

		const TVector3 yhat = (pA_3vec.Cross(pB_3vec)).Unit();
		const TVector3 yhatPrime = (pC_3vec.Cross(pD_3vec)).Unit();
		const TVector3 zhat = pAB_3vec.Unit();
		const TVector3 xhat = (yhat.Cross(zhat)).Unit();

		const double cosPhi = (yhat.Dot(yhatPrime));
		const double sinPhi = (xhat.Dot(yhatPrime));
		const double phi = changeAngleRange_0_to_2pi(TMath::ATan2(sinPhi, cosPhi));

		const double energyAB = pAB_4vec.E();
		TVector3 betaAB;
		betaAB.SetXYZ(pAB_4vec.Px() / energyAB, pAB_4vec.Py() / energyAB, pAB_4vec.Pz() / energyAB);
		const double energyCD = pCD_4vec.E();
		TVector3 betaCD;
		betaCD.SetXYZ(pCD_4vec.Px() / energyCD, pCD_4vec.Py() / energyCD, pCD_4vec.Pz() / energyCD);

		TLorentzVector pAprime_4vec = pA_IN_D0CM;
		pAprime_4vec.Boost(-1. * betaAB);
		const TVector3 pAprime_3vec = pAprime_4vec.Vect();
		const double paPrimeZ = pAprime_3vec.Dot(zhat);
		const double paPrimeMag = pAprime_3vec.Mag();

		TLorentzVector pCprime_4vec = pC_IN_D0CM;
		pCprime_4vec.Boost(-1. * betaCD);
		const TVector3 pCprime_3vec = pCprime_4vec.Vect();
		const double pcPrimeZ = pCprime_3vec.Dot(zhat);
		const double pcPrimeMag = pCprime_3vec.Mag();

		const double cosThetaA = paPrimeZ / paPrimeMag; // cos theta 12
		const double cosThetaC = pcPrimeZ / pcPrimeMag; // cos theta 34

		std::vector<double> vars = {mAB, mCD, cosThetaA, cosThetaC, phi};
		return vars;
	}

	/*
	 * Function to calculate phase space from John's apply_full_selection.py code
	 * returns vector with entries: {m12, m34, cos1, cos2, phi, m13, phiAngleDiff}
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
		bool pi1GoesWithK,
		bool verifyAngles,
		bool printDiff)
	{
		TLorentzVector d1_piGoesWithPi, d2_ssPi, d3_k, d4_piGoesWithK;
		d2_ssPi.SetPtEtaPhiM(Pi_SS_D0Fit_PT, Pi_SS_D0Fit_ETA, Pi_SS_D0Fit_PHI, K3PiStudiesUtils::_PION_MASS);
		d3_k.SetPtEtaPhiM(K_D0Fit_PT, K_D0Fit_ETA, K_D0Fit_PHI, K3PiStudiesUtils::_KAON_MASS);

		// figure out which pi to associate with k
		if (pi1GoesWithK)
		{ // case where kpi1 goes with k
			d1_piGoesWithPi.SetPtEtaPhiM(Pi_OS2_D0Fit_PT, Pi_OS2_D0Fit_ETA, Pi_OS2_D0Fit_PHI, K3PiStudiesUtils::_PION_MASS);
			d4_piGoesWithK.SetPtEtaPhiM(Pi_OS1_D0Fit_PT, Pi_OS1_D0Fit_ETA, Pi_OS1_D0Fit_PHI, K3PiStudiesUtils::_PION_MASS);
		}
		else
		{ // case where kpi2 goes with k
			d1_piGoesWithPi.SetPtEtaPhiM(Pi_OS1_D0Fit_PT, Pi_OS1_D0Fit_ETA, Pi_OS1_D0Fit_PHI, K3PiStudiesUtils::_PION_MASS);
			d4_piGoesWithK.SetPtEtaPhiM(Pi_OS2_D0Fit_PT, Pi_OS2_D0Fit_ETA, Pi_OS2_D0Fit_PHI, K3PiStudiesUtils::_PION_MASS);
		}

		// Boost everything to D0 restframe
		auto mum = d1_piGoesWithPi + d2_ssPi + d3_k + d4_piGoesWithK;
		double m12 = (d1_piGoesWithPi + d2_ssPi).M();
		double m34 = (d3_k + d4_piGoesWithK).M();
		double m13 = (d1_piGoesWithPi + d3_k).M();
		d1_piGoesWithPi.Boost(-mum.BoostVector());
		d2_ssPi.Boost(-mum.BoostVector());
		d3_k.Boost(-mum.BoostVector());
		d4_piGoesWithK.Boost(-mum.BoostVector());

		TLorentzVector d1_piGoesWithPi2, d3_k4;
		d1_piGoesWithPi2 = d1_piGoesWithPi + d2_ssPi;
		d3_k4 = d3_k + d4_piGoesWithK;

		TVector3 d1_piGoesWithPin = d1_piGoesWithPi.Vect().Unit();
		TVector3 d2_ssPin = d2_ssPi.Vect().Unit();
		TVector3 d3_kn = d3_k.Vect().Unit();
		TVector3 d4_piGoesWithKn = d4_piGoesWithK.Vect().Unit();
		TVector3 d1_piGoesWithPi2n = d1_piGoesWithPi2.Vect().Unit();
		TVector3 d3_k4n = d3_k4.Vect().Unit();

		TVector3 n1 = d1_piGoesWithPin.Cross(d2_ssPin);
		TVector3 n2 = d3_kn.Cross(d4_piGoesWithKn);
		TVector3 n3 = n1.Unit().Cross(n2.Unit());

		// Calculation of the angle Phi between the planes, in range -pi to pi
		double cosp = n1.Unit().Dot(n2.Unit());
		double sinp = n3.Dot(d3_k4n);
		double phi = TMath::ATan2(sinp, cosp);

		double phiDiff = (verifyAngles) ? K3PiStudiesUtils::verifyAngle(n1.Unit(), n2.Unit(), phi, true, "phi", false) : 0.0;

		// Vectors in rest fram of their resonance.
		TLorentzVector d1_piGoesWithPir = d1_piGoesWithPi;
		TLorentzVector d3_kr = d3_k;
		d1_piGoesWithPir.Boost(-d1_piGoesWithPi2.BoostVector());
		d3_kr.Boost(-d3_k4.BoostVector());
		TVector3 d1_piGoesWithPirn = d1_piGoesWithPir.Vect().Unit();
		TVector3 d3_krn = d3_kr.Vect().Unit();

		// helicity angle for d1_piGoesWithPi2 and d3_k4 frame
		double cos1 = d1_piGoesWithPi2n.Dot(d1_piGoesWithPirn);
		double cos2 = d3_k4n.Dot(d3_krn);

		std::vector<double> vars = {m12, m34, cos1, cos2, phi, m13, phiDiff};
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