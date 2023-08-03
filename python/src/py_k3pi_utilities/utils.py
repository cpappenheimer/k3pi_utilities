import ROOT

def drawHist(h, statStyleStr, saveName):
    c = ROOT.TCanvas("c")

    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    h.Draw("HIST")
    c.SaveAs(saveName)


def test():
    print("Hello world")


def getAmpGenKName(IS_D0, IS_RS):
    if IS_D0:
        if IS_RS:
            return "K#"
        else:
            return "K~"
    else:
        if IS_RS:
            return "K~"
        else:
            return "K#"
        

def getAmpGenOSPiName(IS_D0, IS_RS):
    if IS_D0:
        if IS_RS:
            return "pi~"
        else:
            return "pi#"
    else:
        if IS_RS:
            return "pi#"
        else:
            return "pi~"


def getAmpGenSSPiName(IS_D0, IS_RS):
    if IS_D0:
        if IS_RS:
            return "pi#"
        else:
            return "pi~"
    else:
        if IS_RS:
            return "pi~"
        else:
            return "pi#"


def loadK3PiCUtils(buildDir, incDir):
    include = '#include "{}/K3PiStudiesUtils.h"'.format(incDir)
    lib = '{}/src/libK3PiStudiesUtils.so'.format(buildDir)

    # First include header, then load C++ library
    ROOT.gInterpreter.ProcessLine(include)
    ROOT.gSystem.Load(lib)


def rootFileToDF(tree, file):
    df = ROOT.RDataFrame(tree, file)
    return df


def csvFileToDF(file):
    df = ROOT.RDF.MakeCsvDataFrame(file, True)
    return df


def getKeysInROOTFile(file):
    print("Keys in {}:", file)
    f = ROOT.TFile.Open(file, 'read')
    for tkey in f.GetListOfKeys():
        key = tkey.GetName()
        print(key)
    f.Close()


def loadHistFromFile(file, histName):
    f = ROOT.TFile.Open(file, 'read')
    h = f.Get(histName)
    #print ("Name: {}, type: {}", histName, type(h))
    if h:
        print ("Found histogram named {} in file".format(histName))
        h.SetDirectory(0) # you have to do this so it doesn't get deleted when you close the file
    else:
        print ("Did not find histogram named {} in file".format(histName))
    f.Close()
    return h


def createHistFromDF(df, col, histName, histTitle, nBins, xMin, xMax):
    h = df.Histo1D((str(histName), str(histTitle), nBins, xMin, xMax), str(col))
    return h.GetValue()


def createHistWithSameXRange(df, getMyXRange, col, histSuffix):
    xMin = getMyXRange.GetXaxis().GetXmin()
    xMax = getMyXRange.GetXaxis().GetXmax() 
    numBins = getMyXRange.GetNbinsX()

    xTitle = getMyXRange.GetXaxis().GetTitle()
    yTitle = getMyXRange.GetYaxis().GetTitle()
    title = getMyXRange.GetTitle()
    titleStr = str( ROOT.K3PiStudies.K3PiStudiesUtils.makeTitleStr(title, xTitle, yTitle) )

    return createHistFromDF(df, str(col), col+histSuffix, titleStr, numBins, xMin, xMax)


def createLeg():
    leg = ROOT.TLegend(0.12, 0.76, 0.45, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.0)
    return leg


def drawSuperimposed(h1, h1LegTitle, h2, h2LegTitle, saveName):
    c = ROOT.TCanvas("c")

    ROOT.K3PiStudies.K3PiStudiesUtils.adjustYAxisForCompare(h1, h2)

    h1.SetLineColor(ROOT.kBlue)
    h1.SetLineWidth(2)
    h1.Draw("HIST")

    h2.SetLineColor(ROOT.kRed + 1)
    h2.SetLineWidth(2)
    h2.Draw("HIST SAME")

    leg = createLeg()
    leg.AddEntry(h1, h1LegTitle, "L")
    leg.AddEntry(h2, h2LegTitle, "L") 
    leg.Draw("SAME") 

    c.SaveAs(saveName)


def printAllCols(df):
    df.Display("").Print()


def printACol(df, colName):
    df.Display([colName]).Print()


def aliasAmpGen4VecComponents(df, pNum, ampGenName, noSymName, IS_D0, IS_RS):
    df = df.Alias( "_{}_{}_E".format(pNum, noSymName), "_{}_{}_E".format(pNum, ampGenName(IS_D0, IS_RS)) )
    df = df.Alias( "_{}_{}_Px".format(pNum, noSymName), "_{}_{}_Px".format(pNum, ampGenName(IS_D0, IS_RS)) )
    df = df.Alias( "_{}_{}_Py".format(pNum, noSymName), "_{}_{}_Py".format(pNum, ampGenName(IS_D0, IS_RS)) )
    df = df.Alias( "_{}_{}_Pz".format(pNum, noSymName), "_{}_{}_Pz".format(pNum, ampGenName(IS_D0, IS_RS)) )
    return df


def trimSpaceColNames(df):
    for c in df.GetColumnNames():
        cNoSpace = str(c).strip()
        if cNoSpace != str(c):
            df = df.Alias(cNoSpace, c)
    return df


def ampGentoTLorentzVector(ampGen4Vec_GeV):
    # AmpGen 4vec ordering: px, py, pz, pE
    return ROOT.K3PiStudies.K3PiStudiesUtils.toTLorentzVector(
        ampGen4Vec_GeV[3] * ROOT.K3PiStudies.K3PiStudiesUtils._GEV_TO_MEV,  # pE
        ampGen4Vec_GeV[0] * ROOT.K3PiStudies.K3PiStudiesUtils._GEV_TO_MEV,  # px
        ampGen4Vec_GeV[1] * ROOT.K3PiStudies.K3PiStudiesUtils._GEV_TO_MEV,  # py
        ampGen4Vec_GeV[2] * ROOT.K3PiStudies.K3PiStudiesUtils._GEV_TO_MEV) # pz