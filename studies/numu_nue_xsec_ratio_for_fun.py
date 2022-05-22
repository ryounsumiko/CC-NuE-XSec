import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
from tools import PlotTools,Utilities
from array import array
ROOT.TH1.AddDirectory(False)

DrawMC = [
    #"MnvTune-v3",
    "MnvTune-v1.2",
    "NuWro-SF",
    #"NuWro-LFG",
    "GENIE3-10a",
    #"GENIE3-10b",
    #"GENIE3-02a",
]

def xsec_ratio(mnvplotter,he,hmu,hemc,hmumc):
    he.AddMissingErrorBandsAndFillWithCV(hmu)
    hmu.AddMissingErrorBandsAndFillWithCV(he)

    he.Divide(he,hmu)
    hemc.Divide(hemc,hmumc)
    mnvplotter.DrawDataMCWithErrorBand(he.GetCVHistoWithError(),hemc.GetCVHistoWithStatError(),1.0,"TR")

def ReadLEResult(filePath):
    xsec = []
    error_matrix = ROOT.TMatrixD(67,67)
    j=0
    with open(filePath) as f:
        lines = f.readlines()
        for line in lines:
            l = line.split()
            if len(l)==2:
                #cross section value
                xsec.append(float(l[1]))
            else:
                #error
                #print(j,len(l))
                for i,v in enumerate(l):
                    error_matrix[j][i]=float(v)
                j+=1
    cov = ROOT.TMatrixD(error_matrix,ROOT.TMatrixD.kTransposeMult,error_matrix)
    #print(cov[0][1])
    #assert(cov[3][2]==cov[2][3])
    return FillHistogram(xsec,cov)

def FillHistogram(xsec,cov):
    q3_bin = [0,0.2,0.3,0.4,0.5,0.6,0.8]
    q3_bin_target = [0.2*i for i in range(5)]
    Eavail_bin = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.8]
    bounds = [    3,   10,   22,   36,   51 ]
    mapping = {}
    x=y=1
    hist = PlotUtils.MnvH2D("Eavail_q3_LE",";E_{avail} (GeV);q_{3} (GeV);(10^{-38})",len(Eavail_bin)-1,array('d',Eavail_bin),len(q3_bin)-1,array('d',q3_bin))
    for i,v in enumerate(xsec):
        if i in bounds:
            x=1
            y+=1
        hist.SetBinContent(x,y,v)
        mapping[(x,y)]=i
        hist.SetBinError(x,y,1e-4*math.sqrt(cov[i][i]))
        x+=1

    hist_target = PlotUtils.MnvH2D("Eavail_q3_LErebin",";E_{avail} (GeV);q_{3} (GeV);(10^{-39})",len(Eavail_bin)-1,array('d',Eavail_bin),len(q3_bin_target)-1,array('d',q3_bin_target))
    for x in range(1,len(Eavail_bin)):
        if (x,1) in mapping:
            i2 = mapping[(x,1)]
            hist_target.SetBinContent(x,1,xsec[i2])
            hist_target.SetBinError(x,1,1e-4*math.sqrt(cov[i2][i2]))
        i1 = mapping[(x,2)] if (x,2) in mapping else None
        i2 = mapping[(x,3)] if (x,3) in mapping else None
        if i1 and i2:
            v = (xsec[i1]+xsec[i2])/2
            err = math.sqrt(sum(cov[a][b] for a in [i1,i2] for b in [i1,i2]))/2
            hist_target.SetBinContent(x,2,v)
            hist_target.SetBinError(x,2,1e-4*err)
        elif i2:
            hist_target.SetBinContent(x,2,xsec[i2]/2)
            hist_target.SetBinError(x,2,1e-4*math.sqrt(cov[i2][i2])/2)
        i1 = mapping[(x,4)] if (x,4) in mapping else None
        i2 = mapping[(x,5)] if (x,5) in mapping else None
        if i1 and i2:
            v = (xsec[i1]+xsec[i2])/2
            err = math.sqrt(sum(cov[a][b] for a in [i1,i2] for b in [i1,i2]))/2
            hist_target.SetBinContent(x,3,v)
            hist_target.SetBinError(x,3,1e-4*err)
        elif i2:
            hist_target.SetBinContent(x,3,xsec[i2]/2)
            hist_target.SetBinError(x,3,1e-4*math.sqrt(cov[i2][i2])/2)
        if (x,6) in mapping:
            i2 = mapping[(x,6)]
            hist_target.SetBinContent(x,4,xsec[i2])
            hist_target.SetBinError(x,4,1e-4*math.sqrt(cov[i2][i2]))
    return hist,hist_target

def RebinHistogram(h,cov):
    q3_bin_target = [0.2*i for i in range(5)]
    Eavail_bin = [0,0.04,0.08,0.12,0.16,0.24,0.32,0.4,0.6,0.8,1,1.2]
    h_new = PlotUtils.MnvH2D("Eavail_q3_MErebin",";E_{avail} (GeV);q_{3} (GeV);(10^{-39})",len(Eavail_bin)-1,array('d',Eavail_bin),len(q3_bin_target)-1,array('d',q3_bin_target))
    for x in range(1,len(Eavail_bin)):
        h_new.SetBinContent(x,1,h.GetBinContent(1,x))
        h_new.SetBinError(x,1,h.GetBinError(1,x))
        h_new.SetBinContent(x,2,(h.GetBinContent(2,x)+h.GetBinContent(3,x))/2)
        #print(h.GetBinContent(2,x),h.GetBinContent(3,x),h_new.GetBinContent(x,2))
        index =x*6-6+1
        if cov:
            err = math.sqrt(sum(cov[a][b] for a in [index,index+1] for b in [index,index+1]))/2
        else:
            err = math.sqrt(h.GetBinError(2,x)**2+h.GetBinError(3,x)**2)/2
        h_new.SetBinError(x,2,err)
        h_new.SetBinContent(x,3,h.GetBinContent(4,x))
        h_new.SetBinError(x,3,h.GetBinError(4,x))
        h_new.SetBinContent(x,4,h.GetBinContent(5,x))
        h_new.SetBinError(x,4,h.GetBinError(5,x))
    h_new.Scale(1e39)
    return h_new




def forfun():
    nuefile = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/xsec_me1D_nx1_electron_MAD.root")
    numufile = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/xsec_me1D_nx_muon_MAD.root")
    name ="Eavail_q3"

    hists = []
    for f in [nuefile,numufile]:
        hists.append(f.Get("{}_dataxsec".format(name)))
        hists.append(f.Get("{}_mcxsec".format(name)))
    
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,xsec_ratio,hists,draw_seperate_legend=True)
    PlotTools.Print("{}_xsecemu_ratio".format(name))

def compXsec():
    def Draw(mnvplotter,he,hmu,hmuME,*mc):
        m = max(map(lambda h: h.GetMaximum(),[he,hmu,hmuME])) * 1.1
        leg = ROOT.TLegend(0.5,0.7,0.8,0.9).Clone()
        he.SetLineColor(ROOT.kRed+2)
        he.SetMaximum(m)
        he.Draw("E1 SAME")
        leg.AddEntry(he,"#nu_{e} ME","LE")
        hmu.SetLineColor(ROOT.kGreen-2)
        #hmu.SetMaximum(m)
        hmu.SetMarkerStyle(1)
        hmu.Draw("E1 SAME")
        leg.AddEntry(hmu,"#nu_{#mu} LE","LE")
        hmuME.SetLineColor(ROOT.kBlue-2)
        #hmu.SetMaximum(m)
        hmuME.SetMarkerStyle(1)
        hmuME.Draw("E1 SAME")
        leg.AddEntry(hmuME,"#nu_{#mu} ME","LE")
        colors = ROOT.MnvColors.GetColors()
        for i,h in enumerate(mc):
            h.SetLineColor(colors[i])
            h.Draw("HIST C SAME")
            leg.AddEntry(h,DrawMC[i],"L")
        leg.Draw()
    nuefile = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/xsec_meFHC_col13.3_MAD.root")
    numufile = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_LE.root")
    numuMEfile = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_ME.root")
    numuMCfile = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_MCME.root")
    he = nuefile.Get("Eavail_q3_dataxsec").GetCVHistoWithError().Clone("nue")
    hmu = numufile.Get("Eavail_q3_LErebin").Clone("numu")
    hmuME = numuMEfile.Get("Eavail_q3_MErebin").Clone("numu")
    hmu.Scale(10) #convert to 10^-39
    hists = [he,hmu,hmuME]
    for i in DrawMC:
        hists.append(numuMCfile.Get(i).Clone(i))
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,Draw,hists,draw_seperate_legend=True)
    PlotTools.Print("Eavail_q3_xseccomp")

def GetLEXsec():
    f = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_LE.root","RECREATE")
    hist,hist_rebined = ReadLEResult("supp.txt")
    f.cd()
    hist.Write()
    hist_rebined.Write()
    f.Close()

def GetMEXsec():
    f = ROOT.TFile.Open("/minerva/data/users/hsu/supplementalData_LowRecoilME.root")
    hmu = f.Get("crossSection")
    cov = f.Get("covarianceMatrix")
    hrebined = RebinHistogram(hmu,cov)
    fp = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_ME.root","RECREATE")
    fp.cd()
    hrebined.Write()
    fp.Close
    f.Close()

def GetMEMCXsec():
    f =  ROOT.TFile.Open("/minerva/data/users/hsu/MC_LowRecoilAnalysis.root")
    fp = ROOT.TFile.Open("/minerva/data/users/hsu/numu_xsec_MCME.root","RECREATE")
    pairs = {"MnvTune-v1.2":"MCmnvtuneV12",
             "MnvTune-v3":"MCmnvtuneV3",
            "NuWro-SF": "MCnuwroSF",
            "NuWro-LFG": "MCnuwroLFG",
            "GENIE3-10a":"MCgenie310a",
            "GENIE3-10b":"MCgenie310b",
            "GENIE3-02a":"MCgenie302a",
            }
    for k,v in pairs.items():
        hmu = f.Get(v)
        hrebined = RebinHistogram(hmu,None)
        fp.cd()
        hrebined.SetTitle(k)
        hrebined.Write(k)
    fp.Close()
    f.Close()


if __name__ == "__main__":
    #GetLEXsec()
    #GetMEXsec()
    #GetMEMCXsec()
    compXsec()
