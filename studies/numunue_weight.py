import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities
from functools import partial


ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx2_muon_MAD.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx2_electron_MAD.root"
nuesignalrichfile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_BigNuE2_electron_MAD.root"
nuebkgtunedfile =  "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx1_electron_MAD.root"
nuedatafile = "/minerva/data/users/hsu/nu_e/kin_dist_datame1D_nx2_electron_MAD.root"
numuweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx2_muon_true_weighted_MAD.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datame1D_nx2_muon_weighted_MAD.root"
}

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]
E_BACKGROUNDS = ["NuEElastic","NonFiducial","NonPhaseSpace","CCNuEAntiNu",
                 "ExcessModel","CCDIS","CCOther","NCCOH","NCDIS","NCRES","NCOther"]
MU_BACKGROUNDS = ["NC","CCNuE","CCAntiNuMu","NonFiducial","NonPhaseSpace","Other"]
PLOTPATH = "/minerva/data/users/hsu/numunueRatioPlots/"

SIGNAL_CHANNEL_PAIR = [
    ("CCQE","CCNuEQE","QE"),
    ("CCDelta","CCNuEDelta","Res"),
    ("CCDIS","CCNuEDIS","DIS"),
    ("CCOther","CCNuE","Others"),
    ("CC2p2h","CCNuE2p2h","2p2h"),
]

COLORS=ROOT.MnvColors.GetColors()


#E_SIGNALS = MU_SIGNALS

def relativePOT(nue_true_signal,nue_true_signal_bignue,nue_POT):
    f = lambda hist:hist.GetEntries()
    #f = lambda hist:hist.Integral(0,-1,0,-1)
    return nue_POT * f(nue_true_signal_bignue)/f(nue_true_signal)

def SubtractPoissonHistograms(h,h1):
    h.AddMissingErrorBandsAndFillWithCV(h1)
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        if errors[i]==0:
            continue
        h.SetBinError(i,errors[i])
    return h

def repeat(f,n):
    def rf(x):
        for i in range(n):
            x = f(x)
        return x
    return rf

def GetSignalHist(f,signals,hist_name):
    h = f.Get(hist_name).Clone()
    h.Reset()
    for i in signals:
        t = f.Get("{}_{}".format(hist_name,i))
        if not t:
            print (hist_name,i)
        h.Add(t)
    return h

def MakeScaleFile(hist_names):
    Fout = ROOT.TFile.Open("emu_scale.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(numufile)
    fesignalrich = ROOT.TFile.Open(nuesignalrichfile)
    esignalrichPOT = relativePOT(fe.Get("Eavail_q3_true_signal"),fesignalrich.Get("Eavail_q3_true_signal"),ePOT)
    print(ePOT,muPOT)
    for hist_name in hist_names:
        he = GetSignalHist(fesignalrich,E_SIGNALS,hist_name)
        hmu = GetSignalHist(fmu,MU_SIGNALS,hist_name)
        he.AddMissingErrorBandsAndFillWithCV(hmu)
        hmu.AddMissingErrorBandsAndFillWithCV(he)
        hmu.Scale(esignalrichPOT/muPOT)
        Fout.cd()
        print([h.GetBinContent(3) for h in [he,hmu]])
        he.Divide(he,hmu)
        he.Write(hist_name)
    fe.Close()
    fmu.Close()
    Fout.Close()

def MakeScaleFile2():
    Fout = ROOT.TFile.Open("emu_scale2.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(numufile)
    fesignalrich = ROOT.TFile.Open(nuesignalrichfile)
    esignalrichPOT = relativePOT(fe.Get("Eavail_q3_true_signal"),fesignalrich.Get("Eavail_q3_true_signal"),ePOT)
    print(ePOT,muPOT,esignalrichPOT)
    he = GetSignalHist(fesignalrich,E_SIGNALS,"tEnu")
    migration = fmu.Get("Enu_migration").Clone()
    he.Scale(muPOT/esignalrichPOT)
    tmp = FindScale2(he,migration)
    scale_hist = fmu.Get("Enu").Clone()
    for i in range(scale_hist.GetNbinsX()+2):
        scale_hist.SetBinContent(i,tmp[i][0])
        scale_hist.SetBinError(i,0)
    Fout.cd()
    scale_hist.Write("Enu")
    fe.Close()
    fmu.Close()
    Fout.Close()


def FindScale2(target,migration):
    m = ROOT.TMatrixD(migration.GetNbinsY()+2,migration.GetNbinsX()+2)
    t = ROOT.TMatrixD(target.GetNbinsX()+2,1)

    for x in range(migration.GetNbinsX()+2):
        for y in range(migration.GetNbinsY()+2):
            m[y][x]=migration.GetBinContent(x,y)

    for y in range(target.GetNbinsX()+2):
        t[y][0] = target.GetBinContent(y)

    SVD = ROOT.TDecompSVD(m)
    tmp = SVD.Invert()
    r = ROOT.TMatrixD(migration.GetNbinsX()+2,1)
    r.Mult(tmp,t)
    return repeat(smooth,0)(r)

def smooth(iput):
    tmp = iput[11][0]
    for i in range(12,iput.GetNrows()-1):
        new = (tmp+iput[i][0]+iput[i+1][0])/3
        tmp = iput[i][0]
        iput[i][0]=new
    iput[iput.GetNrows()-1][0]=(tmp+iput[iput.GetNrows()-1][0])/2
    print (list(iput[i][0] for i in range(iput.GetNrows())))
    return iput

def MakeCompPlot():
    def getFileAndPOT(filename):
        POT = Utilities.getPOTFromFile(filename)
        f = ROOT.TFile.Open(filename)
        return f,POT
    def Draw(mnvplotter,*args,canvas=PlotTools.CANVAS):
        mc_hist1 = args[0]
        mc_hist2 = args[1]
        data_hist1 = args[2] if len(args)>2 else None
        data_hist2 = args[3] if len(args)>3 else None
        leg = ROOT.TLegend(0.5,0.7,0.8,0.9).Clone()
        # if data_hist:
        #     data_hist.Draw("E1 X0")

        # for h in args:
        #     if h is None:
        #         continue
        #     h.GetXaxis().SetRange(1,h.GetNbinsX()+1)

        m = max(map(lambda h: h.GetMaximum(),[h for h in args if h is not None])) * 1.1
        for h in args:
            if h is None:
                continue
            h.SetLineWidth(3)
            h.SetMaximum(m)


        if data_hist1:
            #data_hist1.SetTitle("nue Data")
            data_hist1.SetLineColor(ROOT.kRed+2)
            # data_hist1.SetMaximum(m)
            # data_hist1.SetLineWidth(3)
            #data_hist1.GetXaxis().SetRangeUser(0,data_hist1.GetNbinsX()+1)
            #data_hist1.SetMarkerColor(ROOT.kRed)
            data_hist1.Draw("SAME")
            leg.AddEntry(data_hist1,"nue data")
        if data_hist2:
            #data_hist2.SetTitle("numu Data")
            data_hist2.SetLineColor(ROOT.kGreen-2)
            #data_hist2.SetMarkerColor(ROOT.kGreen)
            # data_hist2.SetMaximum(m)
            # data_hist2.SetLineWidth(3)
            #data_hist2.GetXaxis().SetRangeUser(0,data_hist2.GetNbinsX()+1)
            data_hist2.Draw("SAME")
            leg.AddEntry(data_hist2,"weighted numu data")
        if mc_hist1:
            #mc_hist1.SetTitle("nue MC")
            mc_hist1.SetLineColor(ROOT.kRed)
            # mc_hist1.SetMaximum(m)
            # mc_hist1.SetLineWidth(3)
            #mc_hist1.GetXaxis().SetRangeUser(0,mc_hist1.GetNbinsX()+1)
            mc_hist1.Draw("HIST E1 SAME")
            leg.AddEntry(mc_hist1,"nue MC")
        if mc_hist2:
            #mc_hist2.SetTitle("weighted numu MC")
            mc_hist2.SetLineColor(ROOT.kGreen)
            # mc_hist2.SetMaximum(m)
            # mc_hist2.SetLineWidth(3)
            #mc_hist2.GetXaxis().SetRangeUser(0,mc_hist2.GetNbinsX()+1)
            mc_hist2.Draw("HIST E1 SAME")
            leg.AddEntry(mc_hist2,"weighted numu MC")
        leg.Draw()

    def DrawRatio(mnvplotter,h1,h2,include_systematics = False):
        cast = (lambda x:x.GetCVHistoWithError()) if include_systematics else (lambda x:x.GetCVHistoWithStatError()) 
        mnvplotter.DrawDataMCRatio(cast(h1),cast(h2), 1.0 ,True,0,2)

    def DrawDoubleRatio(mnvplotter,h1,h2,h3,h4,include_systematics = False):
        h_r1 = h1.Clone("{}_ratio1".format(h1.GetName))
        h_r2 = h3.Clone("{}_ratio2".format(h3.GetName))
        h_r1.Divide(h_r1,h2)
        h_r2.Divide(h_r2,h4)
        DrawRatio(mnvplotter,h_r1,h_r2,include_systematics)

    def DrawFraction(mnvplotter,*cates,titles,colors,canvas=PlotTools.CANVAS):
        leg = ROOT.TLegend(0.5,0.7,0.9,0.9).Clone()
        for i,v in enumerate(cates):
            v.SetLineColor(colors[i])
            v.SetLineWidth(3)
            v.SetMaximum(1)
            v.Draw("HIST E1 SAME")
            leg.AddEntry(v,titles[i])
        leg.Draw()

    def GetFractions(histfile,cates,hist_name):
        htotal = GetSignalHist(histfile,cates,hist_name)
        hcomps =[]
        for i in cates:
            hcomps.append(GetSignalHist(histfile,[i],hist_name))
        return PlotTools.Ratio(hcomps,htotal,"B")

    def DrawEfficiency(mnvplotter,henu,hede,hmunu,hmude,canvas=PlotTools.CANVAS):
        henu.Divide(henu,hede,1,1,"B")
        hmunu.Divide(hmunu,hmude,1,1,"B")
        Draw(mnvplotter,henu,hmunu)


    nueMCfile,nueMCPOT = getFileAndPOT(nuefile)
    nueMCsignalfile,nueMCsignalPOT = getFileAndPOT(nuesignalrichfile)
    nueMCsignalPOT = relativePOT(nueMCfile.Get("Eavail_q3_true_signal"),nueMCsignalfile.Get("Eavail_q3_true_signal"),nueMCPOT)
    nueBkgfile,nueBkgPOT = getFileAndPOT(nuebkgtunedfile)
    numuMCfile,numuMCPOT = getFileAndPOT(numuweightedfile["mc"])
    numuDatafile,numuDataPOT = getFileAndPOT(numuweightedfile["data"])
    nueDatafile,nueDataPOT = getFileAndPOT(nuedatafile)
    numuMCrawfile,numuMCrawDataPOT = getFileAndPOT(numufile)
    print (nueMCPOT,numuDataPOT,numuMCPOT,nueDataPOT)

    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","Q3","Q0","tQ0","tQ3","tLepton_Pt","tQ0_tQ3","tQ0_tLepton_Pt","Eavail","tEavail","tEavail_tQ3","tEavail_tLepton_Pt","Lepton_Pt","tEnu","Eavail_q3_true_signal", "Eavail_Lepton_Pt_true_signal", "tEnu_true_signal"]:
        try:
            he = GetSignalHist(nueMCsignalfile,E_SIGNALS,hist_name)
            hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
        except :
            print ("hist: {} not found".format(hist_name))
            continue
        hmu.Scale(numuDataPOT/numuMCPOT)
        he.Scale(numuDataPOT/nueMCsignalPOT)
        hmudata = numuDatafile.Get(hist_name)
        if hmudata:
            hmubkg = GetSignalHist(numuMCfile,MU_BACKGROUNDS,hist_name)
            hmubkg.Scale(numuDataPOT/numuMCPOT)
            SubtractPoissonHistograms(hmudata,hmubkg)
        hedata = nueDatafile.Get(hist_name)
        if hedata:
            hedata.Scale(numuDataPOT/nueDataPOT)
            hebkg = GetSignalHist(nueBkgfile,E_BACKGROUNDS,hist_name)
            hebkg.Scale(numuDataPOT/nueBkgPOT)
            SubtractPoissonHistograms(hedata,hebkg)

        if "true_signal" in hist_name:
            sc = he.Integral()/hmu.Integral()
            hmu.Scale(sc)
        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        if hedata:
            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,he],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_heratio.png".format(PLOTPATH,hist_name))
        if hmudata:
            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hmudata,hmu],draw_seperate_legend = hmu.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_hmuratio.png".format(PLOTPATH,hist_name))

        if hedata and hmudata:
            PlotTools.MakeGridPlot(Slicer,DrawDoubleRatio,[hedata,he,hmudata,hmu],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_doubleratio.png".format(PLOTPATH,hist_name))

            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_dataratio.png".format(PLOTPATH,hist_name))
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu,hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
        else:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}.png".format(PLOTPATH,hist_name))
        PlotTools.MakeGridPlot(Slicer,DrawRatio,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_MCratio.png".format(PLOTPATH,hist_name))

        
        for numu_sig,nue_sig,_ in SIGNAL_CHANNEL_PAIR:
            he = GetSignalHist(nueMCsignalfile,[nue_sig],hist_name)
            hmu = GetSignalHist(numuMCfile,[numu_sig],hist_name)
            hmu.Scale(numuDataPOT/numuMCPOT)
            he.Scale(numuDataPOT/nueMCsignalPOT)
            if "true_signal" in hist_name:
                hmu.Scale(sc)
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_{}MC.png".format(PLOTPATH,hist_name,numu_sig))
            drawer = partial(DrawRatio, include_systematics=False)
            PlotTools.MakeGridPlot(Slicer,drawer,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_{}MCratio.png".format(PLOTPATH,hist_name,numu_sig))

        titles = []
        mu_cates = []
        e_cates = []
        colors = []
        for i,v  in enumerate(SIGNAL_CHANNEL_PAIR):
            (numu_sig,nue_sig,title) = v
            titles.append(title)
            e_cates.append(nue_sig)
            mu_cates.append(numu_sig)
            colors.append(COLORS[i])
        Fraction_e = GetFractions(nueMCsignalfile,e_cates,hist_name)
        Fraction_mu = GetFractions(numuMCfile,mu_cates,hist_name)
        drawer = partial(DrawFraction,titles=titles,colors=colors)
        PlotTools.MakeGridPlot(Slicer,drawer,Fraction_e,draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_MCfractionE.png".format(PLOTPATH,hist_name))
        PlotTools.MakeGridPlot(Slicer,drawer,Fraction_mu,draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_MCfractionMU.png".format(PLOTPATH,hist_name))


    #draw efficiency:
    for i in [("Eavail_q3","tEavail_tQ3","Eavail_q3_true_signal"),
              ("Eavail_Lepton_Pt","tEavail_tLepton_Pt","Eavail_Lepton_Pt_true_signal")]:
        Slicer = PlotTools.Make2DSlice
        hes =list(map(lambda h: GetSignalHist(nueMCsignalfile,E_SIGNALS,h), i))
        hmus = list(map(lambda h: GetSignalHist(numuMCrawfile,MU_SIGNALS,h), i))
        PlotTools.MakeGridPlot(Slicer,DrawEfficiency,[hes[1],hes[2],hmus[1],hmus[2]],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_efficiency.png".format(PLOTPATH,i[0]))
        PlotTools.MakeGridPlot(Slicer,DrawEfficiency,[hes[0],hes[2],hmus[0],hmus[2]],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_recoefficiency.png".format(PLOTPATH,i[0]))

if __name__ == "__main__":
    MakeCompPlot()
    #MakeScaleFile(["Enu","tEnu","tEnu_true_signal"])
    #MakeScaleFile2()

