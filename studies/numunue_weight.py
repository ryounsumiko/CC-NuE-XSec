import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities
from functools import partial
from array import array

# Required input files:
# numufile: unweighted mc numu sample, DO NOT apply background constraint scale factors.
# nuefile: unweighted mc nue sample, DO NOT apply background constraint scale factors.
# nuesignalrichfile: unweighted mc bignue sample, DO NOT apply background constraint scale factors.
# nuebkgtunedfile: unweighted mc nue sample, DO apply background constraint scale factors.
# nuedatafile: data nuesample.
# numuweightedfile: Enu weighted data/mc numu sample.


ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_muon_MAD.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_electron_MAD.root"
nuesignalrichfile = "/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_BigNuE_electron_MAD.root"
nuebkgtunedfile =  "/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_electron_bkgtuned_MAD.root"
nuedatafile = "/minerva/data/users/hsu/nu_e/kin_dist_datameFHC1_electron_bkgtuned_MAD.root"
numuweightedfile =  {
    #"mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC_muon_true_weight_MAD.root",
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_muon_weight_cor_MAD.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datameFHC1_muon_weight_cor_MAD.root"
}

numuTrueweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC1_muon_true_weight_MAD.root",
}

# Defines the signal and background categories used in numu and nue sample. 

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]
E_BACKGROUNDS = ["NuEElastic","NonFiducial","NonPhaseSpace","CCNuEAntiNu",
                 "ExcessModel","CCDIS","CCOther","NCCOH","NCDIS","NCRES","NCOther"]
MU_BACKGROUNDS = ["NC","CCNuE","CCAntiNuMu","NonFiducial","NonPhaseSpace","Other"]

#PLOT OUTPUT PATH 
PLOTPATH = "/minerva/data/users/{}/numunueRatioPlots/".format(os.environ["USER"])

DRAW_SYSTEMATICS = True
fluxbin = array('d',[0.5*i for i in range(4,21)]+[i for i in range(11,21)])

# Matching signal channels. format:
# ("category name used in numu sample", "category name used in nue sample","category name in plot name")

SIGNAL_CHANNEL_PAIR = [
    ("CCQE","CCNuEQE","QE"),
    ("CCDelta","CCNuEDelta","Res"),
    ("CCDIS","CCNuEDIS","DIS"),
    ("CCOther","CCNuE","Others"),
    ("CC2p2h","CCNuE2p2h","2p2h"),
]

COLORS=ROOT.MnvColors.GetColors()

ERROR_BANDS_TO_POP = [
    "MK_model",
    "fsi_weight",
    "SuSA_Valencia_Weight",
    "LowQ2Pi_JOINT",
    "LowQ2Pi_NUPI0",
]

SYSTEMATIC_ERROR_GROUPS = {
    "Interaction model": [
        "GENIE_"+ i for i in [
             "AGKYxF1pi",
            "AhtBY",
            "BhtBY",
            "CCQEPauliSupViaKF",
            "CV1uBY",
            "CV2uBY",
            "EtaNCEL",
            "FrAbs_N",
            "FrAbs_pi",
            "FrCEx_N",
            "FrCEx_pi",
            "FrElas_N",
            "FrElas_pi",
            "FrInel_N",
            #"FrInel_pi",
            "FrPiProd_N",
            "FrPiProd_pi",
            "MFP_N",
            "MFP_pi",
            "MaCCQE",
            "MaNCEL",
            "MaRES",
            "MvRES",
            "NormCCRES",
            "NormDISCC",
            "NormNCRES",
            "RDecBR1gamma",
            "Rvn1pi",
            "Rvp1pi",
            "Rvn2pi",
            "Rvp2pi",
            "Theta_Delta2Npi",
            "VecFFCCQEshape"
        ]
    ],
    "Mnv Tunes" : [
        "LowQ2Pi",
        "Low_Recoil_2p2h_Tune",
        "RPA_HighQ2",
        "RPA_LowQ2"
    ],
    "Detector model" : [
        "eltheta",
        "beam_angle",
        "elE_ECAL","elE_HCAL",
        "Leakage_Uncertainty",
        "Target_Mass_CH",
        "response_p",
        "response_meson",
        "response_em",
        "response_other",
        "response_xtalk",
        "GEANT_Neutron",
        "GEANT_Pion",
        "GEANT_Proton",
        "Muon_Energy_MINERvA",
        "Muon_Energy_MINOS",
        "MuonAngleXResolution",
        "MuonAngleYResolution",
        "Muon_Energy_Resolution",
        "MINOS_Reconstruction_Efficiency",
    ],
    "Alternative Tunning methods":[
        "bkg_tune"
    ]
}


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
        del t
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
    # fe = ROOT.TFile.Open(nuefile)
    # ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(numufile)
    # fesignalrich = ROOT.TFile.Open(nuesignalrichfile)
    # esignalrichPOT = relativePOT(fe.Get("Eavail_q3_true_signal"),fesignalrich.Get("Eavail_q3_true_signal"),ePOT)
    # print(ePOT,muPOT,esignalrichPOT)
    # he = GetSignalHist(fesignalrich,E_SIGNALS,"tEnu")
    migration = fmu.Get("Enu_migration").Clone()
    # he.Scale(muPOT/esignalrichPOT)
    fflux = ROOT.TFile.Open("flux_scale.root")
    ratio = fflux.Get("FHC")
    tmp = FindScale2(ratio,migration)
    scale_hist = fmu.Get("Enu").Clone()
    for i in range(scale_hist.GetNbinsX()+2):
        scale_hist.SetBinContent(i,tmp[i][0])
        scale_hist.SetBinError(i,0)
    Fout.cd()
    scale_hist.Write("Enu")
    fflux.Close()
    fmu.Close()
    Fout.Close()


def FindScale2(target,migration):
    m = ROOT.TMatrixD(migration.GetNbinsY()+2,migration.GetNbinsX()+2)
    t = ROOT.TMatrixD(target.GetNbinsX()+2,1)

    for y in range(migration.GetNbinsY()+2):
        total = 0
        for x in range(migration.GetNbinsX()+2):
            m[y][x]=migration.GetBinContent(x,y)
            total += m[y][x]
        if not total:
            continue
        for x in range(migration.GetNbinsX()+2):
            m[y][x]=m[y][x]/total

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

def GetEfficiencyCorrection():
    Fout = ROOT.TFile.Open("eff_corr_scale.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numuTrueweightedfile["mc"])
    muPOT = Utilities.getPOTFromFile(numuTrueweightedfile["mc"])
    fesignalrich = ROOT.TFile.Open(nuesignalrichfile)
    esignalrichPOT = relativePOT(fe.Get("Eavail_q3_true_signal"),fesignalrich.Get("Eavail_q3_true_signal"),ePOT)
    for i in [("Eavail_q3","tEavail_tQ3","Eavail_q3_true_signal"),
              ("Eavail_Lepton_Pt","tEavail_tLepton_Pt","Eavail_Lepton_Pt_true_signal"),
              ("Enu","tEnu","tEnu_true_signal"),]:

        #Slicer = PlotTools.Make2DSlice
        hes =list(map(lambda h: GetSignalHist(fesignalrich,E_SIGNALS,h), i))
        hmus = list(map(lambda h: GetSignalHist(fmu,MU_SIGNALS,h), i))
        names = {
            "eff_e":hes,
            "eff_mu":hmus,
        }

        for k,v in names.items():
            for j in [0,1]:
                div = v[2]
                v[j].Divide(v[j],div)
                Fout.cd()
                v[j].Write("{}_{}".format(i[j],k))

def GetEfficiency(hist_name,is_muon,is_reco_weight):
    Feff = ROOT.TFile.Open("eff_corr_scale.root")
    name = "{}_eff_{}".format(hist_name,"mu" if is_muon else "e")
    htemp = Feff.Get(name)
    Feff.Close()
    return htemp

def myRebin(hist,binning):
    h2 = ROOT.TH1D("{}_rebin".format(hist.GetName()),hist.GetTitle(),len(binning)-1,array('d',binning))
    for i in range(0,hist.GetSize()):
        h2.Fill(hist.GetXaxis().GetBinCenter(i),hist.GetBinContent(i))
    return h2

def MakeCompPlot():
    def getFileAndPOT(filename):
        POT = Utilities.getPOTFromFile(filename)
        f = ROOT.TFile.Open(filename)
        return f,POT

    def DrawTWRW(mnvplotter,*args,canvas=PlotTools.CANVAS):
        rw = args[0]
        tw = args[1]
        e = args[2] if len(args)>2 else None
        leg = ROOT.TLegend(0.5,0.7,0.8,0.9).Clone()

        m = max(map(lambda h: h.GetMaximum(),[h for h in args if h is not None])) * 1.1
        for h in args:
            if h is None:
                continue
            h.SetLineWidth(3)
            h.SetMaximum(m)

        rw.SetLineColor(ROOT.kGreen)
        tw.SetLineColor(ROOT.kBlue)
        if e:
            e.SetLineColor(ROOT.kRed)
            e.Draw("HIST SAME")
            leg.AddEntry(e,"nue MC")
        tw.Draw("HIST SAME")
        leg.AddEntry(tw,"numu MC tw")
        rw.Draw("HIST SAME")
        leg.AddEntry(rw,"numu MC rw")
        leg.Draw()

    def DrawStages(mnvplotter,*args,canvas=PlotTools.CANVAS):
        ed = args[0]
        en = args[1]
        r = args[2]
        leg = ROOT.TLegend(0.5,0.7,0.8,0.9).Clone()

        ed.Draw("HIST SAME")
        leg.AddEntry(ed,"GENIE event rate (ED)")
        ed.SetLineColor(ROOT.kGreen)
        en.SetLineColor(ROOT.kBlue)
        r.SetLineColor(ROOT.kRed)
        r.Draw("HIST SAME")
        leg.AddEntry(r,"Reco distribution of sample (R)")
        en.Draw("HIST SAME")
        leg.AddEntry(en,"True distribution of sample (EN)")
        leg.Draw()
        pass

    def DrawRecoVsTrue(mnvplotter,*args,canvas=PlotTools.CANVAS):
        en = args[0]
        r = args[1]
        leg = ROOT.TLegend(0.5,0.7,0.8,0.9).Clone()

        m = max(map(lambda h: h.GetMaximum(),[h for h in args if h is not None])) * 1.1
        for h in args:
            if h is None:
                continue
            h.SetLineWidth(3)
            h.SetMaximum(m)
            h.GetXaxis().SetTitle("E_{lep}")

        en.SetLineColor(ROOT.kBlue)
        r.SetLineColor(ROOT.kRed)
        r.Draw("HIST SAME")
        leg.AddEntry(r,"Reco Emu")
        en.Draw("HIST SAME")
        leg.AddEntry(en,"True Emu")
        leg.Draw()
        pass

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



    def DrawRatio(mnvplotter,h1,h2,include_systematics = DRAW_SYSTEMATICS,ratio_title = None):
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        cast = (lambda x:x) if include_systematics else (lambda x:x.GetCVHistoWithStatError())
        mnvplotter.DrawDataMCRatio(cast(h1),cast(h2), 1.0 ,True,True,0,2 , ratio_title or "#nu_{mu}/#nu_{#e}")

    def DrawErrorRatio(mnvplotter,h1,h2,include_systematics = DRAW_SYSTEMATICS,ratio_title = None):
        cast = (lambda x:x) if include_systematics else (lambda x:x.GetCVHistoWithStatError())
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        h = cast(h1)
        h.Divide(h,cast(h2))
        mnvplotter.DrawErrorSummary(h)

    def DrawRatiosOnSamePlot(mnvplotter,h1,h2,h3,h4,include_systematics = DRAW_SYSTEMATICS):
        cast = (lambda x:x) if include_systematics else (lambda x:x.GetCVHistoWithStatError())
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        h =  cast(h1)
        h.Divide(h,cast(h2))
        h.GetYaxis().SetTitle("#nu_{#mu}/#nu_{e}")
        h3.AddMissingErrorBandsAndFillWithCV(h4)
        h4.AddMissingErrorBandsAndFillWithCV(h3)
        hp =  cast(h3)
        hp.Divide(hp,cast(h4))
        hp.GetYaxis().SetTitle("#nu_{#mu}/#nu_{e}")
        mnvplotter.axis_maximum = 2.0
        mnvplotter.DrawDataMCWithErrorBand(h,hp,1.0,"TR",False,ROOT.nullptr,ROOT.nullptr,False,True)
        mnvplotter.axis_maximum = -1111



    def DrawDoubleRatio(mnvplotter,h1,h2,h3,h4,include_systematics = DRAW_SYSTEMATICS):
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        h3.AddMissingErrorBandsAndFillWithCV(h4)
        h4.AddMissingErrorBandsAndFillWithCV(h3)
        h_r1 = h1.Clone("{}_ratio1".format(h1.GetName))
        h_r2 = h3.Clone("{}_ratio2".format(h3.GetName))
        h_r1.Divide(h_r1,h2)
        h_r2.Divide(h_r2,h4)
        DrawRatio(mnvplotter,h_r1,h_r2,include_systematics = DRAW_SYSTEMATICS,ratio_title = "#frac{Data^{#mu}}{Data^{e}}/#frac{MC^{#mu}}{MC^{e}}")
        #h_r1.Divide(h_r1,h_r2)
        #h_r1.GetYaxis().SetTitle("#frac{Data^{#mu}}{Data^{e}}/#frac{MC^{#mu}}{MC^{e}}")
        #h_r1.Draw("")
        #ROOT.SetOwnership(h_r1,False)
        #del h_r2

    def DrawDoubleRatioError(mnvplotter,h1,h2,h3,h4,include_systematics = DRAW_SYSTEMATICS):
        h1.AddMissingErrorBandsAndFillWithCV(h2)
        h2.AddMissingErrorBandsAndFillWithCV(h1)
        h3.AddMissingErrorBandsAndFillWithCV(h4)
        h4.AddMissingErrorBandsAndFillWithCV(h3)
        h_r1 = h1.Clone("{}_ratio1".format(h1.GetName))
        h_r2 = h3.Clone("{}_ratio2".format(h3.GetName))
        h_r1.Divide(h_r1,h2)
        h_r2.Divide(h_r2,h4)
        h_r1.Divide(h_r1,h_r2)
        mnvplotter.DrawErrorSummary(h_r1)
        del h_r1
        del h_r2

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
    numuTrueWeightedfile,numuTrueWeightedPOT = getFileAndPOT(numuTrueweightedfile["mc"])
    print (nueMCPOT,numuDataPOT,numuMCPOT,nueDataPOT)


    #DrawExample

    #for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","Q3","Q0","tQ0","tQ3","tLepton_Pt","tQ0_tQ3","tQ0_tLepton_Pt","Eavail","tEavail","tEavail_tQ3","tEavail_tLepton_Pt","Lepton_Pt","tEnu","Eavail_q3_true_signal", "Eavail_Lepton_Pt_true_signal", "tEnu_true_signal"]:
    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","tEavail_tQ3","tEavail_tLepton_Pt","Lepton_Pt","tEnu","Eavail_q3_true_signal", "Eavail_Lepton_Pt_true_signal", "tEnu_true_signal"]:
        try:
            he = GetSignalHist(nueMCsignalfile,E_SIGNALS,hist_name)
            hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
        except :
            print ("hist: {} not found".format(hist_name))
            continue
        hmu.Scale(numuDataPOT/numuMCPOT,"width")
        he.Scale(numuDataPOT/nueMCsignalPOT,"width")
        hmudata = numuDatafile.Get(hist_name)
        if hmudata:
            hmudata.Scale(numuDataPOT/numuDataPOT,"width")
            hmubkg = GetSignalHist(numuMCfile,MU_BACKGROUNDS,hist_name)
            hmubkg.Scale(numuDataPOT/numuMCPOT,"width")
            SubtractPoissonHistograms(hmudata,hmubkg)
        hedata = nueDatafile.Get(hist_name)
        if hedata:
            hedata.Scale(numuDataPOT/nueDataPOT,"width")
            hebkg = GetSignalHist(nueBkgfile,E_BACKGROUNDS,hist_name)
            hebkg.Scale(numuDataPOT/nueBkgPOT,"width")
            SubtractPoissonHistograms(hedata,hebkg)

        if "true_signal" in hist_name:
            sc = he.Integral()/hmu.Integral()
        #     hmu.Scale(sc)
        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        # if hedata:
        #     PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,he],draw_seperate_legend = he.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_heratio.png".format(PLOTPATH,hist_name))
        # if hmudata:
        #     PlotTools.MakeGridPlot(Slicer,DrawRatio,[hmudata,hmu],draw_seperate_legend = hmu.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_hmuratio.png".format(PLOTPATH,hist_name))

        # if hedata and hmudata:
        #     PlotTools.MakeGridPlot(Slicer,DrawDoubleRatio,[hedata,he,hmudata,hmu],draw_seperate_legend = he.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_doubleratio.png".format(PLOTPATH,hist_name))

        #     PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_dataratio.png".format(PLOTPATH,hist_name))
        #     PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu,hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
        # else:
        #     PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        # PlotTools.CANVAS.Print("{}{}.png".format(PLOTPATH,hist_name))
        # PlotTools.MakeGridPlot(Slicer,DrawRatio,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        # PlotTools.CANVAS.Print("{}{}_MCratio.png".format(PLOTPATH,hist_name))

        PlotTools.updatePlotterErrorGroup(SYSTEMATIC_ERROR_GROUPS)

        if hist_name in ["Eavail_q3","Eavail_Lepton_Pt"]:#,"tEavail_tQ3","tEavail_tLepton_Pt","Enu","tEnu"]:
            thist_name = hist_name
           # print(thist_name)
            for reco_weighted in [True]:#,False]:
                eff_mu = GetEfficiency(thist_name,True,reco_weighted)
                eff_e = GetEfficiency(thist_name,False,reco_weighted)
                hetemp = he.Clone()
                eff_e.AddMissingErrorBandsAndFillWithCV(he)
                hetemp.Divide(hetemp,eff_e)
                if hedata:
                    hedatatemp = hedata.Clone()
                    hedatatemp.Divide(hedatatemp,eff_e)
                hmutemp = hmu.Clone()
                eff_mu.AddMissingErrorBandsAndFillWithCV(hmu)
                hmutemp.Divide(hmutemp,eff_mu)
                if hmudata:
                    hmudatatemp = hmudata.Clone()
                    hmudatatemp.Divide(hmudatatemp,eff_mu)

                for h in [hetemp,hmutemp,hedatatemp,hmudatatemp]:
                    if not h:
                        continue
                    for e in ERROR_BANDS_TO_POP:
                        h.PopVertErrorBand(e)
                if hedata and hmudata:
                    PlotTools.MakeGridPlot(Slicer,Draw,[hetemp,hmutemp,hedatatemp,hmudatatemp],draw_seperate_legend = he.GetDimension()==2)
                else:
                    PlotTools.MakeGridPlot(Slicer,Draw,[hetemp,hmutemp],draw_seperate_legend = he.GetDimension()==2)
                PlotTools.CANVAS.Print("{}{}_{}effcor.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                # PlotTools.MakeGridPlot(Slicer,DrawRatio,[hmutemp,hetemp],draw_seperate_legend = he.GetDimension()==2)
                # PlotTools.CANVAS.Print("{}{}_MC{}effcorratio.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                PlotTools.MakeGridPlot(Slicer,DrawErrorRatio,[hmutemp,hetemp],draw_seperate_legend = he.GetDimension()==2)
                PlotTools.CANVAS.Print("{}{}_MC{}effcorratio_error.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                if hedata and hmudata:
                    # PlotTools.MakeGridPlot(Slicer,DrawRatio,[hmudatatemp,hedatatemp],draw_seperate_legend = he.GetDimension()==2)
                    # PlotTools.CANVAS.Print("{}{}_Data{}effcorratio.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                    PlotTools.MakeGridPlot(Slicer,DrawErrorRatio,[hmudatatemp,hedatatemp],draw_seperate_legend = he.GetDimension()==2)
                    PlotTools.CANVAS.Print("{}{}_Data{}effcorratio_error.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                    PlotTools.MakeGridPlot(Slicer,DrawRatiosOnSamePlot,[hmudatatemp,hedatatemp,hmutemp,hetemp],draw_seperate_legend = he.GetDimension()==2)
                    PlotTools.CANVAS.Print("{}{}_{}effcorratio.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                    PlotTools.MakeGridPlot(Slicer,DrawDoubleRatio,[hmudatatemp,hedatatemp,hmutemp,hetemp],draw_seperate_legend = he.GetDimension()==2)
                    PlotTools.CANVAS.Print("{}{}_{}effcordoubleratio.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))
                    PlotTools.MakeGridPlot(Slicer,DrawDoubleRatioError,[hmudatatemp,hedatatemp,hmutemp,hetemp],draw_seperate_legend = he.GetDimension()==2)
                    PlotTools.CANVAS.Print("{}{}_{}effcordoubleratio_error.png".format(PLOTPATH,hist_name,"reco" if reco_weighted else ""))



        
        # for numu_sig,nue_sig,_ in SIGNAL_CHANNEL_PAIR:
        #     he = GetSignalHist(nueMCsignalfile,[nue_sig],hist_name)
        #     hmu = GetSignalHist(numuMCfile,[numu_sig],hist_name)
        #     hmu.Scale(numuDataPOT/numuMCPOT,"width")
        #     he.Scale(numuDataPOT/nueMCsignalPOT,"width")
        #     if "true_signal" in hist_name:
        #         hmu.Scale(sc)
        #     PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_{}MC.png".format(PLOTPATH,hist_name,numu_sig))
        #     drawer = partial(DrawRatio, include_systematics=False)
        #     PlotTools.MakeGridPlot(Slicer,drawer,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        #     PlotTools.CANVAS.Print("{}{}_{}MCratio.png".format(PLOTPATH,hist_name,numu_sig))

        # titles = []
        # mu_cates = []
        # e_cates = []
        # colors = []
        # for i,v  in enumerate(SIGNAL_CHANNEL_PAIR):
        #     (numu_sig,nue_sig,title) = v
        #     titles.append(title)
        #     e_cates.append(nue_sig)
        #     mu_cates.append(numu_sig)
        #     colors.append(COLORS[i])
        # Fraction_e = GetFractions(nueMCsignalfile,e_cates,hist_name)
        # Fraction_mu = GetFractions(numuMCfile,mu_cates,hist_name)
        # drawer = partial(DrawFraction,titles=titles,colors=colors)
        # PlotTools.MakeGridPlot(Slicer,drawer,Fraction_e,draw_seperate_legend = he.GetDimension()==2)
        # PlotTools.CANVAS.Print("{}{}_MCfractionE.png".format(PLOTPATH,hist_name))
        # PlotTools.MakeGridPlot(Slicer,drawer,Fraction_mu,draw_seperate_legend = he.GetDimension()==2)
        # PlotTools.CANVAS.Print("{}{}_MCfractionMU.png".format(PLOTPATH,hist_name))


    #draw efficiency:
    for i in [("Eavail_q3","tEavail_tQ3","Eavail_q3_true_signal"),
              ("Eavail_Lepton_Pt","tEavail_tLepton_Pt","Eavail_Lepton_Pt_true_signal"),
              ("Enu","tEnu","tEnu_true_signal"),]:
        twoD = "Eavail" in i[0]
        Slicer = PlotTools.Make2DSlice if twoD else lambda x: [x]
        hes =list(map(lambda h: GetSignalHist(nueMCsignalfile,E_SIGNALS,h), i))
        hmus = list(map(lambda h: GetSignalHist(numuTrueWeightedfile,MU_SIGNALS,h), i))
        #hmus[2].Scale(0.01)
        PlotTools.MakeGridPlot(Slicer,DrawEfficiency,[hes[1],hes[2],hmus[1],hmus[2]],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_efficiency.png".format(PLOTPATH,i[0]))
        PlotTools.MakeGridPlot(Slicer,DrawEfficiency,[hes[0],hes[2],hmus[0],hmus[2]],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_recoefficiency.png".format(PLOTPATH,i[0]))

        #draw example:
        hmued = GetSignalHist(numuTrueWeightedfile,MU_SIGNALS,i[2])
        hmued.Scale(numuDataPOT/numuTrueWeightedPOT,"width")
        hmuen = GetSignalHist(numuTrueWeightedfile,MU_SIGNALS,i[1])
        hmuen.Scale(numuDataPOT/numuTrueWeightedPOT,"width")
        hmur = GetSignalHist(numuTrueWeightedfile,MU_SIGNALS,i[0])
        hmur.Scale(numuDataPOT/numuTrueWeightedPOT,"width")
        hmuen_rw = GetSignalHist(numuMCfile,MU_SIGNALS,i[1])
        hmuen_rw.Scale(numuDataPOT/numuMCPOT,"width")
        hmur_rw = GetSignalHist(numuMCfile,MU_SIGNALS,i[0])
        hmur_rw.Scale(numuDataPOT/numuMCPOT,"width")
        her = GetSignalHist(nueMCsignalfile,E_SIGNALS,i[0])
        her.Scale(numuDataPOT/nueMCsignalPOT,"width")
        heen = GetSignalHist(nueMCsignalfile,E_SIGNALS,i[1])
        heen.Scale(numuDataPOT/nueMCsignalPOT,"width")
        heed = GetSignalHist(nueMCsignalfile,E_SIGNALS,i[2])
        heed.Scale(numuDataPOT/nueMCsignalPOT,"width")
        
        PlotTools.MakeGridPlot(Slicer,DrawTWRW,[hmuen_rw,hmuen,heen],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_twrw.png".format(PLOTPATH,i[1]))
        PlotTools.MakeGridPlot(Slicer,DrawTWRW,[hmur_rw,hmur,her],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_twrw.png".format(PLOTPATH,i[0]))
        PlotTools.MakeGridPlot(Slicer,DrawStages,[hmued,hmuen,hmur],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_stages_mu.png".format(PLOTPATH,i[0]))
        PlotTools.MakeGridPlot(Slicer,DrawStages,[heed,heen,her],draw_seperate_legend = twoD)
        PlotTools.CANVAS.Print("{}{}_stages_e.png".format(PLOTPATH,i[0]))

    for i in ["tEnu_tEavail","Enu_tEavail", "tEnu_Eavail","Enu_Eavail",
              "tEnu_tLepton_Pt","Enu_tLepton_Pt","tEnu_Lepton_Pt","Enu_Lepton_Pt",
              "tEel_tLepton_Pt","Eel_tLepton_Pt","tEel_Lepton_Pt","Eel_Lepton_Pt"]:
        hmu_rw = GetSignalHist(numuMCfile,MU_SIGNALS,i)
        hmu_rw.Scale(numuDataPOT/numuMCPOT,"width")
        hmu_tw = GetSignalHist(numuTrueWeightedfile,MU_SIGNALS,i)
        hmu_tw.Scale(numuDataPOT/numuTrueWeightedPOT,"width")
        he = GetSignalHist(nueMCsignalfile,E_SIGNALS,i)
        he.Scale(numuDataPOT/nueMCsignalPOT,"width")
        Slicer = PlotTools.Make2DSlice
        PlotTools.MakeGridPlot(Slicer,DrawTWRW,[hmu_rw,hmu_tw,he],draw_seperate_legend = True)
        PlotTools.CANVAS.Print("{}{}_twrw.png".format(PLOTPATH,i))


    #ad-hoc
    for i in [("Eel_Lepton_Pt","tEel_Lepton_Pt"),
              ("Eel_tLepton_Pt","tEel_tLepton_Pt")]:
         hmu_en = GetSignalHist(numuMCrawfile,MU_SIGNALS,i[1])
         hmu_en.Scale(numuDataPOT/numuMCrawDataPOT,"width")
         hmu_r = GetSignalHist(numuMCrawfile,MU_SIGNALS,i[0])
         hmu_r.Scale(numuDataPOT/numuMCrawDataPOT,"width")
         Slicer = PlotTools.Make2DSlice
         PlotTools.MakeGridPlot(Slicer,DrawRecoVsTrue,[hmu_en,hmu_r],draw_seperate_legend = True)
         PlotTools.CANVAS.Print("{}{}_recovstrue.png".format(PLOTPATH,i[0]))



# MakeScaleFile2 makes reco weights)
# MakeScaleFile makes true weights


def MakeFluxWeight():
    filepath = "/cvmfs/minerva.opensciencegrid.org/minerva/CentralizedFluxAndReweightFiles/MATFluxAndReweightFiles/flux/flux-gen2thin-pdg{}{}-minervame{}_rearrangedUniverses.root"
    FHC = ("","1D1M1NWeightedAve")
    RHC = ("-","6A")
    Fout = ROOT.TFile.Open("flux_scale.root","RECREATE")
    for p in [FHC,RHC]:
        fe = ROOT.TFile.Open(filepath.format(p[0],12,p[1]))
        he = fe.Get("flux_E_cvweighted")
        he = he.Rebin(len(fluxbin)-1,"he_rebin",fluxbin)
        fmu = ROOT.TFile.Open(filepath.format(p[0],14,p[1]))
        hmu = fmu.Get("flux_E_cvweighted")
        hmu = hmu.Rebin(len(fluxbin)-1,"hmu_rebin",fluxbin)
        he.Divide(he,hmu)
        Fout.cd()
        he.Write("FHC" if not p[0] else "RHC")
    Fout.Close()

 
# def Convolution(hin):
#     # <w(E)> = int w(E_reco) * N(E-E_reco,sigma) dE_reco
#     # =sum w(E_reco) * N (E-E_reco,sigma) DelE_reco
#     sigma = 0.7 #GeV
#     h2 = hin.Clone()
#     for i in range(1,hin.GetSize()-1):
#         E_t = hin.GetBinCenter(i)
#         result = 0
#         f = ROOT.TF1("gau","TMath::Gaus(x,[0],[1],true)",0,20)
#         f.SetParameters(E_t,sigma)
#         for j in range(1,hin.GetSize()-1):
#             E_r = hin.GetBinCenter(j)
#             binwidth = hin.GetBinWidth(j)
#             result += hin.GetBinContent(j) * f.Integral(E_r-binwidth/2,E_r+binwidth/2)
#         h2.SetBinContent(i,result)
#     return h2

 
def Convolution(hin,migration):
    # <w(E)> = int w(E_reco) * N(E-E_reco,sigma) dE_reco
    # =sum w(E_reco) * N (E-E_reco,sigma) DelE_reco
    h2 = hin.Clone()
    for i in range(1,hin.GetSize()):
        E_t = hin.GetBinCenter(i)
        migration_row = migration.GetYaxis().FindBin(E_t)
        prob = migration.ProjectionX("prob",migration_row,migration_row)
        prob.Scale(1.0/prob.Integral(0,-1))
        #f = ROOT.TF1("gau","TMath::Gaus(x,[0],[1],true)",0,20)
        #f.SetParameters(E_t,sigma)
        result = 0
        for j in range(1,hin.GetSize()):
            #binwidth = hin.GetBinWidth(j)
            result += hin.GetBinContent(j) * prob.GetBinContent(j)
        h2.SetBinContent(i,result)
        del prob
    return h2


def MakeScaleFile3():
    Fout = ROOT.TFile.Open("emu_scale3.root","RECREATE")
    fflux = ROOT.TFile.Open("flux_scale.root")
    ratio = fflux.Get("FHC")
    Fout.cd()
    ratio.Write("FHC")
    tmp = ratio
    Fmu = ROOT.TFile.Open(numufile)
    migration = Fmu.Get("Enu_migration")
    for i in range(5):
        exp = Convolution(tmp,migration)
        Fout.cd()
        exp.Write("Exp_Enu_weight_{}".format(i+1))
        for i in range(1,exp.GetSize()-1):
            exp.SetBinContent(i,ratio.GetBinContent(i)*tmp.GetBinContent(i)/exp.GetBinContent(i))
        tmp=exp
    exp.Write("Enu")
    fflux.Close()
    Fout.Close()

if __name__ == "__main__":
    MakeCompPlot()
    #MakeScaleFile(["Enu","tEnu","tEnu_true_signal"])
    #GetEfficiencyCorrection()
    #MakeFluxWeight()
    #MakeScaleFile3()

