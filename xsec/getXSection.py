"""
Get Crossection Plots.
Inputs are unfolded data histogram, efficiency histogram(or seperated true signal and selected signal),
"""

import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS,AnaNuPDG
from config import DrawingConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE
from functools import partial

ROOT.TH1.AddDirectory(False)

XSEC_TO_MAKE = [
    #"q0 vs q3",
    "Visible Energy vs q3",
    "Visible Energy vs Lepton Pt"
]
USE_BIGNUE = True
threshold = 100 if USE_BIGNUE else 1
TARGET_UTILS = PlotUtils.TargetUtils.Get()
warping_errorband = ["fsi_weight","SuSA_Valencia_Weight","MK_model","LowQ2Pi_Joint","LowQ2Pi_NUPI0","LowQ2Pi_None"]
#warping_errorband = ["SuSA_Valencia_Weight"]
FLUX="minervame1d1m1nweightedave"
#FLUX="minervame1d"
COLORS=ROOT.MnvColors.GetColors()
MODELS = {
    "MnvTune v2": {
        "errorband":(None,None),
        "color":COLORS[0]
    },
    "2p2h Tune (QE)": {
        "errorband":("Low_Recoil_2p2h_Tune",2),
        "color":COLORS[1]
    },
    "SuSA 2p2h" : {
        "errorband":("SuSA_Valencia_Weight",0),
        "color":COLORS[2]
    },
    "MK model": {
        "errorband": ("MK_model",0),
        "color":COLORS[7]
    },
    # "Low Q2 Pion Joint": {
    #     "errorband" : ("LowQ2Pi_Joint",0),
    #     "color": COLORS[5]
    # },
    "MnvTune v1": {
        "errorband" : ("LowQ2Pi_None",0),
        "color": COLORS[6]
    },
    "FSI bug fix": {
        "errorband" : ("fsi_weight",0),
        "color": COLORS[5]
    },

}

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    #divide by efficiency
    efficiency.AddMissingErrorBandsAndFillWithCV(unfolded)
    unfolded.Divide(unfolded,efficiency)
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N nucleaon
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],is_mc)
    h_nucleon = GetNnucleonError(unfolded,Nnucleon)
    unfolded.Divide(unfolded,h_nucleon)
    return unfolded

def GetNnucleonError(hist,ntargets):
    hist_target = hist.Clone("number_of_targets")
    hist_target.ClearAllErrorBands()
    hist_target.Reset()
    errband_name = "Target_Mass_CH"
    band_err = 0.014
    hist_target.AddVertErrorBand(errband_name,2)

    for i in range(hist_target.GetSize()):
        hist_target.SetBinContent(i,ntargets)
        hist_target.SetBinError(i,0)
        hist_target.GetVertErrorBand(errband_name).SetBinContent(i,ntargets)
        hist_target.GetVertErrorBand(errband_name).GetHist(0).SetBinContent(i,ntargets*(1-band_err))
        hist_target.GetVertErrorBand(errband_name).GetHist(1).SetBinContent(i,ntargets*(1+band_err))
    hist_target.AddMissingErrorBandsAndFillWithCV(hist)
    print ("Target Normalization: {:.4E},{:.4E}".format(ntargets,ntargets*0.014))
    return hist_target


def DivideFlux(unfolded,is_mc):
    frw= PlotUtils.flux_reweighter(FLUX,AnaNuPDG,USE_NUE_CONSTRAINT) #playlist is dummy for now
    flux = frw.GetIntegratedFluxReweighted(AnaNuPDG,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    flux.PopVertErrorBand("Flux_BeamFocus")
    flux.PopVertErrorBand("ppfx1_Total")
    flux.Scale(1e-4*(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
    print ("Flux Normalization: {:.4E},{:.4E}".format(flux.GetBinContent(1,1),flux.GetTotalError(False).GetBinContent(1,1)))
    unfolded.Divide(unfolded,flux)

def ProjectUniverseContent(num, den, i, j,num_o,den_o):
    l = 0
    denNewContent = 0
    numNewContent = 0
    while i-l>=0:
        denNewContent += den_o.GetBinContent(i-l,j)
        numNewContent += num_o.GetBinContent(i-l,j)
        l+=1
        if denNewContent>threshold:
            break

    den.SetBinContent(i,j,denNewContent)
    num.SetBinContent(i,j,numNewContent)
    #need to do this for every universe
    for bandname in den.GetVertErrorBandNames():
        nHists = den.GetVertErrorBand(bandname).GetNHists()
        for k in range(0,nHists):
            denNewContent = sum (den_o.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-_,j) for _ in range(l))
            numNewContent = sum (num_o.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-_,j) for _ in range(l))
            den.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,denNewContent)
            num.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,numNewContent)
    return  i-l>=0

def ProjectBinContent(num, den):
    num_new = num.Clone("{}_smooth".format(num.GetName()))
    den_new = den.Clone("{}_smooth".format(den.GetName()))
    xNBins = num.GetNbinsX()
    yNBins = num.GetNbinsY()
    for j in range(0,yNBins+2):
        for i in range(0,xNBins+2):
            if 0 < den.GetBinContent(i,j) < threshold: #treshold
                if not ProjectUniverseContent(num_new, den_new, i, j,num,den):
                    print ("not enough events in {}-th y bin".format(j))
    return num_new, den_new

def GetEfficiency(ifile, ifile_truth, plot):
    if ifile_truth:
        print ("using separate efficiency file")
        ifile = ifile_truth
    num = Utilities.GetHistogram(ifile,PLOT_SETTINGS[plot+" Migration"]["name"]+"_truth")
    den = Utilities.GetHistogram(ifile,PLOT_SETTINGS["True Signal "+plot]["name"])
    hist_out,den_new = ProjectBinContent(num,den)
    hist_out.Divide(hist_out,den_new,1.0,1.0,"B")
    #plot efficiency
    hist_out.GetYaxis().SetTitle(den.GetYaxis().GetTitle())
    hist_out.GetXaxis().SetTitle("E_{avail} (GeV)") #hard coding for now
    hist_out.GetZaxis().SetTitle("Efficiency")

    plotter = lambda mnvplotter, hist: mnvplotter.DrawMCWithErrorBand(hist.GetCVHistoWithError())
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    #PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist)
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff_err"))
    return hist_out

def GetMCXSectionHistogram(mc_file,plot):
    # mc unfolded is truth_signal:
    unfolded = Utilities.GetHistogram(mc_file,PLOT_SETTINGS["True Signal "+plot]["name"])
    DivideFlux(unfolded,True)
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],True)
    unfolded.Scale(1.0/Nnucleon)
    sig_dep = []
    colors = []
    titles = []
    for chan in DrawingConfig.SignalDecomposition.keys():#["CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"]:
        if chan == "Background":
            continue
        tmp =Utilities.GetHistogram(mc_file,"{}_{}".format(PLOT_SETTINGS["True Signal "+plot]["name"],chan))
        DivideFlux(tmp,True)
        tmp.Scale(1.0/Nnucleon)
        sig_dep.append(tmp)
        colors.append(DrawingConfig.SignalDecomposition[chan]["color"])
        titles.append(DrawingConfig.SignalDecomposition[chan]["title"])
    return unfolded,sig_dep,colors,titles

def DrawModelComparison(data_hist,mc_hist,models=MODELS,band_on_mc=True):
    _cate = []
    _mc_models = []
    _colors = []
    for k,v in models.items():
        _cate.append(k)
        _colors.append(v["color"])
        htmp = mc_hist if band_on_mc else data_hist
        if v["errorband"][0]:
            try:
                _mc_models.append(PlotUtils.MnvH2D(htmp.GetVertErrorBand(v["errorband"][0]).GetHist(v["errorband"][1])))
            except ReferenceError:
                continue
        else:
            _mc_models.append(PlotUtils.MnvH2D(htmp.GetCVHistoWithStatError()))
    plotter = lambda mnvplotter,data_hist, *mc_ints : partial(PlotTools.MakeModelVariantPlot,color=_colors,title=_cate)(data_hist, mc_ints)
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[data_hist,*_mc_models],draw_seperate_legend=True)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_models"))
    # data_hist.GetZaxis().SetTitle("Data/MC")
    # mc_hist.GetZaxis().SetTitle("Data/MC")
    h0 = PlotUtils.MnvH2D(mc_hist.GetCVHistoWithStatError())
    for i in _mc_models:
        i.Divide(i,h0)
        i.GetZaxis().SetTitle("Data/MC")
    h0.AddMissingErrorBandsAndFillWithCV(data_hist)
    data_hist.Divide(data_hist,h0)
    data_hist.GetZaxis().SetTitle("Data/MC")
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[data_hist,*_mc_models],draw_seperate_legend=True)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_models_ratio"))

def DrawClosure(mnvplotter,data_hist,mc_hist):
    mnvplotter.DrawDataMCWithErrorBand(data_hist.GetCVHistoWithError(),mc_hist.GetCVHistoWithStatError(),1.0,"TR")
    leg = PlotTools.GetTLegend(ROOT.gPad)
    l=leg.GetListOfPrimitives();
    l[0].SetLabel("signal rich sample")
    l[1].SetLabel("standard sample")
    leg.Draw()

def setAxisTitles(h,plot):
    xlabel = "E_{avail}"
    ylabel = "q_{3}" if plot.split()[-1]=="q3" else "P^{t}_{lep}"
    zlabel = "d^{{2}}#sigma/d{X}d{Y}".format(X=xlabel,Y=ylabel)
    zlabel += " (#times 10^{-39} cm^{2}/GeV^{2})"
    h.GetXaxis().SetTitle("{} (GeV)".format(xlabel))
    h.GetYaxis().SetTitle("{} (GeV)".format(ylabel))
    h.GetZaxis().SetTitle(zlabel)

def DrawErrors(xsec,plot):
    plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist)
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
    PlotTools.MNVPLOTTER.axis_maximum = 1
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err"))
    # break downs
    PlotTools.MNVPLOTTER.axis_maximum = 0.35
    plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist,"TR",False,True,0.00001,False,"Detector model")
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err_det"))


    #     PlotTools.AdaptivePlotterErrorGroup(xsec,["GENIE_MaCCQE", "GENIE_MaNCEL",
    # "GENIE_MaRES",
    # "GENIE_MvRES",])
    PlotTools.updatePlotterErrorGroup(DETAILED_ERROR_GROUPS)
    for i in ["GENIE","GENIE-FSI","MnvTunes"]:
        plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist,"TR",False,True,0.00001,False,i)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err_{}".format(i)))

    PlotTools.MNVPLOTTER.axis_maximum = -1111

if __name__ == "__main__":
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_reco_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,AnalysisConfig.bkgTune_tag,False))
    #mc_reco_file = ROOT.TFile.Open(AnalysisCmc_reco_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist,False,False))
    mc_truth_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist+"-BigNuE",False)) if USE_BIGNUE else None 
    xsec_file = ROOT.TFile.Open(AnalysisConfig.XSecHistoPath(playlist),"RECREATE")
    for plot in XSEC_TO_MAKE:
        unfolded = Utilities.GetHistogram(unfolded_file,PLOT_SETTINGS[plot]["name"]+"_bkg_unfolding")
        efficiency = GetEfficiency(mc_reco_file,mc_truth_file,plot)
        xsec = GetXSectionHistogram(unfolded,efficiency,False)
        
        mc_xsec,sig_dep,colors,titles = GetMCXSectionHistogram(mc_reco_file,plot)
        setAxisTitles(xsec,plot)
        setAxisTitles(mc_xsec,plot)
        # ylabel = "q_{3}" if plot.split()[-1] == "q3" else
        # ylabel = xsec.GetZaxis().GetTitle().replace("NEvents","d^{{2}}#sigma/d{X}d{Y}".format(X="E_{avail}",Y=plot.split()[-1]))
        # ylabel += "(#times 10^{39} cm^{2}/GeV^{2}/Nucleon)"
        # xsec.GetZaxis().SetTitle(ylabel)
        # mc_xsec.GetZaxis().SetTitle(ylabel)
        # mc_xsec.GetXaxis().SetTitle("Eavail (GeV)") #hard coding for now

        for h in [xsec,mc_xsec,*sig_dep]:
            h.Scale(1e39,"width")
        for e in warping_errorband:
            xsec.PopVertErrorBand(e)

        xsec_file.cd()
        xsec.Write(PLOT_SETTINGS[plot]["name"]+"_dataxsec")
        mc_xsec.Write(PLOT_SETTINGS[plot]["name"]+"_mcxsec")

        plotter = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist.GetCVHistoWithError(True),mc_hist.GetCVHistoWithStatError(),1.0,"TR")
        #plotter = DrawClosure
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec"))
        plotter = lambda mnvplotter, data_hist,mc_hist : mnvplotter.DrawDataMCRatio(data_hist.GetCVHistoWithError(),mc_hist.GetCVHistoWithStatError(),1.0,True,0.9,1.1)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_ratio"))

        DrawErrors(xsec,plot)
       

        plotter = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(PlotTools.MakeSignalDecomposePlot,color=colors,title=titles)(data_hist,PlotUtils.MnvH1D(mc_hist.GetCVHistoWithStatError()),mc_ints)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec,*sig_dep],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_sigdep"))
        DrawModelComparison(xsec,mc_xsec,band_on_mc=True)
        print("done")
    xsec_file.Close()
