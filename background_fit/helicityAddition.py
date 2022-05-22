#python 2/3 compatible lines
from __future__ import print_function
if hasattr(__builtins__, 'raw_input'):
    input = raw_input

import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict
from tools.PlotLibrary import PLOT_SETTINGS,VariantPlotsNamingScheme, HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools


# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

fitPath = '/minerva/data/users/hsu/nu_e/'
nuebar_mcfile = 'kin_dist_mcmeRHC_wrong_sign_p_MAD.root'
nuebar_datafile = 'kin_dist_datameRHC_wrong_sign_p_MAD.root'
PLOTPATH = "/minerva/data/users/hsu/WrongSign/"

#nue_file = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_datafile =  'kin_dist_dataFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_mcfile = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)

# nue_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_file))
# data_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_datafile))

# nuedata_map_path = "{}/{}".format(fitPath,nue_datafile)
# nuemc_map_path = "{}/{}".format(fitPath,nue_mcfile)

HISTOGRAMS_TO_NUE_ADD = [
  "Visible Energy vs q3",
  "Visible Energy vs Lepton Pt",
#  "Lepton Energy",
#  "Lepton Theta",
#  "Q3",
#  "Visible Energy",
#  "Lepton Pt",
]


def BackgroundHistAddition(mc_hists): #used incase you want to add to existing hist

    out_h1 = mc_hists.hists["NCOther"].Clone() 
    out_h2 = mc_hists.hists["NuEElastic"].Clone() #nue events
    out_h3 = mc_hists.hists["NonPhaseSpace"].Clone()
    out_h4 = mc_hists.hists["NonFiducial"].Clone()
    out_h5 = mc_hists.hists["CCDIS"].Clone()
    out_h6 = mc_hists.hists["CCOther"].Clone()
    out_h7 = mc_hists.hists["NCCOH"].Clone()
    out_h8 = mc_hists.hists["NCHDIS"].Clone()
    out_h9 = mc_hists.hists["NCDIS"].Clone()
    out_h10 = mc_hists.hists["NCRES"].Clone()
    out_h11 = mc_hists.hists["CCNuEAntiNu"].Clone()
 
    out_h1.Add(out_h2) #add the nue events to the Other background
    out_h1.Add(out_h3)
    out_h1.Add(out_h4)
    out_h1.Add(out_h5)
    out_h1.Add(out_h6)
    out_h1.Add(out_h7)
    out_h1.Add(out_h8)
    out_h1.Add(out_h9)
    out_h1.Add(out_h10)
    out_h1.Add(out_h11)
    return out_h1


def SignalHistAddition(mc_hists):

    out_h1 = mc_hists.hists["CCNuEQE"].Clone()
    out_h2 = mc_hists.hists["CCNuEDelta"].Clone()
    out_h3 = mc_hists.hists["CCNuEDIS"].Clone()
    out_h4 = mc_hists.hists["CCNuE2p2h"].Clone()
    out_h5 = mc_hists.hists["CCNuE"].Clone()

    out_h1.Add(out_h2) #add the nue events to the Other background
    out_h1.Add(out_h3)
    out_h1.Add(out_h4)
    out_h1.Add(out_h5)

    return out_h1

def DataSubtraction(data_hist,nuebar_hist,background_hist,pot_scale):
    out_data = data_hist.Clone()
    out_nuebar = nuebar_hist.Clone()
    out_background = background_hist.Clone()  
     
    out_nuebar.AddMissingErrorBandsAndFillWithCV(nuebar_hist)
    out_background.AddMissingErrorBandsAndFillWithCV(background_hist)
    
    out_nuebar.Scale(pot_scale) #FHC data pot normalize
    out_background.Scale(pot_scale) #FHC data pot normalize
    out_data.AddMissingErrorBandsAndFillWithCV(background_hist) 
    out_data.Add(out_nuebar,-1)
    out_data.Add(out_background,-1)
   
    return out_data

def DataSubtractionFHC(data_hist,mc_hist,wrong_sign_hist,ref):
    data_hist.AddMissingErrorBandsAndFillWithCV(ref)
    mc_hist.AddMissingErrorBandsAndFillWithCV(ref)
    wrong_sign_hist.AddMissingErrorBandsAndFillWithCV(ref)
    data_hist.Add(mc_hist,-1)
    data_hist.Add(wrong_sign_hist)
    return data_hist

def DrawWrongSign(mnvplotter,data_hist,mc_hist):
    data_hist.SetLineColor(ROOT.kRed-2)
    data_hist.Draw("HIST")
    mc_hist.SetLineColor(ROOT.kBlue-2)
    mc_hist.Draw("HIST SAME")
    leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
    leg.AddEntry(data_hist,"Constrained wrong sign")
    leg.AddEntry(mc_hist,"GENIE wrong sign")
    ROOT.SetOwnership(leg,False)
    leg.Draw()
    

if __name__ == "__main__":
    
    #input knobs
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_reco_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)

    nuebar_files = []
    nuebar_pots = []
    for i in [nuebar_datafile,nuebar_mcfile]:
        path = "{}/{}".format(fitPath,i)
        nuebar_files.append(ROOT.TFile.Open(path))
        nuebar_pots.append(Utilities.getPOTFromFile(path))

    
    nuebar_pot_scale  = nuebar_pots[0] / nuebar_pots[1] #RHC_data_pot / FHC_MC_pot
    #data_pot_scale = nuepots[0] / nuepots[1] #FHC_data_pot / FHC_MC_pot needed for subtraction
    #datanue_pot_scale = pots[0] / nuepots[0] #RHC_data_pot / FHC_data_pot
    
    regions = ["Signal"] + AnalysisConfig.sidebands
    newplaylist=ROOT.TFile.Open('/minerva/data/users/hsu/nu_e/BkgWrongSignPred-meFHC.root',"RECREATE")
 
    for region in regions:
        for config in HISTOGRAMS_TO_NUE_ADD: 
            data_nuebar = HistHolder(config["name"] if "name" in config else config,nuebar_files[0],region,False,nuebar_pots[0],nuebar_pots[0])
            mc_nuebar = HistHolder(config["name"] if "name" in config else config,nuebar_files[1],region,True,nuebar_pots[1],nuebar_pots[0])
            ref = HistHolder(config["name"] if "name" in config else config,mc_reco_file,region,True,mc_pot,mc_pot)
            data_nuebar.POTScale(False)
            mc_nuebar.POTScale(False)
            pred = DataSubtractionFHC(data_nuebar.GetHist(),mc_nuebar.GetHist(),mc_nuebar.hists["CCNuEAntiNu"],ref.GetHist())
            newplaylist.cd()
            pred.Scale(mc_pot/nuebar_pots[0])
            pred.Write("{}_{}".format(pred.GetName(),"CCNuEAntiNu"))
            # background = BackgroundHistAddition(nue_hists)

            # nuebar = nue_hists.hists['CCNuEAntiNu']

            # nue_prediction =  DataSubtraction(data,nuebar,background,data_pot_scale)

            # #nue_prediction.Scale(1/nue_pot_scale) 
            # print(nue_prediction.GetBinContent(1,1), 'nue pred')
            # nue_prediction.Write("{}_OtherNue".format(nuebar.GetName())) #write as a new hist in the current tuple being used


            #drawing
            Slicer = PlotTools.Make2DSlice
            #plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
            PlotTools.MakeGridPlot(Slicer,DrawWrongSign,[pred,ref.hists["CCNuEAntiNu"]],draw_seperate_legend = True)
            PlotTools.CANVAS.Print("{}{}.png".format(PLOTPATH,pred.GetName()))
    if input("Warning: Will update mc file, continue? (y/N)").upper()!="Y":
        newplaylist.Close()
        exit(1)
    mc_reco_file.Close()
    mc_reco_file=ROOT.TFile.Open(type_path_map["mc"],"UPDATE")
    newplaylist.ReOpen("READ")
    region="Signal"
    for config in HISTOGRAMS_TO_NUE_ADD:
        ref = HistHolder(config["name"] if "name" in config else config,mc_reco_file,region,True,1)
        pred = newplaylist.Get(ref.plot_name+"_CCNuEAntiNu")
        ref.hists["CCNuEAntiNu"]=pred.Clone(ref.hists["CCNuEAntiNu"].GetName())
        print (ref.hists["CCNuEAntiNu"].GetName())
        ref.ResumTotal()
        mc_reco_file.cd()
        for i in ["Total","CCNuEAntiNu"]:
            h = ref.hists[i]
            h.Write(h.GetName(),ROOT.TObject.kOverwrite)
        print(ref.hists["CCNuEAntiNu"].GetBinContent(1,1))

    newplaylist.Close()
    mc_reco_file.Write()
    mc_reco_file.Close()
    data_file.Close()



