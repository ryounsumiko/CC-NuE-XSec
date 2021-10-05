#!/usr/bin/env python
"""
Contains methods for getting weights from a event.
Includes flux, genie, minerva_genie_tunes(non-res-pion, rpa, and 2p2h)
"""
import os, ROOT, sys, uuid

from array import *
import math
from . import KinematicsCalculator
import PlotUtils
#from PlotUtils import FluxReweighter
import random
from functools import partial
from config.AnalysisConfig import AnalysisConfig

def getWeight(Ee, Theta, ETh, pc, channel):
    pcspline_hist = ROOT.TFile.Open("/minerva/data/users/shenry/antinu_e/kin_dist_datame5a6a_25mmExclusion_nx_PCPolynomialSF_new_fspline.root")

    Ee_hist = pcspline_hist.Get('Electron Energy_me6a_%s_FlatFinal_nx'%(pc))
    Theta_hist = pcspline_hist.Get('Electron Theta_me6a_%s_FlatFinal_nx'%(pc))
    ETheta_hist = pcspline_hist.Get('EeTheta_me6a_%s_FlatFinal_nx'%(pc))
    
    Ee_bin = Ee_hist.GetXaxis().FindBin(Ee)
    Ee_value = Ee_hist.GetBinContent(Ee_bin)
    Theta_bin = Theta_hist.GetXaxis().FindBin(Theta)
    
   
    Theta_value = Theta_hist.GetBinContent(Theta_bin)
    ETheta_bin = ETheta_hist.GetXaxis().FindBin(ETh) 
    ETheta_value = ETheta_hist.GetBinContent(ETheta_bin)

    Ee_pcval = Ee_value 
    Theta_pcval = Theta_value
    ETh_pcval = ETheta_value

    excessspline_hist = ROOT.TFile.Open("/minerva/data/users/shenry/antinu_e/kin_dist_datame5a6a_25mmExclusion_nx_ExcessSplineFits_new_fspline.root")

    Ee_Exhist = excessspline_hist.Get('Electron Energy')
    Theta_Exhist = excessspline_hist.Get('Electron Theta')
    ETheta_Exhist = excessspline_hist.Get('EeTheta') 
    
    EeEx_bin = Ee_Exhist.GetXaxis().FindBin(Ee)
    EeEx_value = Ee_Exhist.GetBinContent(EeEx_bin)
    ThetaEx_bin = Theta_Exhist.GetXaxis().FindBin(Theta)
    ThetaEx_value = Theta_Exhist.GetBinContent(ThetaEx_bin)
    EThetaEx_bin = ETheta_Exhist.GetXaxis().FindBin(ETh) 
    EThetaEx_value = ETheta_Exhist.GetBinContent(EThetaEx_bin)

    Ee_exval = EeEx_value
    Theta_exval = ThetaEx_value
    ETh_exval = EThetaEx_value

    sf_Ee = Ee_exval / Ee_pcval
    sf_Theta = Theta_exval / Theta_pcval
    sf_ETh = ETh_exval / ETh_pcval
   
    return eval('sf_%s'%(channel))

def getNormalization(pc, channel):
    weight_hist = ROOT.TFile.Open("/minerva/data/users/shenry/antinu_e/kin_dist_datame5a6a_25mmExclusion_nx_PCPolynomialSF_new_fspline.root")

    Ee_hist = weight_hist.Get('Electron Energy_me6a_%s_FlatFinal_nx'%(pc))
    Theta_hist = weight_hist.Get('Electron Theta_me6a_%s_FlatFinal_nx'%(pc))
    ETheta_hist = weight_hist.Get('EeTheta_me6a_%s_FlatFinal_nx'%(pc))

    Ee_normhist = Ee_hist.Clone()
    Theta_normhist = Theta_hist.Clone()
    ETheta_normhist = ETheta_hist.Clone()
    
    Ee_norm = 1/Ee_normhist.Integral()
    Theta_norm = 1/Theta_normhist.Integral()
    ETh_normhist = 1/Theta_normhist.Integral()
 
    return eval('%s_norm'%(channel))


def GetModelWeight(event,channel,pc):
    Ee = event.kin_cal.reco_E_e 
    theta = event.kin_cal.reco_theta_e_rad
    Theta = event.kin_cal.reco_theta_e
    Eth = Ee*theta #Gev*rad 
    
    return getWeight(Ee, Theta, Eth, pc, channel)#, getNormalization(pc, channel) #returns scale fns for Ee, theta, eth

class CustomizedWeighter(object):
    def __init__ (self):
        RHC_Factors = {
            "Range":[2.5,4,6,9,12,15,20],
            "COH": [1,1.88,1.94,2.44,3.41,3.69,3.36,1],
            "DFR": [1,3.56,5.91,6.25,8.68,4.41,0.004,1],
        }
        FHC_Factors1 = {
            "Range":[2.5,4,6,9,12,15,20],
            "COH": [1,1.96,2.02,1.76,1.18,1.46,1.46,1.54],
            "DFR": [1,5.59,9.04,4.68,1.06,1,1,1],
        }
        FHC_Factors = {
            "Range":[2.5,4,6,9,12,15,20],
            "COH": [1,1.51,1.94,2.57,4.91,5.12,1.01,1],
            "DFR": [1,5.33,12.0,14.8,20.7,24.4,12.23,1],
        }
        Default = {
            "Range":[100],
            "COH": [1],
            "DFR": [1],
        }
        self.weight = Default

    def GetWeight(self,universe):
        return self.GetCohWeight(universe)*self.GetDFRWeight(universe)

    def GetCohWeight(self,universe):
        if universe.mc_current == 2 and universe.mc_intType == 4:
            Eel = universe.kin_cal.reco_E_lep
            for i in range(len(self.weight["Range"])):
                if Eel <= self.weight["Range"][i]:
                    return self.weight["COH"][i]
                else:
                    continue
            return self.weight["COH"][-1]
        else:
            return 1

    def GetDFRWeight(self,universe):
        if universe.mc_intType==10:
            Eel = universe.kin_cal.reco_E_lep
            for i in range(len(self.weight["Range"])):
                if Eel <= self.weight["Range"][i]:
                    return self.weight["DFR"][i]
                else:
                    continue
            return self.weight["DFR"][-1]
        else:
            return 1

class MyWeighterBase(object):
    def __init__(self,fvalue = None):
        self.fvalue = fvalue
        self.cate_map={}

    def rangeBasedWeight(self,universe,ran,weight):
        v =self.fvalue(universe)
        if v is None:
            return 1
        for i in range(len(ran)):
            if v <= ran[i]:
                return weight[i]
            else:
                continue
        return weight[-1]

    def fileBasedWeight(self,universe,hist):
        v = self.fvalue(universe)
        if v is None:
            return 1
        if universe.ShortName() != "cv":
            errorhist = hist.GetVertErrorBand(universe.ShortName()).GetHist(universe.ithInWrapper)
        else:
            errorhist = hist
        return errorhist.GetBinContent(hist.FindBin(v))

    def GetWeight(self,universe):
        if universe.classifier.truth_class in self.cate_map:
            return self.cate_map[universe.classifier.truth_class](universe)
        else:
            return 1

class FileBasedWeighter(MyWeighterBase):
    def __init__(self,file_path,fvalue=None):
        super(FileBasedWeighter,self).__init__(fvalue)
        self.f = ROOT.TFile.Open(file_path)

class EelTuningWeight(FileBasedWeighter):
    def __init__(self):
        super(EelTuningWeight,self).__init__("{}/background_fit/bkgfit_scale.root".format(os.environ["CCNUEROOT"]), lambda universe:universe.kin_cal.reco_E_lep)
        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.f.Get("DIS"))
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.f.Get("DIS"))
        self.cate_map["ExcessModel"] = partial(self.fileBasedWeight,hist=self.f.Get("Excess"))
        self.cate_map["NCCOH"] = partial(self.fileBasedWeight,hist=self.f.Get("NCCoh"))
        self.cate_map["CCNuEQE"] = partial(self.fileBasedWeight,hist=self.f.Get("Signal"))
        self.cate_map["CCNuEDelta"] = partial(self.fileBasedWeight,hist=self.f.Get("Signal"))
        self.cate_map["CCNuE2p2h"] = partial(self.fileBasedWeight,hist=self.f.Get("Signal"))
        self.cate_map["CCNuEDIS"] = partial(self.fileBasedWeight,hist=self.f.Get("Signal"))
        self.cate_map["CCNuE"] = partial(self.fileBasedWeight,hist=self.f.Get("Signal"))

if AnalysisConfig.extra_weighter is None:
    MyWeighter = MyWeighterBase()
elif AnalysisConfig.extra_weighter =="Eel_tune":
    MyWeighter = EelTuningWeight()
else:
    raise ValueError("Unknown extra weighter")
