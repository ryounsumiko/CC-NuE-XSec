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
        if hist.HasVertErrorBand(universe.ShortName()):
            errorhist = hist.GetVertErrorBand(universe.ShortName()).GetHist(universe.ithInWrapper)
        else:
            errorhist = hist
        return errorhist.GetBinContent(hist.FindBin(v))

    def GetWeight(self,universe):
        if universe.nsigma is None:
            return 1
        if universe.classifier.truth_class in self.cate_map:
            return self.cate_map[universe.classifier.truth_class](universe)
        else:
            return 1

    def weightTruth(self):
        return False

class DataWeight(MyWeighterBase):
    def __init__(self,fvalue=None):
        super(DataWeight,self).__init__(fvalue)
        self.weighter = lambda v :1

    def GetWeight(self,universe):
        return self.weighter(universe)

class EelTuningWeight(MyWeighterBase):
    def __init__(self):
        super(EelTuningWeight,self).__init__( lambda universe:universe.kin_cal.reco_E_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/bkgfit_scale.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:self.f.Get(i) for i in ["DIS","Excess","NCCoh","Signal"]}
        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["DIS"])
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["DIS"])
        self.cate_map["ExcessModel"] = partial(self.fileBasedWeight,hist=self.hist_dict["Excess"])
        self.cate_map["NCCOH"] = partial(self.fileBasedWeight,hist=self.hist_dict["NCCoh"])
        self.cate_map["CCNuEQE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDelta"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE2p2h"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])

class RHCNewWeight(MyWeighterBase):
    def __init__(self):
        super(RHCNewWeight,self).__init__( lambda universe:universe.kin_cal.reco_E_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/bkgfit_NoScaleCollab.root".format(os.environ["CCNUEROOT"]))
        # self.fFHC = ROOT.TFile.Open("{}/background_fit/bkgfit_pi0Nrhc.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:self.f.Get(i) for i in ["Excess","NCCoh","Signal"]}
        # self.hist_dict["Signal"] = self.fFHC.Get("Signal")
        # self.hist_dict["DIS"] = self.fFHC.Get("Pi0")

        # self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["DIS"])
        # self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["DIS"])
        self.cate_map["ExcessModel"] = partial(self.fileBasedWeight,hist=self.hist_dict["Excess"])
        self.cate_map["NCCOH"] = partial(self.fileBasedWeight,hist=self.hist_dict["NCCoh"])
        #self.cate_map["CCNuEAntiNu"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        # self.cate_map["CCNuEQE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        # self.cate_map["CCNuEDelta"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        # self.cate_map["CCNuE2p2h"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        # self.cate_map["CCNuEDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        # self.cate_map["CCNuE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])

class EnuElectronMuonWeight(DataWeight):
    def __init__(self):
        super(EnuElectronMuonWeight,self).__init__(lambda universe:universe.kin_cal.reco_E_nu_cal)
        # self.f = ROOT.TFile.Open("{}/studies/flux_scale.root".format(os.environ["CCNUEROOT"]))
        # self.h = self.f.Get("FHC")
        self.f = ROOT.TFile.Open("{}/studies/emu_scale3.root".format(os.environ["CCNUEROOT"]))
        self.h = self.f.Get("Enu")
        self.h.ClearAllErrorBands()
        self.weighter = partial(self.fileBasedWeight,hist=self.h)


class EnuElectronMuonWeightTrue(DataWeight):
    def __init__(self):
        super(EnuElectronMuonWeightTrue,self).__init__(lambda universe:universe.kin_cal.true_enu_genie)
        self.f = ROOT.TFile.Open("{}/studies/flux_scale.root".format(os.environ["CCNUEROOT"]))
        self.h = self.f.Get("FHC")
        self.h.ClearAllErrorBands()
        self.weighter = partial(self.fileBasedWeight,hist=self.h)

    def weightTruth(self):
        return True

class EnuElectronMuonWeightRaw(DataWeight):
    def __init__(self):
        super(EnuElectronMuonWeightRaw,self).__init__(lambda universe:universe.kin_cal.reco_E_nu_cal)
        self.f = ROOT.TFile.Open("{}/studies/flux_scale.root".format(os.environ["CCNUEROOT"]))
        self.h = self.f.Get("FHC")
        self.h.ClearAllErrorBands()
        self.weighter = partial(self.fileBasedWeight,hist=self.h)

class RHCWrongSignEstimatorWeight(DataWeight):
    def __init__(self):
        super(RHCWrongSignEstimatorWeight,self).__init__(lambda universe:universe.kin_cal.reco_visE+universe.kin_cal.reco_E_lep)
        self.f = ROOT.TFile.Open("{}/studies/FHC_RHC_tEnu_Scale.root".format(os.environ["CCNUEROOT"]))
        self.h = self.f.Get("tEnu")
        self.weighter = partial(self.fileBasedWeight,hist=self.h)

# class SarahMagicWeight(MyWeighterBase):
#     def __init__(self):
#         super(SarahMagicWeight,self).__init__(lambda universe:universe.kin_cal.reco_E_lep)
#         RHC_Factors = {
#             "Range":[2.5,4,6,9,12,15,20],
#             "COH": [1,1.88,1.94,2.44,3.41,3.69,3.36,1],
#             "DFR": [1,3.56,5.91,6.25,8.68,4.41,0.004,1],
#         }
#         self.cate_map["ExcessModel"]=partial(self.rangeBasedWeight,ran=RHC_Factors["Range"],weight=RHC_Factors["DFR"])
#         self.cate_map["NCCOH"]=partial(self.rangeBasedWeight,ran=RHC_Factors["Range"],weight=RHC_Factors["COH"])

class PtTuningWeight(MyWeighterBase):
    def __init__(self):
        super(PtTuningWeight,self).__init__( lambda universe:universe.kin_cal.reco_Pt_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/tune2-Eel.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:PlotUtils.MnvH1D(self.f.Get(i).GetCVHistoWithStatError()) for i in ["Pi0"]}

        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        
class RHCPtTuningWeight(MyWeighterBase):
    def __init__(self):
        super(RHCPtTuningWeight,self).__init__( lambda universe:universe.kin_cal.reco_Pt_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/cv_tune-Eel.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:self.f.Get(i) for i in ["Pi0","Signal"]}

        #self.cate_map["ExcessModel"] = partial(self.fileBasedWeight,hist=self.hist_dict["Excess"])
        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCNuEQE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDelta"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE2p2h"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])

class WrongSignPtTuningWeight(MyWeighterBase):
    def __init__(self):
        super(WrongSignPtTuningWeight,self).__init__( lambda universe:universe.kin_cal.reco_Pt_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/bkgfit_EePtScaleCollab.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:self.f.Get(i) for i in ["Pi0"]}

        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])

class PtTuningWeightAlt(MyWeighterBase):
    def __init__(self):
        super(PtTuningWeightAlt,self).__init__( lambda universe:universe.kin_cal.reco_Pt_lep)
        self.f = ROOT.TFile.Open("{}/background_fit/tune1-Eel.root".format(os.environ["CCNUEROOT"]))
        self.hist_dict = {i:PlotUtils.MnvH1D(self.f.Get(i).GetCVHistoWithStatError()) for i in ["Excess","Pi0","Signal","NCCoh"]}

        self.cate_map["ExcessModel"] = partial(self.fileBasedWeight,hist=self.hist_dict["Excess"])
        self.cate_map["NCCoh"] = partial(self.fileBasedWeight,hist=self.hist_dict["NCCoh"])
        self.cate_map["NCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Pi0"])
        self.cate_map["CCNuEQE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDelta"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE2p2h"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuEDIS"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])
        self.cate_map["CCNuE"] = partial(self.fileBasedWeight,hist=self.hist_dict["Signal"])


class RHCEelPtTuningWeight(MyWeighterBase):
    def __init__(self):
        self.EelWeight = RHCNewWeight()
        self.PtWeight = RHCPtTuningWeight()

    def GetWeight(self,universe):
        return self.EelWeight.GetWeight(universe)*self.PtWeight.GetWeight(universe)

class FHCPtTuningWeight(MyWeighterBase):
    def __init__(self):
        self.CVWeight = RHCEelPtTuningWeight()
        self.FHCWeight = PtTuningWeight()

    def cut(self,universe):
        return 0.2<universe.kin_cal.reco_visE <0.5 and universe.kin_cal.reco_Pt_lep<1.0

    def GetWeight(self,universe):
        return self.CVWeight.GetWeight(universe)* (self.FHCWeight.GetWeight(universe) if self.cut(universe) else 1)

class FHCPtTuningWeightAlt(MyWeighterBase):
    def __init__(self):
        self.RHCWeight = RHCEelPtTuningWeight()
        self.FHCWeight = PtTuningWeightAlt()

    def GetWeight(self,universe):
        return self.RHCWeight.GetWeight(universe)*self.FHCWeight.GetWeight(universe)

class RHCWrongSignBkgWeight(MyWeighterBase):
    def __init__(self):
        self.Wrongsign = RHCWrongSignEstimatorWeight() #rhc->fhc weight
        self.RHCPtWeight = WrongSignPtTuningWeight() #pi^0 weight
        self.RHCEelWeight = RHCNewWeight() # Excess weight

    def GetWeight(self,universe):
        return self.Wrongsign.GetWeight(universe)* self.RHCPtWeight.GetWeight(universe)* self.RHCEelWeight.GetWeight(universe)


if AnalysisConfig.extra_weighter is None:
    MyWeighter = MyWeighterBase()
elif AnalysisConfig.extra_weighter =="Eel_tune":
    MyWeighter = EelTuningWeight()
elif AnalysisConfig.extra_weighter =="emu_weight":
    MyWeighter = EnuElectronMuonWeight()
elif AnalysisConfig.extra_weighter =="emu_true_weight":
    MyWeighter = EnuElectronMuonWeightTrue()
elif AnalysisConfig.extra_weighter =="emu_raw_weight":
    MyWeighter = EnuElectronMuonWeightRaw()
elif AnalysisConfig.extra_weighter =="rhc_weight":
    MyWeighter = RHCNewWeight()
elif AnalysisConfig.extra_weighter =="CV_tune":
    MyWeighter = RHCEelPtTuningWeight()
elif AnalysisConfig.extra_weighter =="FHCPt_tune2":
    MyWeighter = FHCPtTuningWeight()
elif AnalysisConfig.extra_weighter =="FHCPt_tune1":
    MyWeighter = FHCPtTuningWeightAlt()
elif AnalysisConfig.extra_weighter =="wrong_sign":
    MyWeighter = RHCWrongSignBkgWeight()

else:
    raise ValueError("Unknown extra weighter")
