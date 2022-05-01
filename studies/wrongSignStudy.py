import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities
from functools import partial

ROOT.TH1.AddDirectory(False)
FHCfile = "/minerva/data/users/hsu/nu_e/kin_dist_mcmeFHC_wrong_sign_MAD.root"
RHCfile = "/minerva/data/users/hsu/nu_e/kin_dist_mcmeRHC_wrong_sign_MAD.root"
hist_name = "tEnu_CCNuEAntiNu"

def getFileAndPOT(filename):
    POT = Utilities.getPOTFromFile(filename)
    f = ROOT.TFile.Open(filename)
    return f,POT

def makeRatio(hist_name):
    Fout = ROOT.TFile.Open("FHC_RHC_tEnu_Scale.root","RECREATE")
    f_fhc,POT_fhc = getFileAndPOT(FHCfile)
    f_rhc,POT_rhc = getFileAndPOT(RHCfile)
    hFHC = f_fhc.Get(hist_name)
    hRHC = f_rhc.Get(hist_name)
    hRHC.Scale(POT_fhc/POT_rhc)
    hFHC.Divide(hFHC,hRHC)
    Fout.cd()
    hFHC.Write("tEnu")
    Fout.Close()
    f_fhc.Close()
    f_rhc.Close()

if __name__ =="__main__":
    makeRatio(hist_name)


