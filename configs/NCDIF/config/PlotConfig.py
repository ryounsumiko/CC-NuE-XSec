
"""
  PlotConfig.py:
   Centralization of common plotting configurations.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2014
"""

import math
import ROOT
from collections import OrderedDict

#from AnalysisConfig import AnalysisConfig

# constants needed in calculations.
M_e = 0.511 # MeV
M_p = 938.272 # MeV
M_n = 939.565 # MeV
QE_binding_E = 34 # MeV  (same as nu_mu PRL for FHC)

# these will be used repeatedly in calculations.
# better just do the arithmetic once here.
M_n_star = M_n - QE_binding_E
M_e_sqr = M_e**2
M_n_star_sqr = M_n_star**2
M_p_sqr = M_p**2

NQ3 = 8
NQ0 = 19
LOW_RECOIL_BIN_Q3 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
BACKGROUND_FIT_Q3_BIN = [0.0, 0.6, 0.8, 1.0, 1.2, 1.6, 2]
LOW_RECOIL_BIN_Q0 = [0.0, 0.04, 0.08,
		             0.12, 0.16, 0.2,
		             0.25, 0.30, 0.35, 0.40,
		             0.50, 0.60, 0.80, 1.00, 1.2]

PT_BINNING = [0.2*i for i in range(9)]

RESOLUTION_BINNING = [-1.0+0.1*i for i in range(41)]




MC_EM_ENERGY_SCALE = 1.05  # this from Trung: Doc 9370, slides 2-5, 15; Doc 9911; Doc 10102, slides 18-20
MC_MICHEL_ELECTRON_ENERGY_SCALE = 1.03  # from the NIM
MC_DEDX_ENERGY_SCALE = 1.015  # Based on eyeballing the dE/dx distribution for this analysis in the extra energy sideband.


# WARNING -- change this at your peril.  it affects the accepted region for true events
# via the 'truth' hacks in SelectedSamplePlots,
# which consequently affects the range that the efficiency correction corrects into!
"""
if AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.CCQE_LIKE:
	NEUTRINO_ENERGY_RANGE = [0, 10] # in GeV.  
elif AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.GENIE_CCQE:
	NEUTRINO_ENERGY_RANGE = [1.5, 10]
elif AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.EXCESS:
	NEUTRINO_ENERGY_RANGE = [0, float("inf")]  # this is presumptively an NC process, so anything's fair game
"""

# used in making plots of various particle types
PID_LABELS = {
	0: "other",
	11: "e^{#pm}",
	13: "#mu^{#pm}",
	22: "#gamma",
	111: "#pi^{0}",
	211: "#pi^{#pm}",
	2212: "p^{+}",
}

# for plotting atomic number.
# want bins for H, C, O, Si, Cl, Ti, Fe, Pb 
A_LABELS = {
	0: "other",
	1: "H",
	12: "C",
#	14: "N",
	16: "O",
	27: "Al",
	28: "Si",
	35: "Cl",
	48: "Ti",
	56: "Fe",
	207: "Pb", 
}

#ELECTRON_ANGLE_BINNING =  range(10) + [10, 12, 15, 20, 27, 35]
ELECTRON_ANGLE_BINNING = [2*i for i in range(25)] 
#ELECTRON_ANGLE_BINNING =  [1.3 * i for i in range(32)]
ELECTRON_ANGLE_RESIDUAL_BINNING =  [-1.+0.1* i for i in range(0,21)]
#ELECTRON_ANGLE_2D_BINNING = [-x for x in reversed(ELECTRON_ANGLE_BINNING)] + ELECTRON_ANGLE_BINNING[1:]
#EXCESS_ANGLE_BINNING = range(0, 15, 3) + [15, 20, 27, 35]
EXCESS_ANGLE_BINNING = list(range(0, 15, 3)) + [20]
LOW_RECOIL_BIN_HIGH_Q0 = [0.2, 0.25, 0.30, 0.35, 0.40,
                     0.50, 0.60, 0.80, 1.00, 1.5]
LOW_RECOIL_BIN_LOW_Q0 = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1,
                     0.12, 0.14, 0.16, 0.2]
NUECONE_BINNING = [10* i for i in range(0,51)]

#Changed ee binning because it was going up too high
#ELECTRON_ENERGY_BINNING_VS = [1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 10, 14, 24]
ELECTRON_ENERGY_BINNING_VS = [1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 10, 14, 24]
ELECTRON_ENERGY_BINNING = [1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 10, 14, 24]
#ELECTRON_ENERGY_BINNING = [0.5 * i for i in range(30)]
VISIBLE_ENERGY_RESIDUAL_BINNING = [-1+0.04* i for i in range(0,51)]
ELECTRON_ENERGY_RESIDUAL_BINNING = [-1+0.05* i for i in range(0,41)]
#EXCESS_ENERGY_BINNING = [0, 3, 6, 9, 12, 15, 20]  # ELECTRON_ENERGY_BINNING + [12, 15, 20]
#ELECTRON_ENERGY_BINNING = [0.75, 2, 3, 5, 7, 9, 20]  # Jaewon's bins
#NEUTRINO_ENERGY_BINNING = [i for i in range(6)] + [7, 10] # , 13, 18, 25]
NEUTRINO_ENERGY_BINNING = range(0,21)
#NEUTRINO_ENERGY_BINNING_BIGGER = [i for i in range(6)] + [7, 10, 13, 18, 25]

VISIBLE_ENERGY_BINNING = [ 0.02*i for i in range(9) ] + [0.05*i for i in range(4,9)]+ [0.1*i for i in range(5, 7)] + [0.2*i for i in range(4, 8)] 
APOTHEM_BINNING = [50*i for i in range(17)]
VERTEX_Y_BINNING = [100*i -1000 for i in range(21)]

#ETH2_BINNING = [0, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05]
ETH2_BINNING = [0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.5]

QSQUARED_BINNING_CCQE_LIKE = [0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0, 1.25, 2.0]
QSQUARED_BINNING_GENIE_CCQE = [0, 0.1, 0.2, 0.4, 0.8, 1.2, 2.0]
QSQUARED_BINNING_INC = [0.1 *i for i in range(21) ]
EXCESS_ENERGY_BINNING = [0, 3, 6, 9, 12, 15, 20]
SUM_VISIBLE_ENERGY_BINNING = [1 * i for i in range(1,11)]
RECO_W_BINNING = [0.1 * i for i in range(51)]
VERTEX_Z_BINNING = [5000+ 200*i for i in range(21)]

BACKGROUND_FIT_Q3_BIN =[0.0, 0.6, 0.8, 1.0, 1.2, 1.6, 2]
LOW_RECOIL_BIN_Q3_Truth = [0.0, 0.2, 0.4,0.6, 0.8, 1.0, 1.2]
LOW_RECOIL_BIN_Q3 =  LOW_RECOIL_BIN_Q3_Truth
LOW_RECOIL_BIN_Q0_Truth = [0,0.05]+[ 0.1*i for i in range(1,7)]+[0.2*i for i in range(4,7)]
#LOW_RECOIL_BIN_Q0 = [0.05 * i for i in range(17)] + [0.1 * i for i in range(9,21)]#,0.05,0.1,0.12,0.16,0.24,0.32,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
#LOW_RECOIL_BIN_Q0 = [0.05 * i for i in range(25)]#,0.05,0.1,0.12,0.16,0.24,0.32,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
# LOW_RECOIL_BIN_Q0_Truth = [0.0, 0.08, 0.16,0.32,0.60, 1.00, 1.2]
LOW_RECOIL_BIN_Q0 = [0.0,0.04,0.08,0.12,0.16,0.24,0.32,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
#LOW_RECOIL_BIN_Q0 =  [0.04*i for i in range(31)]

PT_BINNING_Truth = PT_BINNING

#QSQUARED_BINNING = QSQUARED_BINNING_CCQE_LIKE if AnalysisConfig.signal_defn != AnalysisConfig.SIGNAL_DEFNS.GENIE_CCQE else QSQUARED_BINNING_GENIE_CCQE
QSQUARED_BINNING = QSQUARED_BINNING_INC

QSQUARED_SLICE_EDGES = [0.2, 0.75, 2,]

W_BINNING = [0.5*i for i in range(21)]

PSI_BINNING = [0.01*i for i in range(40)]
#PSI_BINNING2 = [0.04*i for i in range(27)]
PSI_BINNING2 = [0.01*i for i in range(30)] #Ellie's binning for psi vs Ee
PSI_TAIL_BINNING = [0.05*i for i in range(4,15)]
PSI_FRONT_BINNING = [0.01*i for i in range(0,11)]
BIN_AVAIL_OPTIMIZATION = [0.0, 0.080,0.260,0.540,1.2,2.0]
TRANS_BINNING = [0.1*i for i in range(0,21)]
LONG_BINNING = [0.1*i for i in range(-20,51)]
LONG_2D_BINNING = [-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]
LONG_NEG_BINNING = [0.1*i for i in range(-20,1)]
LONG_POS_BINNING = [0.1*i for i in range(0,51)]

INLINE_UPSTREAM_ENERGY_BINS_VS = [5*i for i in range(2, 25)]
INLINE_UPSTREAM_ENERGY_BINS = [0, 10, 15, 25, 35, 50, 150]
UPSTREAM_INLINE_ENERGY_BINS_LOG = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500] 
UPSTREAM_INLINE_ENERGY_BINS_VS_VTEX = [10 * i for i in range (0, 13)]

VERTEX_ENERGY_BINS = [0, 5, 10, 15, 20, 30, 40,] + list(range(50, 350,30))

WHOLE_DET_Z_BINS = range(4000, 10000, 100)

VERTEX_Z_RESOLUTION_BINNING = [0.1*i for i in range(41)]


BIRK_MIGRATION_PID_CV_BINNING=[0.49 + 0.05*i for i in range(0,11)]
BIRK_MIGRATION_PID_SHIFTED_BINNING=[0.49 + 0.05*i for i in range(0,11)]

BIRK_MIGRATION_DEDX_BINNING=[2+0.05* i for i in range(0,21)]

FRACTION_BINNING = [0.02 * i for i in range(51) ]

LLR_BINNING = [i for i in range(-20,21)]

T_BINNING_FINE = [0.03*i for i in range(0, 16)]
T_BINNING = [0.06*i  for i in range(0, 8)]
T_BINNING_VS = [0, 0.05, 0.1, 0.15,  0.2, 0.25, 0.3, 0.35, 0.4]
T_BINS_MIG_VS = [0.025*i for i in range(0, 26)]

T_BINS = [0, 0.04, 0.07, 0.1, 0.15, 0.5]

ENERGY_BALANCE_BINNING = [0.02* i for i in range(26)]

LOG_ENERGY_BINNING = [0.1*i for i in range(41)]
FRONTDEDX_BINNING= [0.4 * i for i in range(26)]
FRONTDEDX_TAIL_BINNING = [0.4 * i for i in range(6,26)]
FRONTDEDX_FRONT_BINNING = [0.4 * i for i in range(7)]
#FRONTDEDX_BINNING = [0.2 * i for i in range(51)]
FRONTDEDX_POSITION_BINNING = [5 * i for i in range(-3, 83)]
DEDX_BINNING = [0.2 * i for i in range(26)]
FRONTDEDX_MEDIAN_BINNING = [0.4 * i for i in range(26)]
FRONTDEDX_2D_BINNING = [0,2.401,10]

PION_ENERGY_BINNING = [0.4 * i for i in range(26)]
PROTON_ENERGY_BINNING = [0.1 * i for i in range(21)]
PROTON_KE_BINNING = [0.05 * i for i in range(21)]
#VERTEX_DIFF_BINNING = range(-20,40)
VERTEX_DIFF_BINNING = [10*i for i in range(-5,40)]
VERTEX_DIFF_BINNING_VS_INLINE = [45*i for i in range(-1, 10)]
PT_BINNING = [0.07 * i for i in range(19)]
PI0_ENERGY_BINNING = [0.5 * i for i in range(21)]
#EXTRA_ENERGY_BINNING = [0.05*i for i in range(41)]
EXTRA_ENERGY_BINNING = [0.05*i for i in range(15)]
OUTSIDE_ENERGY_BINNING = [0.02*i for i in range(30)]
#OUTSIDE_ENERGY_BINNING = [0.1*i for i in range(50)]

#ADDED BELOW TO FIX ERROR
PROTON_ANGLE_BINNING = [0.4 * i for i in range(26)]
PROTON_ANGLE_RESIDUAL_BINNING = [0.4 * i for i in range(26)]
APOTHEM_BINNING = [50*i for i in range(17)]

PARTICLE_PDG_BINNING = range(-11,10)
GENIE_EVENT_BINNING = [0.5 * i for i in range(0,17)]
FRACTION_BINNING = [0.1 * i for i in range(-0,12)]
CLUSTER_PDG_FRAC = [0.05 * i for i in range(0,23)]
CLUSTER_PDG = [5* i for i in range(-45, 45)]
CLUSTER_PDG_HIGH = [5* i for i in range(400, 460)]
#CLUSTER_PDG_FULL = [5*i for i in range(-6,460)]
CLUSTER_PDG_FULL = [5*i for i in range(-50,460)]

POSITION_BINNING_LOW = [0, 25.01, 50.01, 75.01, 100.01, 125.01, 150.01, 175.01, 200.01, 225.01]
POSITION_BINNING_HIGH = [225.01, 275.01, 300.01, 325.01, 350.01, 375.01, 400.01]
POSITION_BINNING_AVG = [0, 100.01, 200.01, 300.01, 400.01]
DIFFERENCE_BINNING = [0.01 * i for i in range(-20, 30)]
DIFFERENCE_BINNING = [0.01 * i for i in range(-20, 30)]
SHOWER_WIDTH_BINNING = [0.1 * i for i in range(0, 21)]

#Jeremy comparison plots
ZCOORDINATE_BINNING = [10*i for i in range(50, 90)]
UPSTREAM_BINNING = [10*i for i in range(0, 50)]
UPSTREAM_ENERGY_BINS_VS = [5*i for i in range(2, 25)]
INLINE_BINNING = [1*i for i in range(0,26)]
INLINE_LOW_BINNING = [0.01*i for i in range(0,11)]
INLINE_HIGH_BINNING = [0.1*i for i in range(1,51)]
INLINE_W_BINNING = [0.1 * i for i in range(0,26)]
INLINE_W_BINNING_VS_VTX = [0.07*i for i in range(0, 25)]
VERTEXR_BINNING = [0.5*i for i in range(0,21)]
RADIAL_BINNING = [0.5*i for i in range(-10,11)]
#TRANSVERSE_SHOWER_BINNING = [0.05*i for i in range(0,41)]
TRANSVERSE_SHOWER_BINNING = [0.1*i for i in range(0,21)]
#EXCESS_ENERGY_FIT = [1.5, 3, 4.5, 6, 9, 12, 15, 20]
EETHETA_BINNING = [0.2*i for i in range(0,7)]
EXCESS_ENERGY_FIT = [1.5, 2, 3.5, 5, 7, 10]
LEAKAGE_BINNING  = [0.01 * i for i in range(-20,20)]

EOUTSIDE_VS_EE_BINS = [1.5, 2.5, 3.5, 4.5, 6, 8, 12, 16, 22]
#EOUTSIDE_VS_EE_BINS = [0.5, 1.5, 2.5, 3.5, 4.5, 6, 7, 12, 18]
PSI_EE_VS_EE_BINS = [1.5, 2.5, 3.5, 4.5, 6, 8, 12, 16, 22]
EE_VS_BINS = [2*i for i in range(0, 12)]
UPSTREAM_INLINE_ENERGY_BINS = [0.4, 1, 2, 4, 6.5, 10, 20, 40, 65, 100,150,250]  # designed to look ok on log10 axis



def FindQ2Bin(q2_val):
	q2_bin = None
	if q2_val >= QSQUARED_BINNING[-1]:
		return len(QSQUARED_BINNING)-1
	
	if q2_val < QSQUARED_BINNING[0]:
		print ("underflow")
		return 1
	
	for bin_num in range(1, len(QSQUARED_BINNING)):
		if q2_val >= QSQUARED_BINNING[bin_num-1] and q2_val < QSQUARED_BINNING[bin_num]:
			q2_bin = bin_num
			break
	assert q2_bin is not None
	return q2_bin

# the order and color of each band in stacked plots defined here.
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
# Categories["ExcessModel"] = {
#     "title":"Ad hoc excess model",
#     "color" : COLORS[2]
# }
#Categories["CCNuE"]={
#    "title" : "CC #nu_{e}",
#    "cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE","NonPhaseSpace"},
#    "color": COLORS[0]
#}
#Categories["CCNuMu"]= {
#    "title": "CC #nu_{#mu}",
#    "cate" : {"CCPi0","CCOther"},
#    "color" : COLORS[3]
#}
#Categories["NCCOH"]= {
#    "title": "NC Coh",
#    "color" : COLORS[4]
#}
#Categories["NC"]= {
#    "title": "Other NC",
#    "cate" : {"NCPi0","NCOther"},
#    "color" : COLORS[5]
#}

#Categories["NuEElastic"] = {
#    "title":"#nu + e elastic",
#    "color" : COLORS[6]
#}
#Categories["Other"] = {
#    "title":"Others",
#    "cate" : {"Other","NonFiducial","CCNuEAntiNu"},
#    "color" : COLORS[7]
#}


Categories["CCQElike"]={
    "title" : "CCQElike #nu_{e}",
    #"cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE","NonPhaseSpace"},
    "color": COLORS[0]
}

Categories["notCCQElike"]= {
    "title": "notCCQElike #nu_{e}",
    #"cate" : {"CCPi0","CCOther"},
    "color" : COLORS[1]
}

Categories["NCDiff"] = {
     "title":"NC Diffractive #pi^{0}",
     "color" : COLORS[2]
}

Categories["CCPi0"]= {
    "title": "Not #nu_{e} CC #pi^{0}", 
    "color" : COLORS[3]
}

Categories["NCCohPi0"]= {
    "title": "NC Coh #pi^{0}",
    "color" : COLORS[4]
}

Categories["NCPi0"]= {
    "title": "NC #pi^{0}",
    "color" : COLORS[5]
}

Categories["NuEElastic"] = {
    "title":"#nu + e elastic",
    "color" : COLORS[6]
}

Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other","NonFiducial","NonPhaseSpace"},
    "color" : COLORS[7]
}


#SignalDecomposition = {
#    "CCNuEQE" :
#    {
#        "title" : "CC #nu_{e}-QE",
#        "color": COLORS[0]
#    },
#    "CCNuEDelta" : {
#        "title" : "CC #nu_{e}-Delta",
#        "color": COLORS[1]
#    },
#    "CCNuEDIS" : {
#        "title" : "CC #nu_{e}-DIS",
#        "title" : "CC #nu_{e}-DIS",
#        "color": COLORS[2]
#    },
#    "CCNuE2p2h" : {
#        "title" : "CC #nu_{e}-2p2h",
#        "color": COLORS[3]
#    },
#    "Background" : {
#        "title" : "Backgrounds",
#        "cate": {"CCPi0","CCOther","NCCOH","NCPi0","NCOther","NuEElastic","Other","CCNuEAntiNu"},
#        "color": COLORS[4]
#    }
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
#}

#make sure all categories have a color and a name
#assert set(INT_COLORS.keys())==set(INT_NAMES.keys())

HISTS_TO_MAKE = [
    "Neutrino Energy",
    #"Visible Energy",
    "Energy Outside Cone",
    #"Psi*Ee",
    "Inline Upstream Energy Weighted Position",
    #"tDef2",
    "tDef1",
    "tDef1 Diffractive",
    "PiZeroE",
    "PiZero E Diffractive",
    #"Pion Transverse Momentum",
    #"Pion Transverse Momentum vs Electron Energy",
    "Inline Upstream Energy vs t",
    #"Lepton Transverse Momentum",
    #"Proton Transverse Momentum",
    #"Q0",
    #"Q3",
    #{"variables":["Lepton Pt"],
    #"tags": {"sideband","truth_class"}
    #},
    #"Visible Energy vs q3",
    #"Visible Energy vs q3 Migration",
    #"True Signal Visible Energy vs q3",
    #"True Signal q0 vs q3",
    #"Q0 Migration",
    #"Q3 Migration",
    #"Visible Energy Migration",
    #"Psi",
    #"Front dEdX",
    #"FrontdEdX_position",
    "Lepton Energy",
    #"Electron Theta",
    "Inline Upstream Energy",
    #"Psi vs Electron Energy",
    #"Inline Upstream Energy Weighted Position",
    #{"variables":["Visible Energy","Electron Energy Course Binning"],
    #"tags": {"sideband","truth_class"}
    #},
    #{"variables":["Visible Energy","Lepton Pt"],
    #"tags": {"sideband","truth_class"}
    #},
    #{"variables":["Print Arachne"],
    # "tags": {"sideband"},
    # "cuts": [lambda event: event.mc_intType == 4, lambda event: event.UpstreamInlineEnergy > 50]
    #},
    #"Long vs Trans",
    #"Extra Energy Long",
    #"Extra Energy Trans",
    #"Extra Energy Long Negative",
    #"Inline Upstream Energy Weighted Position Diffractive",
    #"Inline Upstream Energy Weighted Position NCPiZero",
    #"Vertex Difference vs Inline Upstream Energy Diffractive",
    #"Vertex Difference vs Inline Upstream Energy Coherent",
    #"Vertex Difference vs Inline Upstream Energy NCPiZero",
    #"Psi*Ee vs Ee",
    #"Energy Outside Cone vs Ee",
    #"Inline Upstream Arachne",
    #"Inline Upstream Energy Weighted Position Sarah's Events",
    #"Inline Weighted vs Vertex Difference",
    #"True Theta vs t",
    #"Reco Electron Theta",
    "t vs Electron Energy",
    "Upstream Energy",
    "Inline Upstream vs Electron Energy",
    "Weighted Inline Upstream vs Electron Energy",
    #"Energy Outside vs Electron Energy",
    "Electron Energy vs Pion Energy",
    #"Total Upstream vs Inline Upstream",
    "Inline Upstream vs t",
    "Ee and Eout vs Pion Energy",
    "t vs Upstream Energy",
    "Median Plane Shower Width",
    "Median Shower Width vs Ee",
    "Shower Width strip",
    "t vs PiZero E",
    "True Signal Pi Zero Energy",
    "True Signal t",
]
