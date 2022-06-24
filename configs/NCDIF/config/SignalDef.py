"""
   Define your signal and background here.

   author: H. Su
"""
from collections import OrderedDict
from config.CutConfig import KINEMATICS_CUTS
from tools.CutLibrary import CUTS
from tools import TruthTools

#helper functions for signal defination

def IsFiducial(event):
    try:
        return event.truth_is_fiducial == 1
    except RuntimeError:
        # this ntuple uses old variable name
        return event.truth_IsFiducial == 1

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2

Is2p2h = lambda event: event.mc_intType == 8
IsElastic = lambda event: event.mc_intType == 7
IsCoherent = lambda event: event.mc_intType == 4
IsPC = lambda event: event.mc_processType ==5
IsUnknown  = lambda event : event.mc_intType == 10
IsDiffractive = lambda event: event.mc_intType == 10
IsQE = lambda event: event.mc_intType == 1
IsRes = lambda event: event.mc_intType == 2


IsNuE = lambda event: abs(event.mc_incoming) == 12
IsAntiNu = lambda event: event.mc_incoming < 0 # assuming only neutrino incomes
IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG

IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsDeexcitationPhoton =  lambda event: event.mc_FSPartPDG[0] == 22 and event.mc_FSPartE[0] < 10
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

def IsInKinematicPhaseSpace(event):
    passed = True
    #electron_candidate_index = TruthTools.MostEnergeticParticle(event,11)
    #if electron_candidate_index is None:
    #    return False
    for cut in KINEMATICS_CUTS:
        passed = passed and CUTS["True{}".format(cut)].DoesEventPass(event)
    return passed


# In case a event satisfy multiple definations, the first takes priority.
TRUTH_CATEGORIES = OrderedDict()

TRUTH_CATEGORIES["NCDiff"] = lambda event: IsDiffractive(event)
TRUTH_CATEGORIES["CCQElike"] = lambda event: IsCC(event) and IsNuE(event) and not IsPi0InFinalState(event) and not IsMeson(event) and not IsHeavyBaryon(event) and not IsPhoton(event)
TRUTH_CATEGORIES["notCCQElike"] = lambda event: IsNuE(event) and (IsPi0InFinalState(event) or IsMeson(event) or IsPhoton(event) or IsHeavyBaryon(event))

TRUTH_CATEGORIES["NuEElastic"] = lambda event: IsElastic(event)
TRUTH_CATEGORIES["NonPhaseSpace"] = lambda event: IsCC(event) and IsNuE(event) and not IsInKinematicPhaseSpace(event)
TRUTH_CATEGORIES["NonFiducial"] = lambda event: IsCC(event) and IsNuE(event) and not IsFiducial(event)

TRUTH_CATEGORIES["CCPi0"] = lambda event: IsCC(event) and IsPi0InFinalState(event)
TRUTH_CATEGORIES["NCCohPi0"] = lambda event: IsCoherent(event) and IsNC(event) and IsPi0InFinalState(event)
TRUTH_CATEGORIES["NCPi0"] = lambda event: IsNC(event) and IsPi0InFinalState(event)
TRUTH_CATEGORIES["Other"] = lambda event: (not Is2p2h(event) and not IsQE(event) and not IsRes(event)) or not IsFiducial or (IsCC(event) and IsNuE(event) and not IsInKinematicPhaseSpace(event))


#TRUTH_CATEGORIES["ExcessModel"] = lambda event: IsPC(event) or IsUnknown(event)
#TRUTH_CATEGORIES["CCNuEQELike"] = lambda event: IsCC(event) and IsNuE(event) and not IsPi0InFinalState(event) and not IsMeson(event) and not IsHeavyBaryon(event) and not IsPhoton(event)
#TRUTH_CATEGORIES["CCNuE"] = lambda event: IsCC(event) and IsNuE(event)
#TRUTH_CATEGORIES["CCOther"] = lambda event: IsCC(event)
#TRUTH_CATEGORIES["NCCOH"] = lambda event: IsNC(event) and IsCoherent(event)
#TRUTH_CATEGORIES["NCOther"] = lambda event: IsNC(event)




# My signal is one or more of the listed categories.
SIGNAL_DEFINATION = [
    #"CCNuE",
    "NCDiff",
    #"CCQElike",
    #"notCCQElike"
    #"CCNuEQELike",
    #"NonPhaseSpace",
]

#Dump some categories to other to make plots easier to read:

EXTRA_OTHER = [
#    "NonFiducial",
#    "NonPhaseSpace",
]
