import random
import math


def Gillespie_step(structure, C, k=None):
    ''' This functions takes concentrations of the enzymes and the structure info of a glycogen granules and
    return what is the next reaction to occurs and which time has been spent. (Following a gillespie algorithm)
    '''
    # propensity assuming mass action kinetics
    if k:
        h_gs = k["k_GS"]*C["GS"]*len(structure.Find_chain_for_gs())
        h_gp = k["k_GP"]*C["GP"]*len(structure.Find_chain_for_gp())
        h_gbe = k["k_GBE"]*C["GBE"]*len(structure.Find_chain_for_gbe())
        h_gde = k["k_GDE"]*C["GDE"]*len(structure.Find_chain_for_gde())
    else:
        h_gs = C["GS"]*len(structure.Find_chain_for_gs())
        h_gp = C["GP"]*len(structure.Find_chain_for_gp())
        h_gbe = C["GBE"]*len(structure.Find_chain_for_gbe())
        h_gde = C["GDE"]*len(structure.Find_chain_for_gde())

    a = h_gs + h_gp + h_gbe + h_gde

    if a == 0:
        return "no reaction can be proceed, all propensities are zero", 0
    r2 = random.uniform(0, a)
    r1 = random.uniform(0, 1)

    d_t = (1/a)*math.log(1/r1)
    if r2 < h_gs:
        return "Act_gs()", d_t
    if r2 >= h_gs and r2 < h_gs + h_gp:
        return "Act_gp()", d_t
    if r2 >= h_gs + h_gp and r2 < h_gs + h_gp + h_gbe:
        return "Act_gbe()", d_t
    if r2 >= h_gs + h_gp + h_gbe and r2 < h_gs + h_gp + h_gbe + +h_gde:
        return "Act_gde()", d_t
