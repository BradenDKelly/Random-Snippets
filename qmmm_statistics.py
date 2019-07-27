# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:07:29 2019

@author: Zarathustra
"""

import pandas as pd
import numpy as np
import functools

import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def MakePlot(x1,y1,x2,y2, name, marker='o'):
    
    fig = plt.figure()

    #MP2_mbis_polar_error= np.sqrt(df.expError[conjunction(c_2,c_3,c_9,c_13, c_16)] + df.calcError[conjunction(c_2,c_3,c_9,c_13, c_16)])
    #expErrors = np.sqrt(df.expError[conjunction(c_2,c_3,c_9,c_13, c_16)])
    ax = fig.add_subplot(1,2,1,aspect='equal')
    ax.text(sp1x, sp1y, pbanner,style='italic',fontsize=pbsize,bbox={'facecolor': bfc, 'alpha': balpha,'edgecolor':bec, 'pad': pad, 'lw':lw})
    ax.scatter(x1,y1,marker=marker,edgecolors=mec,linewidths=mew,c=mfc,s=sm,alpha=alpha) #,s=expErrors )
    ax.plot( [minl, maxl],[minl, maxl],c='k',linewidth=linewidth )
    ax.set_ylabel(ylabel) 
    ax.set_xlabel(xlabel) 
    ax.set_xlim(minl, maxl)
    ax.set_ylim(minl, maxl)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(minorTicks))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(minorTicks))
    ax.tick_params(direction='in', which='both',bottom=True, top=True, left=True, right=True)
    
    
    
    #MP2_mbis_fixed_error= np.sqrt(df.expError[conjunction(c_2,c_3,c_9,c_13, c_16)] + df.calcError[conjunction(c_2,c_3,c_9,c_13, c_16)])
    ax = fig.add_subplot(1,2,2,aspect='equal')
    ax.text(sp2x, sp2y, fbanner,style='italic',fontsize=pbsize,bbox={'facecolor': bfc, 'alpha': balpha,'edgecolor':bec, 'pad': pad, 'lw':lw})
    ax.scatter(x2,y2,marker='o',edgecolors=mec,linewidths=mew,c=mfc,s=sm,alpha=alpha)
    ax.plot( [minl, maxl],[minl, maxl],c='k',linewidth=linewidth )
    ax.set_xlabel(xlabel) 
    ax.set_xlim(minl, maxl)
    ax.set_ylim(minl, maxl)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(minorTicks))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(minorTicks))
    ax.tick_params(direction='in', which='both',bottom=True, top=True, left=True, right=True)
    
    fig.tight_layout()
    
    plt.show()
    fig.set_size_inches(width, height)
    fig.savefig(name,dpi=dpi)

###################################
#   Read in Excel file into Pandas
###################################    
    
data = pd.read_excel (r'Hydration_Energies_For_Plots.xlsx')

###############################
#   Cleanup up all empty cells
###############################
data.replace('', np.nan, inplace=True)
data.dropna( inplace=True)

df = pd.DataFrame(data)

pd.to_numeric(df.calc)
pd.to_numeric(df.experiment)
pd.to_numeric(df.calcError)
pd.to_numeric(df.expError)
pd.to_numeric(df.FreeSolv)
pd.to_numeric(df.Riquelme)

"""This function is for filtering dataframe results"""
def conjunction(*conditions):
    return functools.reduce(np.logical_and, conditions)

"""These are the many conditions that can be used for filtering"""
c_1 = df.qqType == "fixed"
c_2 = df.qqType == "polar"
c_3 = df.stepsPerUpdate = 10000
c_4 = df.nhAtoms >= 6
c_4b = df.nhAtoms < 6
c_5 = df.structure == "rotatable"
c_6 = df.structure == "cyclic"
c_7 = df.N == True
c_8 = df.O == True
c_9 = df.chargeMethod == "mbis"
c_10 = df.chargeMethod == "resp"
c_11 = df.chargeMethod == "bcc"
c_12 = df.qmTheory == "B3LYP"
c_13 = df.qmTheory == "MP2"
c_14 = df.qmTheory == "AM1"
c_15 = df.qmTheory == "HF"
c_16 = df.qmBasis == "cc-pVTZ"
c_17 = df.qmBasis == "6-31G*"
c_18 = abs(df.FreeSolv) > 0.0001
c_19 = abs(df.Riquelme) > 0.0001


B3LYP_mbis_polar_rmse   = ((df.experiment[conjunction(c_2,c_3,c_9,c_12, c_16)]            - df.calc[conjunction(c_2,c_3,c_9,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_mbis_fixed_rmse   = ((df.experiment[conjunction(c_1,c_3,c_9,c_12, c_16)]            - df.calc[conjunction(c_1,c_3,c_9,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_mbis_polar_small_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4b,c_9,c_12, c_16)] - df.calc[conjunction(c_2,c_3,c_4b,c_9,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_mbis_polar_big_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4,c_9,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_4,c_9,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_mbis_polar_rot_rmse   = ((df.experiment[conjunction(c_2,c_3,c_5,c_9,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_5,c_9,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_mbis_polar_cyc_rmse   = ((df.experiment[conjunction(c_2,c_3,c_6,c_9,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_6,c_9,c_12, c_16)]) ** 2).mean() ** .5

B3LYP_resp_polar_rmse   = ((df.experiment[conjunction(c_2,c_3,c_10,c_12, c_16)]            - df.calc[conjunction(c_2,c_3,c_10,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_resp_fixed_rmse   = ((df.experiment[conjunction(c_1,c_3,c_10,c_12, c_16)]            - df.calc[conjunction(c_1,c_3,c_10,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_resp_polar_small_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4b,c_10,c_12, c_16)] - df.calc[conjunction(c_2,c_3,c_4b,c_10,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_resp_polar_big_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4,c_10,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_4,c_10,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_resp_polar_rot_rmse   = ((df.experiment[conjunction(c_2,c_3,c_5,c_10,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_5,c_10,c_12, c_16)]) ** 2).mean() ** .5
B3LYP_resp_polar_cyc_rmse   = ((df.experiment[conjunction(c_2,c_3,c_6,c_10,c_12, c_16)]    - df.calc[conjunction(c_2,c_3,c_6,c_10,c_12, c_16)]) ** 2).mean() ** .5

HF_resp_polar_rmse       = ((df.experiment[conjunction(c_2,c_3,c_10,c_12, c_16)]      - df.calc[conjunction(c_2,c_3,c_10,c_12, c_16)]) ** 2).mean() ** .5
HF_resp_polar_small_rmse = ((df.experiment[conjunction(c_2,c_3,c_4b,c_10,c_12, c_16)] - df.calc[conjunction(c_2,c_3,c_4b,c_10,c_12, c_16)]) ** 2).mean() ** .5
HF_resp_polar_big_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4,c_10,c_12, c_16)]  - df.calc[conjunction(c_2,c_3,c_4,c_10,c_12, c_16)]) ** 2).mean() ** .5
HF_resp_polar_rot_rmse   = ((df.experiment[conjunction(c_2,c_3,c_5,c_10,c_12, c_16)]  - df.calc[conjunction(c_2,c_3,c_5,c_10,c_12, c_16)]) ** 2).mean() ** .5
HF_resp_polar_cyc_rmse   = ((df.experiment[conjunction(c_2,c_3,c_6,c_10,c_12, c_16)]  - df.calc[conjunction(c_2,c_3,c_6,c_10,c_12, c_16)]) ** 2).mean() ** .5

HF_resp_fixed_rmse       = ((df.experiment[conjunction(c_1,c_3,c_10,c_15, c_17)]      - df.calc[conjunction(c_1,c_3,c_10,c_15, c_17)]) ** 2).mean() ** .5
HF_resp_fixed_small_rmse = ((df.experiment[conjunction(c_1,c_3,c_4b,c_10,c_15, c_17)] - df.calc[conjunction(c_1,c_3,c_4b,c_10,c_15, c_17)]) ** 2).mean() ** .5
HF_resp_fixed_big_rmse   = ((df.experiment[conjunction(c_1,c_3,c_4,c_10,c_15, c_17)]  - df.calc[conjunction(c_1,c_3,c_4,c_10,c_15, c_17)]) ** 2).mean() ** .5
HF_resp_fixed_rot_rmse   = ((df.experiment[conjunction(c_1,c_3,c_5,c_10,c_15, c_17)]  - df.calc[conjunction(c_1,c_3,c_5,c_10,c_15, c_17)]) ** 2).mean() ** .5
HF_resp_fixed_cyc_rmse   = ((df.experiment[conjunction(c_1,c_3,c_6,c_10,c_15, c_17)]  - df.calc[conjunction(c_1,c_3,c_6,c_10,c_15, c_17)]) ** 2).mean() ** .5
                              
MP2_mbis_polar_rmse     = ((df.experiment[conjunction(c_2,c_3,c_9,c_13, c_16)]          - df.calc[conjunction(c_2,c_3,c_9,c_13, c_16)]) ** 2).mean() ** .5
MP2_mbis_fixed_rmse     = ((df.experiment[conjunction(c_1,c_3,c_9,c_13, c_16)]          - df.calc[conjunction(c_1,c_3,c_9,c_13, c_16)]) ** 2).mean() ** .5
MP2_mbis_polar_small_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4b,c_9,c_13, c_16)] - df.calc[conjunction(c_2,c_3,c_4b,c_9,c_13, c_16)]) ** 2).mean() ** .5
MP2_mbis_polar_big_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4,c_9,c_13, c_16)]    - df.calc[conjunction(c_2,c_3,c_4,c_9,c_13, c_16)]) ** 2).mean() ** .5
MP2_mbis_polar_rot_rmse   = ((df.experiment[conjunction(c_2,c_3,c_5,c_9,c_13, c_16)]    - df.calc[conjunction(c_2,c_3,c_5,c_9,c_13, c_16)]) ** 2).mean() ** .5
MP2_mbis_polar_cyc_rmse   = ((df.experiment[conjunction(c_2,c_3,c_6,c_9,c_13, c_16)]    - df.calc[conjunction(c_2,c_3,c_6,c_9,c_13, c_16)]) ** 2).mean() ** .5

BCC_polar_rmse          = ((df.experiment[conjunction(c_2,c_3,c_11)]      - df.calc[conjunction(c_2,c_3,c_11)]) ** 2).mean() ** .5
BCC_polar_small_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4b,c_11)]  - df.calc[conjunction(c_2,c_3,c_4b,c_11)]) ** 2).mean() ** .5
BCC_polar_big_rmse   = ((df.experiment[conjunction(c_2,c_3,c_4,c_11)]     - df.calc[conjunction(c_2,c_3,c_4,c_11)]) ** 2).mean() ** .5
BCC_polar_rot_rmse   = ((df.experiment[conjunction(c_2,c_3,c_5,c_11)]     - df.calc[conjunction(c_2,c_3,c_5,c_11)]) ** 2).mean() ** .5
BCC_polar_cyc_rmse   = ((df.experiment[conjunction(c_2,c_3,c_6,c_11)]     - df.calc[conjunction(c_2,c_3,c_6,c_11)]) ** 2).mean() ** .5

BCC_fixed_rmse          = ((df.experiment[conjunction(c_1,c_3,c_11)]      - df.calc[conjunction(c_1,c_3,c_11)]) ** 2).mean() ** .5
BCC_fixed_small_rmse   = ((df.experiment[conjunction(c_1,c_3,c_4b,c_11)]  - df.calc[conjunction(c_1,c_3,c_4b,c_11)]) ** 2).mean() ** .5
BCC_fixed_big_rmse   = ((df.experiment[conjunction(c_1,c_3,c_4,c_11)]     - df.calc[conjunction(c_1,c_3,c_4,c_11)]) ** 2).mean() ** .5
BCC_fixed_rot_rmse   = ((df.experiment[conjunction(c_1,c_3,c_5,c_11)]     - df.calc[conjunction(c_1,c_3,c_5,c_11)]) ** 2).mean() ** .5
BCC_fixed_cyc_rmse   = ((df.experiment[conjunction(c_1,c_3,c_6,c_11)]     - df.calc[conjunction(c_1,c_3,c_6,c_11)]) ** 2).mean() ** .5

B3LYP_resp_polar_rmse   = ((df.experiment[conjunction(c_2,c_3,c_10,c_12)] - df.calc[conjunction(c_2,c_3,c_10,c_12)]) ** 2).mean() ** .5

FreeSolv_rmse           = ((df.experiment[conjunction(c_18)]    - df.FreeSolv[conjunction(c_18)]) ** 2).mean() ** .5
FreeSolv_small_rmse   = ((df.experiment[conjunction(c_18,c_4b)] - df.FreeSolv[conjunction(c_18,c_4b)]) ** 2).mean() ** .5
FreeSolv_big_rmse   = ((df.experiment[conjunction(c_18,c_4)]    - df.FreeSolv[conjunction(c_18,c_4)]) ** 2).mean() ** .5
FreeSolv_rot_rmse   = ((df.experiment[conjunction(c_18,c_5)]    - df.FreeSolv[conjunction(c_18,c_5)]) ** 2).mean() ** .5
FreeSolv_cyc_rmse   = ((df.experiment[conjunction(c_18,c_6)]    - df.FreeSolv[conjunction(c_18,c_6)]) ** 2).mean() ** .5

Riquelme_rmse           = ((df.experiment[conjunction(c_19)]     - df.Riquelme[conjunction(c_19)]) ** 2).mean() ** .5
Riquelme_small_rmse   = ((df.experiment[conjunction(c_19,c_4b)]  - df.Riquelme[conjunction(c_19,c_4b)]) ** 2).mean() ** .5
Riquelme_big_rmse   = ((df.experiment[conjunction(c_19,c_4)]     - df.Riquelme[conjunction(c_19,c_4)]) ** 2).mean() ** .5
Riquelme_rot_rmse   = ((df.experiment[conjunction(c_19,c_5)]     - df.Riquelme[conjunction(c_19,c_5)]) ** 2).mean() ** .5
Riquelme_cyc_rmse   = ((df.experiment[conjunction(c_19,c_6)]     - df.Riquelme[conjunction(c_19,c_6)]) ** 2).mean() ** .5                          
                          
                          
                          
                          
print("===================================================================================")
print("B3LYP/cc-pVTZ mbis polar :", B3LYP_mbis_polar_rmse, " kJ/mol ", B3LYP_mbis_polar_rmse / 4.184, "  kcal/mol"  )
print("MP2/cc-pVTZ mbis polar   :", MP2_mbis_polar_rmse, "  kJ/mol ", MP2_mbis_polar_rmse / 4.184, " kcal/mol"  )
print("B3LYP/cc-pVTZ resp polar :", B3LYP_resp_polar_rmse, " kJ/mol ", B3LYP_resp_polar_rmse / 4.184, "  kcal/mol"  )
print("HF/6-31G* resp polar     :", HF_resp_polar_rmse, " kJ/mol ", HF_resp_polar_rmse / 4.184, "  kcal/mol"  )
print("HF/6-31G* resp fixed     :", HF_resp_fixed_rmse, " kJ/mol ", HF_resp_fixed_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC fixed             :", BCC_fixed_rmse, "  kJ/mol ", BCC_fixed_rmse / 4.184, " kcal/mol"  )
print("FreeSolv database        :", FreeSolv_rmse, " kJ/mol ", FreeSolv_rmse / 4.184, " kcal/mol" )
print("Riquelme results         :", Riquelme_rmse, " kJ/mol ", Riquelme_rmse / 4.184, "  kcal/mol" )
print("B3LYP/cc-pVTZ mbis fixed :", B3LYP_mbis_fixed_rmse, "  kJ/mol ", B3LYP_mbis_fixed_rmse / 4.184, "  kcal/mol"  )
print("MP2/cc-pVTZ mbis fixed   :", MP2_mbis_fixed_rmse, " kJ/mol ", MP2_mbis_fixed_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("Best case Ratios: ")
print("B3LYP: ", (B3LYP_mbis_polar_rmse / 4.184) * 2.0 / (Riquelme_rmse / 4.184), " kcal/mol" )
print("MP2:   ", (MP2_mbis_polar_rmse / 4.184) * 2.0 / (Riquelme_rmse / 4.184), " kcal/mol" )
print("====================================")
print("counting mbis b3lyp: ", df.calc[conjunction(c_2,c_3,c_9,c_12, c_16)].count()   )
print("counting mbis MP2: ", df.calc[conjunction(c_2,c_3,c_9,c_13, c_16)].count()   )
print("counting FreeSolv: ", df.calc[conjunction(c_18)].count()   )
print("counting Riquelme: ", df.calc[conjunction(c_19)].count()   )
print("counting AM1BCC spce: ", df.calc[conjunction(c_1,c_11)].count()   )

print("===================================================================================")
print("B3LYP/cc-pVTZ mbis polar small :", B3LYP_mbis_polar_small_rmse, " kJ/mol ", B3LYP_mbis_polar_small_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ mbis big  :", B3LYP_mbis_polar_big_rmse, " kJ/mol ", B3LYP_mbis_polar_big_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ mbis rot  :", B3LYP_mbis_polar_rot_rmse, " kJ/mol ", B3LYP_mbis_polar_rot_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ mbis cyclic :", B3LYP_mbis_polar_cyc_rmse, " kJ/mol ", B3LYP_mbis_polar_cyc_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("MP2/cc-pVTZ mbis polar small :", MP2_mbis_polar_small_rmse, " kJ/mol ", MP2_mbis_polar_small_rmse / 4.184, "  kcal/mol"  )
print("MP2/cc-pVTZ mbis big  :", MP2_mbis_polar_big_rmse, " kJ/mol ", MP2_mbis_polar_big_rmse / 4.184, "  kcal/mol"  )
print("MP2/cc-pVTZ mbis rot  :", MP2_mbis_polar_rot_rmse, " kJ/mol ", MP2_mbis_polar_rot_rmse / 4.184, "  kcal/mol"  )
print("MP2/cc-pVTZ mbis cyclic :", MP2_mbis_polar_cyc_rmse, " kJ/mol ", MP2_mbis_polar_cyc_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("B3LYP/cc-pVTZ resp polar small  :", B3LYP_resp_polar_small_rmse, " kJ/mol ", B3LYP_resp_polar_small_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ resp polar big    :", B3LYP_resp_polar_big_rmse, " kJ/mol ", B3LYP_resp_polar_big_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ resp polar rot    :", B3LYP_resp_polar_rot_rmse, " kJ/mol ", B3LYP_resp_polar_rot_rmse / 4.184, "  kcal/mol"  )
print("B3LYP/cc-pVTZ resp polar cyclic :", B3LYP_resp_polar_cyc_rmse, " kJ/mol ", B3LYP_resp_polar_cyc_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("AM1BCC polar small  :", BCC_polar_small_rmse, " kJ/mol ", BCC_polar_small_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC polar big    :", BCC_polar_big_rmse, " kJ/mol ", BCC_polar_big_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC polar rot    :", BCC_polar_rot_rmse, " kJ/mol ", BCC_polar_rot_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC polar cyclic :", BCC_polar_cyc_rmse, " kJ/mol ", BCC_polar_cyc_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("AM1BCC fixed small  :", BCC_fixed_small_rmse, " kJ/mol ", BCC_fixed_small_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC fixed big    :", BCC_fixed_big_rmse, " kJ/mol ", BCC_fixed_big_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC fixed rot    :", BCC_fixed_rot_rmse, " kJ/mol ", BCC_fixed_rot_rmse / 4.184, "  kcal/mol"  )
print("AM1BCC fixed cyclic :", BCC_fixed_cyc_rmse, " kJ/mol ", BCC_fixed_cyc_rmse / 4.184, "  kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("HF/6-31G* resp polar small  :", HF_resp_polar_small_rmse, " kJ/mol ", HF_resp_polar_small_rmse / 4.184, " kcal/mol"  )
print("HF/6-31G* resp polar big    :", HF_resp_polar_big_rmse, " kJ/mol ", HF_resp_polar_big_rmse / 4.184, " kcal/mol"  )
print("HF/6-31G* resp polar rot    :", HF_resp_polar_rot_rmse, " kJ/mol ", HF_resp_polar_rot_rmse / 4.184, "  kcal/mol"  )
print("HF/6-31G* resp polar cyclic :", HF_resp_polar_cyc_rmse, " kJ/mol ", HF_resp_polar_cyc_rmse / 4.184, " kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("HF/6-31G* resp fixed small  :", HF_resp_fixed_small_rmse, " kJ/mol ", HF_resp_fixed_small_rmse / 4.184, " kcal/mol"  )
print("HF/6-31G* resp fixed big    :", HF_resp_fixed_big_rmse, " kJ/mol ", HF_resp_fixed_big_rmse / 4.184, "  kcal/mol"  )
print("HF/6-31G* resp fixed rot    :", HF_resp_fixed_rot_rmse, " kJ/mol ", HF_resp_fixed_rot_rmse / 4.184, "  kcal/mol"  )
print("HF/6-31G* resp fixed cyclic :", HF_resp_fixed_cyc_rmse, " kJ/mol ", HF_resp_fixed_cyc_rmse / 4.184, "kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("FreeSolv resp small  :", FreeSolv_small_rmse, " kJ/mol ", FreeSolv_small_rmse / 4.184, " kcal/mol"  )
print("FreeSolv resp big    :", FreeSolv_big_rmse, " kJ/mol ", FreeSolv_big_rmse / 4.184, "  kcal/mol"  )
print("FreeSolv resp rot    :", FreeSolv_rot_rmse, " kJ/mol ", FreeSolv_rot_rmse / 4.184, "   kcal/mol"  )
print("FreeSolv resp cyclic :", FreeSolv_cyc_rmse, " kJ/mol ", FreeSolv_cyc_rmse / 4.184, " kcal/mol"  )
print("===================================================================================")
print("===================================================================================")
print("Riquelme small  :", Riquelme_small_rmse, " kJ/mol ", Riquelme_small_rmse / 4.184, " kcal/mol"  )
print("Riquelme big    :", Riquelme_big_rmse, " kJ/mol ", Riquelme_big_rmse / 4.184, "  kcal/mol"  )
print("Riquelme rot    :", Riquelme_rot_rmse, " kJ/mol ", Riquelme_rot_rmse / 4.184, "  kcal/mol"  )
print("Riquelme cyclic :", Riquelme_cyc_rmse, "kJ/mol ", Riquelme_cyc_rmse / 4.184, "kcal/mol"  )
print("===================================================================================")

#pl.grid(which='both', color='w', lw=0.25, axis='y', zorder=12)
#pl.ylabel(r'$\mathrm{\langle{\frac{ \partial U } { \partial \lambda }}\rangle_{\lambda}\/%s}$' % P.units, fontsize=20, color='#151B54')
#pl.annotate('$\mathit{\lambda}$', xy=(0, 0), xytext=(0.5, -0.05), size=18, textcoords='axes fraction', va='top', ha='center', color='#151B54')

plt.rc('font', family='sans-serif', serif='Garamond')
#plt.rc('text', usetex=True)
#plt.rcParams["text.usetex"] = False
plt.rc('xtick', labelsize=4)
plt.rc('ytick', labelsize=4)
plt.rc('axes', labelsize=6)
plt.rcParams['xtick.major.size']   = 4
plt.rcParams['xtick.major.width']  = 0.4
plt.rcParams['xtick.minor.size']   = 2
plt.rcParams['xtick.minor.width']  = 0.3
plt.rcParams['ytick.major.size']   = 4
plt.rcParams['ytick.major.width']  = 0.4
plt.rcParams['ytick.minor.size']   = 2
plt.rcParams['ytick.minor.width']  = 0.3
plt.rcParams['axes.linewidth']     = 0.4
#plt.rcParams['scatter.edgecolors'] = 'm'

width =5 #3.487
height = width / 1.618
ylabel = r'Experiment $\mathrm{ \Delta G_{hyd} }$ (kJ/mol)'
xlabel = r'Calculated $\mathrm{ \Delta G_{hyd} }$ (kJ/mol)'
linewidth=0.5
markersize=6
title_font=5
maxl=25
minl=-110
minorTicks = 10
#######################
# marker options
#######################
sm=12       # marker size
mfc='m'    # face color
alpha=1.0  # transparency
mec='k'    # edge color
mew=0.2    # edge width
#########################
#  banner inside figures
#########################
pbanner='Fluctuating & Polarized'
fbanner='Fixed'
fbsize=6  # fixed banner font size
pbsize=6  # polarized banner font size
sp1x=-70  # polarized banner x position
sp1y=10   # both banners y position
sp2x=-47.5  # fixed banner x position
sp2y=sp1y
lw=0.4     # lineweight of outside box
pad=2
bec='k'
bfc='cyan'
balpha=0.9
dpi=800   # resolution of pdf image

#fig, ax = plt.subplots()
#fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

###############################################################################
###############################################################################
#                 MP2/cc-pVTZ Polar and Fixed
###############################################################################
############################################################################### 
x1=df.calc[conjunction(c_2,c_3,c_9,c_13, c_16)]          #MP2 polar
y1=df.experiment[conjunction(c_2,c_3,c_9,c_13, c_16)]
x2=df.calc[conjunction(c_1,c_3,c_9,c_13, c_16)]          #MP2 fixed
y2=df.experiment[conjunction(c_1,c_3,c_9,c_13, c_16)]
outputName="MP2_MBIS.pdf"
    
MakePlot(x1,y1,x2,y2,outputName)    

###############################################################################
###############################################################################
#                B3LYP/cc-pVTZ Polar and Fixed
###############################################################################
############################################################################### 
x1=df.calc[conjunction(c_2,c_3,c_9,c_12, c_16)]          
y1=df.experiment[conjunction(c_2,c_3,c_9,c_12, c_16)]
x2=df.calc[conjunction(c_1,c_3,c_9,c_12, c_16)]         
y2=df.experiment[conjunction(c_1,c_3,c_9,c_12, c_16)]
outputName="B3LYP_MBIS.pdf"
    
MakePlot(x1,y1,x2,y2,outputName) 

###############################################################################
###############################################################################
#               FreeSolv AM1BCC and Riquelme MBIS
###############################################################################
############################################################################### 
x1=df.FreeSolv[conjunction(c_18)]          
y1=df.experiment[conjunction(c_18)]
x2=df.Riquelme[conjunction(c_19)]        
y2=df.experiment[conjunction(c_19)]
outputName="FreeSolv_Riquelme.pdf"
    
MakePlot(x1,y1,x2,y2,outputName)

###############################################################################
###############################################################################
#                       AM1-BCC Polar and Fixed
###############################################################################
############################################################################### 
x1=df.calc[conjunction(c_2,c_3,c_11)]         
y1=df.experiment[conjunction(c_2,c_3,c_11)]
x2=df.calc[conjunction(c_1,c_3,c_11)]        
y2=df.experiment[conjunction(c_1,c_3,c_11)]
outputName="AM1-BCC"
    
MakePlot(x1,y1,x2,y2,outputName)

###############################################################################
###############################################################################
#               HF/6-31G* RESP Polar and Fixed
###############################################################################
############################################################################### 
x1=df.calc[conjunction(c_2,c_3,c_10,c_15, c_17)]         
y1=df.experiment[conjunction(c_2,c_3,c_10,c_15, c_17)]
x2=df.calc[conjunction(c_1,c_3,c_10,c_15, c_17)]        
y2=df.experiment[conjunction(c_1,c_3,c_10,c_15, c_17)]
outputName="HF_RESP"
    
MakePlot(x1,y1,x2,y2,outputName)

###############################################################################
###############################################################################
#               B3LYP/cc-pVTZ RESP Polar and Fixed
###############################################################################
############################################################################### 
x1=df.calc[conjunction(c_2,c_3,c_10,c_12, c_16)]        
y1=df.experiment[conjunction(c_2,c_3,c_10,c_12, c_16)]
x2=df.calc[conjunction(c_1,c_3,c_10,c_12, c_16)]       
y2=df.experiment[conjunction(c_1,c_3,c_10,c_12, c_16)]
outputName="B3LYP_RESP"
    
MakePlot(x1,y1,x2,y2,outputName)


###############################################################################
###############################################################################
