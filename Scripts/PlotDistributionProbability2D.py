########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     23/05/2025                                                                                #
#  Purpose:  Plot                                                                                      #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import mpltern
import picos
from collections import Counter
from picos.modeling.problem import SolutionFailure  


NumberOfMatrices              = 11
MatrixDimension               = 11
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
IterationsRandomSets          = 100000
DataMinimumError              = []
EigenValue0                   = []
EigenValue1                   = []
SDPMinimumErrorList           = []
SDPZeroErrorList              = []
LowerBound                    = []
LowerBoundZE                  = []
EigenValues = [1 for i in range (NumberOfMatrices-1)]
EigenValues.append(210)
EigenValuesNormalized   = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]

mode = "ZE"

for iIteration in tqdm(range(IterationsRandomSets),f"Computing points for the Minimum Error plot"):
    Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    if (iIteration == 0):
        print("Minimum Error Lower Bound: ",round(Conditions.getPerfectExlusionLowerBoundMinimumError(),4))
        print("Minimum Error random gess: ",(NumberOfMatrices-1)/NumberOfMatrices)
        print("Zero Error Lower Bound:    ",round(Conditions.getPerfectExlusionLowerBoundZeroError(),4))

    if Conditions.getPerfectExlusion():
        continue
    if mode == "ME" or mode == "both":
        MinimumErrorResult = round(QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)['SDPSolution'],4)
        if (MinimumErrorResult < round(Conditions.getPerfectExlusionLowerBoundMinimumError(),4)):
            raise ValueError ("SDP Minimum Error result", MinimumErrorResult,"Lower Bound Minimum Error", Conditions.getPerfectExlusionLowerBoundMinimumError(),"Conditions",Conditions)
        SDPMinimumErrorList.append(MinimumErrorResult)

    if mode == "ZE" or mode == "both":
        try: 
            ZeroErrorResult = round(QSExSetUp.SDPSolver.SolveSDPExclusionZeroError(Conditions)['SDPSolution'],4)
            if (ZeroErrorResult < round(Conditions.getPerfectExlusionLowerBoundZeroError(),4)):
                raise ValueError ("SDP Zero Error result", ZeroErrorResult,"Lower Bound Zero Error", Conditions.getPerfectExlusionLowerBoundMinimumError(),"Conditions",Conditions)
            SDPZeroErrorList.append(ZeroErrorResult)
        except SolutionFailure as e:
            print(f"SDP solver failed: {e}")
            continue

if mode == "ME" or mode == "both":
    plt.figure(figsize=(10, 6)) 
    plt.hist(SDPMinimumErrorList, bins=100, color = "orange",  label  = f"Number of entries: {len(SDPMinimumErrorList)}") 
    plt.axvline(x = Conditions.getPerfectExlusionLowerBoundMinimumError(), color='blue', linestyle='--', linewidth=2, label="LowerBound")
    plt.xlabel('Success probability minimum error')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f"../Plots/ExclusionMinimumErrorRandomDistributionZ{NumberOfMatrices}Prob{round(Conditions.getPerfectExlusionLowerBoundMinimumError(),3)}.pdf")
    plt.close()

if mode == "ZE" or mode == "both":
    plt.figure(figsize=(10, 6)) 
    plt.hist(SDPZeroErrorList, bins=100, color = "orange",  label  = f"Number of entries: {len(SDPZeroErrorList)}") 
    plt.axvline(x = Conditions.getPerfectExlusionLowerBoundZeroError(),color='blue', linestyle='--', linewidth=2, label="LowerBound")
    plt.xlabel('Success probability minimum error')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(f"../Plots/ExclusionZeroErrorRandomDistributionZ{NumberOfMatrices}Prob{round(Conditions.getPerfectExlusionLowerBoundZeroError(),3)}.pdf")
    plt.close()

