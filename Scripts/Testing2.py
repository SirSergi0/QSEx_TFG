########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import picos
import numpy as np
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt

NumberOfMatrices              = 10
MatrixDimension               = 10
MatrixGenerationMethod        = "Zn"
ProbabliltiesGenerationMethod = "Equal"
Overlap                       = 0
Acuracy                       = 100
SuccesProbabilityList         = []
OverlapsList                  = []

for i in range (1,Acuracy):
    Overlap += 1/Acuracy
    Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod, Overlap = Overlap)
    SuccesProbabilityList.append(Conditions.getSRMSuccessProbability())
    OverlapsList.append(Overlap)

plt.scatter(OverlapsList, SuccesProbabilityList)
plt.xlabel("Overlap")
plt.ylabel("Success Probability")
plt.title(f"Group generated Z{NumberOfMatrices}")
plt.savefig(f"../Plots/OverlapVSSucessProbabilityZ{NumberOfMatrices}.pdf")


"""
print("------------------------\n\
Computing the primal SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPMinimumError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))
"""
