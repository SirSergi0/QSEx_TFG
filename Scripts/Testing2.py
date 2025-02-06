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

NumberOfMatrices              = 4
MatrixDimension               = 4
MatrixGenerationMethod        = "Zn"
ProbabliltiesGenerationMethod = "Equal"
Overlap                       = 0.5
Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod, Overlap = Overlap)

print(Conditions)

# print("------------------------\n\
# Computing the primal SDP\n\
# ------------------------")
#
# Solution = QSExSetUp.SDPSolver.SolveSDPMinimumError(Conditions)
#
# print("The given problem has been:\n", Solution['SDPSolution'])
#
# print("The success probability is: ", round(Solution['SDPSolution'],4))

print("------------------------\n\
Computing the Dual SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("----------------------------\n\
Computing the Zero Error SDP\n\
----------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPZeroError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))


