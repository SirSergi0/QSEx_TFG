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

NumberOfMatrices = 3
MatrixDimension  = 2
MatrixGenerationMethod = "RandomPureStates"
ProbabliltiesGenerationMethod = "Equal"

Conditions = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabliltiesGenerationMethod)

print(Conditions)

DensityMatrices = Conditions.getDensityMatrices()

for iMatix in range(len(DensityMatrices)):
    print(f"Matrix number {iMatix}:\n",DensityMatrices[iMatix])
"""
print("------------------------\n\
Computing the primal SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDP(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("------------------------\n\
Computing the Dual SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPDual(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))
"""
print("----------------------\n\
Pretty Good Measurement\n\
-----------------------")

GramMatrix = Conditions.getGramMatrix()

GramEigenValues, GraEigenVectors = np.linalg.eigh(GramMatrix)

SquareRoot = picos.Constant(sqrtm(GramMatrix))

SquareRootDiagonal = np.diagonal(SquareRoot.value)

sum = 0
for iElement in SquareRootDiagonal:
    sum += abs(iElement)**2

print(SquareRoot)

print(SquareRootDiagonal)

print(sum)
