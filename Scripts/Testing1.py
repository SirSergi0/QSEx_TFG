########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import picos
import numpy as np

NumberOfMatrices = 4 
MatrixDimension  = 3
GenerationMethod = 'Random'
Overlap          = 0

Conditions = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension,GenerationMethod, Overlap = Overlap)

print(Conditions)

DensityMatrices = Conditions.getDesityMatrices()

for iMatix in range(len(DensityMatrices)):
    print(f"Matrix number {iMatix}:\n",DensityMatrices[iMatix])

Solution = QSExSetUp.SDPSolver.SolveSDP(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("The found POVMs are:\n")

for iPOVM in range(len(Solution['POVMs'])):
    print(f"POVM_{iPOVM}:\n",Solution['POVMs'][iPOVM])
