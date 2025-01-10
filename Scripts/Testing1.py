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

NumberOfMatrices = 3
MatrixDimension  = 2
MatrixGenerationMethod = "Random"
ProbabliltiesGenerationMethod = "Equal"

Conditions = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabliltiesGenerationMethod)

print(Conditions)

DensityMatrices = Conditions.getDesityMatrices()

for iMatix in range(len(DensityMatrices)):
    print(f"Matrix number {iMatix}:\n",DensityMatrices[iMatix])

Solution = QSExSetUp.SDPSolver.SolveSDP(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("The found POVMs are:\n")

# POVMsum = picos.Constant([[0 for iDimension in range(MatrixDimension)] for jDimension in range(MatrixDimension)])

for iPOVM in range(len(Solution['POVMs'])):
    # POVMsum +=Solution['POVMs'][iPOVM]
    print(f"POVM_{iPOVM}:\n",Solution['POVMs'][iPOVM])

# print("POVMs Sum:\n",POVMsum)
