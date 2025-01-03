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

NumberOfMatrices       = 10
MatrixDimension        = 3
RandomConditions       = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension)

print(RandomConditions)

Solution = QSExSetUp.SDPSolver.SolveSDP(RandomConditions)

print("The given problem has been", Solution['SDPSolution'])

print("The POVMs found are:\n")

for iPOVM in range(len(Solution['POVMs'])):
    print(f"POVM_{iPOVM}:\n",Solution['POVMs'][iPOVM])
