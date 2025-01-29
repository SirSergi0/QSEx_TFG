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

NumberOfMatrices              = 3
MatrixDimension               = 2
MatrixGenerationMethod        = "GroupGeneratedStates"
ProbabliltiesGenerationMethod = "Equal"
initialState                  = picos.Constant([1,0])
sigmaX                        = picos.Constant([[0,1],[1,0]])
Conditions                    = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabliltiesGenerationMethod, seedState = initialState, involutionalMatrix = sigmaX)
DensityMatrices               = Conditions.getDensityMatrices()

print(Conditions)

for iMatix in range(len(DensityMatrices)):
    print(f"Matrix number {iMatix}:\n",DensityMatrices[iMatix])


