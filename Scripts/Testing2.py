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

MyDensityMatrices      = Conditions.getDensityMatrices()
NumberOfMatrices       = Conditions.getNumberOfMatrices()
NumberOfDimensions     = Conditions.getMatrixDimension()
GramMatrix             = Conditions.getGramMatrix()


print(GramMatrix)
