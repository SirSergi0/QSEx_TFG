#######################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################


import picos
import QSExSetUp
import numpy as np
from tqdm import tqdm

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [0.1 for i in range(NumberOfMatrices - 1)]
EigenValues.append(100)
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
GroupGeneratedConditions      = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

print("Initial\n",GroupGeneratedConditions.getGramMatrix())
GroupGeneratedConditions.variationUnitaryGramMatrix(1)
print("Initial\n",GroupGeneratedConditions.getGramMatrix())
# for Iteration in range (10):
#     print(f"Iteration{Iteration}:\n",variationGramMatrix(GroupGeneratedConditions,0.1*Iteration))
#     print("Trace", picos.trace(variationGramMatrix(GroupGeneratedConditions,0.1*Iteration)))



