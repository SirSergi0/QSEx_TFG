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
import matplotlib.pyplot as plt
from tqdm import tqdm

NumberOfMatrices              = 11
MatrixDimension               = 11
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [0.1 for i in range(NumberOfMatrices - 1)]
EigenValues.append(100)
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
Iterations                    = 10
StaticIterations              = 50
Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
GroupGeneratedConditions      = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

epsilon = 0.001

DistanceFrobenius             = []
Distance1                     = []
DistanceInfinite              = []
DistanceSpectral              = []
DistanceNuclear               = []
ExlusionSuccessProbability    = []

for iIteration in tqdm(range(Iterations),"Computing SDP"):
    for jIteration in range(StaticIterations):
        Conditions  = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
        Conditions.variationUnitaryGramMatrix(epsilon*iIteration)
        Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)
        ExlusionSuccessProbability.append(round(Solution['SDPSolution'],10))
        DistanceFrobenius.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),2))
        Distance1.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),1))
        DistanceInfinite.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),float("inf")))
        DistanceSpectral.append(picos.SpectralNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))
        DistanceNuclear.append(picos.NuclearNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))

Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(GroupGeneratedConditions)
print(round(Solution['SDPSolution'],10))

plt.figure(figsize=(10, 7))
plt.scatter(DistanceFrobenius,ExlusionSuccessProbability, label = "Frobenius")
# plt.scatter(Distance1,ExlusionSuccessProbability, label = "D1")
# plt.scatter(DistanceInfinite,ExlusionSuccessProbability, label = "Dinf")
# plt.scatter(DistanceSpectral,ExlusionSuccessProbability, label = "Spectral")
# plt.scatter(DistanceNuclear,ExlusionSuccessProbability, label = "Nuclear")
plt.xlabel("Norm Distance")
plt.ylabel("Exclusion success probability")
plt.legend(loc='upper right')
plt.savefig(f"../Plots/ExclusionNormPlots{NumberOfMatrices}.pdf")
plt.show()
plt.close()
