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

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [1 for i in range(NumberOfMatrices - 1)]
EigenValues.append(100)
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
Iterations                    = 10
StaticIterations              = 50
Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
GroupGeneratedConditions      = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

epsilon = 0.01

DistanceFrobenius             = []
Distance1                     = []
DistanceInfinite              = []
DistanceSpectral              = []
DistanceNuclear               = []
ExlusionSuccessProbability    = []

MinimumErrorLowerBound = (round(Conditions.getPerfectExlusionLowerBoundMinimumError(),3))

PossibleErrors = open(f"../Plots/PlotsExclusionNormPlots{NumberOfMatrices}.txt", "w")

for iIteration in tqdm(range(Iterations),"Computing SDP"):
    for jIteration in range(StaticIterations):
        Conditions  = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
        Conditions.variationUnitaryGramMatrix(epsilon*iIteration)
        Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)
        if (MinimumErrorLowerBound > round(Solution['SDPSolution'],3)):
            PossibleErrors.write("SDP Minimum Error result"+str(round(Solution['SDPSolution'],3))+"\nLowerBound"+str(MinimumErrorLowerBound)+"\nConditions"+str(Conditions))
            PossibleErrors.write("##############################################################################################\n")
        ExlusionSuccessProbability.append(round(Solution['SDPSolution'],4))
        DistanceFrobenius.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),2))
        Distance1.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),1))
        DistanceInfinite.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),float("inf")))
        DistanceSpectral.append(picos.SpectralNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))
        DistanceNuclear.append(picos.NuclearNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))

PossibleErrors.close()
Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(GroupGeneratedConditions)
print(round(Conditions.getPerfectExlusionLowerBoundMinimumError(),3))

plt.figure(figsize=(10, 7))
plt.scatter(DistanceFrobenius,ExlusionSuccessProbability, label = "Frobenius")
# plt.scatter(Distance1,ExlusionSuccessProbability, label = "D1")
# plt.scatter(DistanceInfinite,ExlusionSuccessProbability, label = "Dinf")
# plt.scatter(DistanceSpectral,ExlusionSuccessProbability, label = "Spectral")
# plt.scatter(DistanceNuclear,ExlusionSuccessProbability, label = "Nuclear")
plt.axhline(y = MinimumErrorLowerBound,color='blue', linestyle='--', linewidth=2, label="LowerBound")
plt.xlabel("Norm Distance")
plt.ylabel("Exclusion success probability")
plt.legend(loc='upper right')
# plt.show()
plt.savefig(f"../Plots/ExclusionNormPlots{NumberOfMatrices}.pdf")
plt.close()
