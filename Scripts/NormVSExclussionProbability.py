########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     05/04/2025                                                                                #
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
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [1 for i in range(NumberOfMatrices - 1)]
EigenValues.append(100)
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
Iterations                    = 1000
GroupGeneratedConditions      = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, "ZnEigenValues", QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)


DistanceFrobenius             = []
Distance1                     = []
DistanceInfinite              = []
DistanceSpectral              = []
DistanceNuclear               = []
ExlusionSuccessProbability    = []

for iIteration in tqdm(range(Iterations),"Computing SDP"):
    Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)
    ExlusionSuccessProbability.append(round(Solution['SDPSolution'],4))
    DistanceFrobenius.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),2))
    Distance1.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),1))
    DistanceInfinite.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix(),float("inf")))
    DistanceSpectral.append(picos.SpectralNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))
    DistanceNuclear.append(picos.NuclearNorm(GroupGeneratedConditions.getUnitaryMatrix() - Conditions.getUnitaryMatrix()))

plt.figure(figsize=(10, 7))
plt.scatter(DistanceFrobenius,ExlusionSuccessProbability, label = "Frobenius")
plt.scatter(Distance1,ExlusionSuccessProbability, label = "D1")
plt.scatter(DistanceInfinite,ExlusionSuccessProbability, label = "Dinf")
plt.scatter(DistanceSpectral,ExlusionSuccessProbability, label = "Spectral")
plt.scatter(DistanceNuclear,ExlusionSuccessProbability, label = "Nuclear")
plt.xlabel("Norm Distance")
plt.ylabel("Exclusion success probability")
plt.legend(loc='upper right')
plt.show()
plt.savefig(f"../Plots/ExclusionNormPlots{NumberOfMatrices}.pdf")
plt.close()
