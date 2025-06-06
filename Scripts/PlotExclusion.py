########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     05/03/2025                                                                                #
#  Purpose:  Plot                                                                                      #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
IterationsPerfectExclusion    = 50000
DataPerfectExlusionProbability= []
DataPerfectExlusionZone       = []

for iIteration in tqdm(range(IterationsPerfectExclusion),"Computing plots for the perfect exclusion zone"):
    EigenValues                   = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    Overlap,Phase                 = Conditions.getOverlapsAndPhases()
    DataPerfectExlusionProbability.append((Overlap[1],Phase[1],Conditions.getPerfectExlusionLowerBoundMinimumError()))
    DataPerfectExlusionZone.append((Overlap[1],Phase[1],int(Conditions.getPerfectExlusion())))

Overlap_values, Phase_values, SuccessProb_values = zip(*DataPerfectExlusionProbability)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(Overlap_values, Phase_values, SuccessProb_values, c=SuccessProb_values, cmap='jet')

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
ax.set_zlabel("Exclusion probability")
# ax.set_title(f"Overlap vs. Phase vs. Probability of excluding correctly with minimum error procedure for Z{NumberOfMatrices} group generated ensambles")
plt.colorbar(sc, ax=ax, label="Exclusion probability")
plt.savefig(f"../Plots/ExclusionOverlapVSPhaseVSMinimumErrorProbabilityZ{NumberOfMatrices}3Dplot.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(12, 7))

sc = ax.scatter(Overlap_values, Phase_values, c=SuccessProb_values, cmap='jet', s=50)

plt.colorbar(sc, ax=ax, label="Exclusion probability")

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
# plt.show()
plt.savefig(f"../Plots/ExlusionOverlapVSPhaseVSMinimumErrorProbabilityZ{NumberOfMatrices}HeatMap.pdf")
plt.close()

# Overlap_values, Phase_values, SuccessProb_values = zip(*DataPerfectExlusionZone)
#
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
#
# sc = ax.scatter(Overlap_values, Phase_values, SuccessProb_values, c=SuccessProb_values, cmap='jet')
#
# ax.set_xlabel("Overlap")
# ax.set_ylabel("Phase (radians)")
# ax.set_zlabel("Perfect Exclusion Probability")
# ax.set_title(f"Z{NumberOfMatrices} - 3D Plot of Overlap vs. Phase vs. Perfect Exclusion Zone")
# plt.colorbar(sc, ax=ax, label="Perfect Exclusion")
# plt.savefig(f"../Plots/ExclusionOverlapVSPhaseVSPerfectExclussionZoneZ{NumberOfMatrices}3Dplot.pdf")
# plt.close()
#
# fig, ax = plt.subplots(figsize=(8, 6))
#
# sc = ax.scatter(Overlap_values, Phase_values, c=SuccessProb_values, cmap='jet', s=50)
#
# plt.colorbar(sc, ax=ax, label="Perfect Exclusion Zone")
#
# ax.set_xlabel("Overlap")
# ax.set_ylabel("Phase (radians)")
# ax.set_title(f"Z{NumberOfMatrices} - Heat Map Overlap vs. Phase vs. Perfect Exclusion Zone")
# plt.savefig(f"../Plots/OverlapVSPhaseVSPerfectExclussionZoneZ{NumberOfMatrices}HeatMap.pdf")
# plt.close()
#
# plt.show()
#################################################
Overlap_values, Phase_values, SuccessProb_values = zip(*DataPerfectExlusionZone)

plt.figure(figsize=(10, 7))

Overlap_values0 = [Overlap_values[iDot] for iDot in range(IterationsPerfectExclusion) if SuccessProb_values[iDot] == 0]
Phase_values0   = [Phase_values[iDot]   for iDot in range(IterationsPerfectExclusion) if SuccessProb_values[iDot] == 0]
Overlap_values1 = [Overlap_values[iDot] for iDot in range(IterationsPerfectExclusion) if SuccessProb_values[iDot] == 1]
Phase_values1   = [Phase_values[iDot]   for iDot in range(IterationsPerfectExclusion) if SuccessProb_values[iDot] == 1]


plt.scatter(Overlap_values0, Phase_values0, c='darkblue',label = 'Non-perfect exlusion zone')
plt.scatter(Overlap_values1, Phase_values1, c='darkred', label = 'Perfect exclusion zone')
plt.xlabel("Overlap")
plt.ylabel("Phase (radians)")
plt.legend()
plt.savefig(f"../Plots/ExlusionOverlapVSPhaseVSPerfectExclussionZoneZ{NumberOfMatrices}BinaryHeatMap.pdf")
# plt.show()
plt.close()
