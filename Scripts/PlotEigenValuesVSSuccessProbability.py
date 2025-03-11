########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     22/02/2025                                                                                #
#  Purpose:  Plot                                                                                      #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

######################## Minimum Error Anlaysis ######################## 

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
IterationsMinimumError        = 50000
IterationsZeroError           = 50000
DataMinimumError              = []
DataZeroError                 = []

for iIteration in tqdm(range(IterationsMinimumError),"Computing points for the Minimum Error plot"):
    EigenValues                   = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    Overlap,Phase                 = Conditions.getOverlapsAndPhases()
    DataMinimumError.append((Overlap[1],Phase[1],Conditions.getSRMSuccessProbability()))

Overlap_values, Phase_values, SuccessProb_values = zip(*DataMinimumError)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(Overlap_values, Phase_values, SuccessProb_values, c=SuccessProb_values, cmap='jet')

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
ax.set_zlabel("Discrimination probability")
# ax.set_title(f"Discrimination - Z{NumberOfMatrices} - 3D Plot of Overlap vs. Phase vs. Success Probability Minimum Error")
plt.colorbar(sc, ax=ax, label="Discrimination probability")
plt.savefig(f"../Plots/DiscriminationOverlapVSPhaseVSMinimumErrorProbabilityZ{NumberOfMatrices}3Dplot.pdf")

fig, ax = plt.subplots(figsize=(12, 7))

sc = ax.scatter(Overlap_values, Phase_values, c=SuccessProb_values, cmap='jet', s=50)

plt.colorbar(sc, ax=ax, label="Discrimination probability")

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
# ax.set_title(f"Discrimination - Z{NumberOfMatrices} - Heat Map Overlap vs. Phase vs. Success Probability Minimum Error")

plt.savefig(f"../Plots/DiscriminationOverlapVSPhaseVSMinimumErrorProbabilityZ{NumberOfMatrices}HeatMap.pdf")

# target_phase = 0
# tolerance = 0.01
#
# Phase_values = np.array(Phase_values)
# fileteredValues = (Phase_values > target_phase - tolerance) & (Phase_values < target_phase + tolerance)
#
# Overlap_filtered = np.array(Overlap_values)[fileteredValues]
# SuccessProb_filtered = np.array(SuccessProb_values)[fileteredValues]
#
# fig, ax = plt.subplots(figsize=(8, 6))
#
# ax.scatter(Overlap_filtered, SuccessProb_filtered)
#
# ax.set_xlabel("Overlap")
# ax.set_ylabel("Success Probability")
# ax.set_title(f"Cross-Section for Phase {target_phase:.2f} radians")
#
# plt.savefig(f"../Plots/OverlapVSPhaseVSSuccessProbabilityZ{NumberOfMatrices}CorossSectionPhase{target_phase}.pdf")
#
############################# Zero Error VS Minimum Error ###########################

for iIteration in tqdm(range(IterationsZeroError),"Computing points for the Zero Error plot"):
    EigenValues                   = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    Overlap,Phase                 = Conditions.getOverlapsAndPhases()
    SolutionZeroError             = QSExSetUp.SDPSolver.SolveSDPZeroError(Conditions)
    DataZeroError.append((Overlap[1],Phase[1],round(SolutionZeroError['SDPSolution'],8)))

Overlap_values, Phase_values, SuccessProb_values = zip(*DataZeroError)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(Overlap_values, Phase_values, SuccessProb_values, c=SuccessProb_values, cmap='jet')

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
ax.set_zlabel("Discrimination probability")
# ax.set_title(f"Z{NumberOfMatrices} - 3D Plot of Overlap vs. Phase vs. Success Probability Zero Error")
plt.colorbar(sc, ax=ax, label="Discrimination probability")
plt.savefig(f"../Plots/DiscriminationOverlapVSPhaseVSSuccessProbabilityZeroErrorZ{NumberOfMatrices}3Dplot.pdf")

fig, ax = plt.subplots(figsize=(12, 7))

sc = ax.scatter(Overlap_values, Phase_values, c=SuccessProb_values, cmap='jet', s=50)

plt.colorbar(sc, ax=ax, label="Discrimination probability")

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
# ax.set_title(f"Z{NumberOfMatrices} - Heat Map Overlap vs. Phase vs. Success Probability Zero Error")

plt.savefig(f"../Plots/DiscriminationOverlapVSPhaseVSSuccessProbabilityZeroErrorZ{NumberOfMatrices}HeatMap.pdf")

# target_phase = 0
# tolerance = 0.01
#
# Phase_values = np.array(Phase_values)
# fileteredValues = (Phase_values > target_phase - tolerance) & (Phase_values < target_phase + tolerance)
#
# Overlap_filtered = np.array(Overlap_values)[fileteredValues]
# SuccessProb_filtered = np.array(SuccessProb_values)[fileteredValues]
#
# fig, ax = plt.subplots(figsize=(8, 6))
#
# ax.scatter(Overlap_filtered, SuccessProb_filtered)
#
# ax.set_xlabel("Overlap")
# ax.set_ylabel("Success Probability")
# ax.set_title(f"Cross-Section for Phase {target_phase:.2f} radians")
#
# plt.savefig(f"../Plots/OverlapVSPhaseVSSuccessProbabilityZ{NumberOfMatrices}CorossSectionPhase{target_phase}.pdf")
#
