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
    DataPerfectExlusionProbability.append((Overlap[1],Phase[1],Conditions.getPerfectExlusionLowerBoundZeroError()))
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
plt.savefig(f"../Plots/ExclusionOverlapVSPhaseVSZeroErrorProbabilityZ{NumberOfMatrices}3Dplot.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(12, 7))

sc = ax.scatter(Overlap_values, Phase_values, c=SuccessProb_values, cmap='jet', s=50)

plt.colorbar(sc, ax=ax, label="Exclusion probability")

ax.set_xlabel("Overlap")
ax.set_ylabel("Phase (radians)")
# plt.show()
plt.savefig(f"../Plots/ExlusionOverlapVSPhaseVSZeroErrorProbabilityZ{NumberOfMatrices}HeatMap.pdf")
plt.close()
