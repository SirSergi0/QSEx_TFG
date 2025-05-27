########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     23/05/2025                                                                                #
#  Purpose:  Plot                                                                                      #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import mpltern
import picos
from collections import Counter

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
IterationsRandomSets          = 10
DataMinimumError              = []
EigenValue0                   = []
EigenValue1                   = []
SDPMinimumErrorList           = []
SDPZeroErrorList              = []
LowerBound                    = []
LowerBoundZE                  = []
data                          = []
for iIteration in range (500):
    EigenValues = [1 for i in range (NumberOfMatrices-1)]
    EigenValues.append(iIteration^3)
    # EigenValues             = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized   = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    # while (Conditions.getPerfectExlusion()):
    #     EigenValues             = np.random.rand(NumberOfMatrices)
    #     EigenValuesNormalized   = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    #     print("UwU")
    #     Conditions              = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

    for iIteration in tqdm(range(IterationsRandomSets),f"{iIteration}Computing points for the Minimum Error plot"):
        Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
        if Conditions.getPerfectExlusion():
            continue

        UwU = round(QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)['SDPSolution'],4)
        if (UwU < round(Conditions.getPerfectExlusionLowerBoundMinimumError(),4)):
            raise ValueError ("SDP result", UwU,"LowerBound", Conditions.getPerfectExlusionLowerBoundMinimumError(),"Conditions",Conditions)
        data.append([Conditions.getPerfectExlusionLowerBoundMinimumError(),UwU])
        # SDPZeroErrorList.append(round(QSExSetUp.SDPSolver.SolveSDPExclusionZeroError(Conditions)['SDPSolution'],4))
    # print((Conditions.getPerfectExlusionLowerBoundMinimumError(),SDPMinimumErrorList))
    # data.append((SDPMinimumErrorList,Conditions.getPerfectExlusionLowerBoundMinimumError()))

# plt.figure(figsize=(10, 6)) 
# plt.hist(SDPMinimumErrorList, bins=100, color = "orange",  label  = f"Number of entries: {len(SDPMinimumErrorList)}") 
# plt.axvline(x = Conditions.getPerfectExlusionLowerBoundMinimumError(), color='blue', linestyle='--', linewidth=2, label="LowerBound")
# plt.xlabel('Success probability minimum error')
# plt.ylabel('Frequency')
# plt.legend()
# plt.savefig(f"../Plots/ExclusionMinimumErrorRandomDistributionZ{NumberOfMatrices}Prob{round(Conditions.getPerfectExlusionLowerBoundMinimumError(),3)}.pdf")
# plt.close()
# plt.figure(figsize=(10, 6)) 
# plt.hist(SDPZeroErrorList, bins=100, color = "orange",  label  = f"Number of entries: {len(SDPZeroErrorList)}") 
# plt.axvline(x = Conditions.getPerfectExlusionLowerBoundZeroError(),color='blue', linestyle='--', linewidth=2, label="LowerBound")
# plt.xlabel('Success probability minimum error')
# plt.ylabel('Frequency')
# plt.legend()
# plt.savefig(f"../Plots/ExclusionZeroErrorRandomDistributionZ{NumberOfMatrices}Prob{round(Conditions.getPerfectExlusionLowerBoundZeroError(),3)}.pdf")
# plt.close()



# Convert to NumPy array
data = np.array(data)
x = data[:, 0]
y = data[:, 1]

# Create 2D histogram
x_bins = 100
y_bins = 100
heatmap, xedges, yedges = np.histogram2d(x, y, bins=[x_bins, y_bins])

# Plot
plt.figure(figsize=(8, 6))
plt.imshow(
    heatmap.T,
    origin='lower',
    aspect='auto',
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    cmap='jet'
)
plt.colorbar(label='Frequency')
plt.xlabel('lower bound success probability')
plt.ylabel('SDP Success probability')
plt.show()
