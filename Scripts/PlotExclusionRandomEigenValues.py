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
import mpltern
import picos

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
IterationsRandomSets          = 10000
DataMinimumError              = []
EigenValue0                   = []
EigenValue1                   = []
SDPMinimumErrorList           = []
LowerBound                    = []
LowerBoundZE                  = []

for iIteration in tqdm(range(IterationsRandomSets),"Computing points for the Minimum Error plot"):
    EigenValues             = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized   = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    Conditions              = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    print("UwU\n",Conditions.getGramMatrix())
    print("Yay\n",picos.trace(Conditions.getGramMatrix()))
    if Conditions.getPerfectExlusion():
        continue
    SDPMinimumError = round(QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)['SDPSolution'],8)
    EigenValue0.append(EigenValuesNormalized[0])
    EigenValue1.append(EigenValuesNormalized[1])
    SDPMinimumErrorList.append(SDPMinimumError)
    LowerBound.append(Conditions.getPerfectExlusionLowerBoundMinimumError())
    LowerBoundZE.append(Conditions.getPerfectExlusionLowerBoundZeroError())
    DataMinimumError.append((EigenValuesNormalized[0],EigenValuesNormalized[1],SDPMinimumError,Conditions.getPerfectExlusionLowerBoundMinimumError()))

# fig, ax = plt.subplots(figsize=(12, 7))
# sc = ax.scatter(EigenValue0, EigenValue1, c=LowerBound, cmap='jet', s=50)
# plt.colorbar(sc, ax=ax, label="Exclusion probability")
#
# ax.set_xlabel("$\lambda_1$")
# ax.set_ylabel("$\lambda_2$")
#
# plt.savefig(f"../Plots/ExclusionEignevaluesGroupGeneratedLowerBoundZ{NumberOfMatrices}.pdf")
# plt.close()
#
# fig, ax = plt.subplots(figsize=(12, 7))
# sc = ax.scatter(EigenValue0, EigenValue1, c=LowerBoundZE, cmap='jet', s=50)
# plt.colorbar(sc, ax=ax, label="Exclusion probability")
#
# ax.set_xlabel("$\lambda_1$")
# ax.set_ylabel("$\lambda_2$")
#
# plt.savefig(f"../Plots/ExclusionEignevaluesGroupGeneratedLowerBoundZ{NumberOfMatrices}ZeroError.pdf")
# plt.close()
#
# fig, ax = plt.subplots(figsize=(12, 7))
# sc = ax.scatter(EigenValue0, EigenValue1, c=SDPMinimumErrorList, cmap='jet', s=50)
# plt.colorbar(sc, ax=ax, label="Exclusion probability")
#
# ax.set_xlabel("$\lambda_1$")
# ax.set_ylabel("$\lambda_2$")
#
# plt.savefig(f"../Plots/ExclusionEignevalueRandom{NumberOfMatrices}.pdf")
# plt.close()
#
# for iEnsamble in DataMinimumError:
#     if round(iEnsamble[2],10)<round(iEnsamble[3],6):
#         raise ValueError(f"Something went wrong, {iEnsamble}")

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(EigenValue0, EigenValue1, SDPMinimumErrorList, color='forestgreen', marker='o', label='SDP Minimum Error', alpha=0.25)

ax.scatter(EigenValue0, EigenValue1, LowerBound, color='dodgerblue', marker='o', label='Group generated (Lower bound)', alpha=1)

ax.set_xlabel("$\lambda_1$")
ax.set_ylabel("$\lambda_2$")
ax.set_zlabel('Exclusion probability')

ax.legend()

plt.savefig(f"../Plots/ExclusionEignevalueRandom{NumberOfMatrices}3D.pdf")
plt.show()
plt.close()



import QSExSetUp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
import mpltern

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
IterationsRandomSets          = 50000
SDPMinimumError               = []
LowerBound                    = []
EigenValuesList               = []


for iIteration in tqdm(range(IterationsRandomSets),"Computing points for the Minimum Error plot"):
    EigenValues             = np.random.rand(NumberOfMatrices)
    EigenValuesNormalized   = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
    Conditions              = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    if Conditions.getPerfectExlusion():
        continue
    SDPMinimumError.append(round(QSExSetUp.SDPSolver.SolveSDPExlusionDualMinimumError(Conditions)['SDPSolution'],8))
    LowerBound.append(Conditions.getPerfectExlusionLowerBoundMinimumError())
    EigenValuesList.append(EigenValuesNormalized)

ax = plt.subplot(projection='ternary')
# ax.set_ternary_min(0,0,0)
# ax.set_ternary_max(NumberOfMatrices,NumberOfMatrices,NumberOfMatrices)

EigenValue0, EigenValue1, EigenValue2 = zip(*EigenValuesList)
scatter = ax.scatter(EigenValue0, EigenValue1, EigenValue2, marker="o",c=SDPMinimumError, cmap='jet')
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Exclusion probability")
ax.set_tlabel(f'$\lambda_0$')
ax.set_llabel(f'$\lambda_1$')
ax.set_rlabel(f'$\lambda_2$')
plt.savefig("../Plots/ExclusionTerniayMinimumErrorSDP.pdf")
plt.close()

ax = plt.subplot(projection='ternary')
# ax.set_ternary_min(0,0,0)
# ax.set_ternary_max(NumberOfMatrices,NumberOfMatrices,NumberOfMatrices)
EigenValue0, EigenValue1, EigenValue2 = zip(*EigenValuesList)
scatter = ax.scatter(EigenValue0, EigenValue1, EigenValue2, marker="o",c=LowerBound, cmap='jet')
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Exclusion probability")
ax.set_tlabel(f'$\lambda_0$')
ax.set_llabel(f'$\lambda_1$')
ax.set_rlabel(f'$\lambda_2$')
plt.savefig("../Plots/ExclusionTerniayMinimumErrorLoweBound.pdf")
plt.close()
