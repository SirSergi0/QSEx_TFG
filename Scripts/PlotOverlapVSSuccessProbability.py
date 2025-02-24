########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Plot Overlap VS Succes Probability for group generated states                             #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


############################ SDP vs SRM #################################

Group                         = 3
MatrixGenerationMethod        = "ZnOverlap"
ProbabliltiesGenerationMethod = "Equal"
Acuracy                       = 100
PhaseList                     = [0]
OverlapLimit                  = 1
NumberOfMatrices              = Group
MatrixDimension               = Group
Overlap                       = 0
SuccesProbabilitySRM          = []
SuccesProbabilityMinimumError = []
SuccesProbabilityZeroError    = []
OverlapsList                  = []

for iPoint in tqdm(range(1,Acuracy),"Computing the points"):
    Overlap     += OverlapLimit/Acuracy
    OverlapList = [Overlap]
    GenerationConditions = QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsOverlap(OverlapList, PhaseList, NumberOfMatrices)
    Conditions  = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,
                                                                         MatrixDimension,
                                                                         MatrixGenerationMethod, GenerationConditions, 
                                                                         ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
    SolutionMinimumError = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)
    SolutionZeroError    = QSExSetUp.SDPSolver.SolveSDPDualZeroError(Conditions)

    SuccesProbabilitySRM.append(Conditions.getSRMSuccessProbability())
    OverlapsList.append(Overlap)
    SuccesProbabilityMinimumError.append(round(SolutionMinimumError['SDPSolution'],8))
    SuccesProbabilityZeroError.append(round(SolutionZeroError['SDPSolution'],8))

plt.plot(OverlapsList, SuccesProbabilitySRM, label = "SRM", color = "darkorange")
plt.scatter(OverlapsList, SuccesProbabilityMinimumError, label = "Minimum Error SDP",marker = "o", color = "mediumslateblue")
plt.scatter(OverlapsList, SuccesProbabilityZeroError, label = "Zero Error SDP",marker = "o", color = "forestgreen")
plt.xlabel("Overlap")
plt.ylabel("Success Probability")
plt.title(f"SDP vs SRM Zn{Group}")
plt.legend(loc = 3)
plt.savefig(f"../Plots/OverlapVSSucessProbabilitySDPvsSRM{MatrixGenerationMethod}{Group}Phase{PhaseList[0]}p.pdf")

############################ 3D Zn #################################
#
# import QSExSetUp
# import matplotlib.pyplot as plt
# import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# 
# NumberOfMatrices              = 3
# MatrixDimension               = 3
# MatrixGenerationMethod        = "Zn"
# ProbabliltiesGenerationMethod = "Equal"
#
# OverlapPrecission = 50
# PhasePrecission = 50
#
# Overlap = 0
#
# Data = []
#
# # for i, iGroup in enumerate(tqdm(range(FirstGroup, LastGroup + 1), desc="Calculating...")):
# for iOverlapIteration in tqdm(range(1,OverlapPrecission), desc="Computating..."):
#     Overlap     += 1/OverlapPrecission
#     OverlapList  = [Overlap]
#     Phase        = 0
#     for iPhaseIteration in range(1,PhasePrecission):
#         Phase               += 2*np.pi/PhasePrecission
#         PhaseList            = [Phase]
#         GenerationConditions = QSExSetUp.GramGeneratedStates.ZnGramMatrixConditions(OverlapList,PhaseList, NumberOfMatrices)
#         Conditions           = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, GenerationConditions, ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)
#         Solution = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)
#         Data.append((Overlap,Phase,round(Solution['SDPSolution'],4)))
#
# Overlap_values, Phase_values, SuccessProb_values = zip(*Data)
#
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
#
# sc = ax.scatter(Overlap_values, Phase_values, SuccessProb_values, c=SuccessProb_values, cmap='viridis')
#
# ax.set_xlabel("Overlap")
# ax.set_ylabel("Phase (radians)")
# ax.set_zlabel("Success Probability")
# ax.set_title("3D Plot of Overlap vs. Phase vs. Success Probability")
# plt.colorbar(sc, ax=ax, label="success probability")
# plt.show()
#
