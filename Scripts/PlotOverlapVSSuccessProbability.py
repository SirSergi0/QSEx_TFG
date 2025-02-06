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

############################## Zn Plot ##############################
FirstGroup                    = 2
LastGroup                     = 10
MatrixGenerationMethod        = "Zn"
ProbabliltiesGenerationMethod = "Equal"
Acuracy                       = 100

colors = [plt.cm.rainbow(i / (LastGroup-FirstGroup)) for i in range(LastGroup+1-FirstGroup)]

for i, iGroup in enumerate(tqdm(range(FirstGroup, LastGroup + 1), desc="Calculating...")):
    NumberOfMatrices              = iGroup
    MatrixDimension               = iGroup
    Overlap                       = 0
    SuccesProbabilitySRM          = []
    SuccesProbabilitySDP          = []
    OverlapsList                  = []

    for iPoint in range (1,Acuracy):
        Overlap     += 1/Acuracy
        Conditions  = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,
                                                                             MatrixDimension,
                                                                             MatrixGenerationMethod, 
                                                                             ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod, 
                                                                             Overlap = Overlap)
        Solution    = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)
        SuccesProbabilitySRM.append(Conditions.getSRMSuccessProbability())
        SuccesProbabilitySDP.append(round(Solution['SDPSolution'],8))
        OverlapsList.append(Overlap)
    if i == 0: plt.scatter(OverlapsList, SuccesProbabilitySDP, label = "SDP",marker = ".", color = "grey")
    else     : plt.scatter(OverlapsList, SuccesProbabilitySDP,marker = ".", color = "grey")
    plt.plot(OverlapsList, SuccesProbabilitySRM, label = f'SRM{MatrixGenerationMethod[0]}{iGroup}', color = colors[i])

plt.xlabel("Overlap")
plt.ylabel("Success Probability")
plt.title(f"Group generated {MatrixGenerationMethod}")
plt.legend()
plt.savefig(f"../Plots/OverlapVSSucessProbability{MatrixGenerationMethod}from{FirstGroup}to{LastGroup}.pdf")

############################ SDP vs SRM #################################
#
# Group                         = 3
# MatrixGenerationMethod        = "Zn"
# ProbabliltiesGenerationMethod = "Equal"
# Acuracy                       = 100
# NumberOfMatrices              = Group
# MatrixDimension               = Group
# Overlap                       = 0
# SuccesProbabilitySRM          = []
# SuccesProbabilitySDP          = []
# OverlapsList                  = []
#
# for iPoint in range (1,Acuracy):
#     Overlap     += 1/Acuracy
#     Conditions  = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,
#                                                                          MatrixDimension,
#                                                                          MatrixGenerationMethod, 
#                                                                          ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod, 
#                                                                          Overlap = Overlap)
#     Solution    = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)
#     SuccesProbabilitySRM.append(Conditions.getSRMSuccessProbability())
#     OverlapsList.append(Overlap)
#     SuccesProbabilitySDP.append(round(Solution['SDPSolution'],8))
#
# plt.plot(OverlapsList, SuccesProbabilitySRM, label = "SRM")
# plt.scatter(OverlapsList, SuccesProbabilitySDP, label = "SDP",marker = "2", color = "black")
# plt.xlabel("Overlap")
# plt.ylabel("Success Probability")
# plt.title(f"SDP vs SRM {MatrixGenerationMethod}{Group}")
# plt.legend()
# plt.savefig(f"../Plots/OverlapVSSucessProbabilitySDPvsSRM{MatrixGenerationMethod}{Group}.pdf")
#
#
