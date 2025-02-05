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

FirstGroup                    = 2
LastGroup                     = 10
MatrixGenerationMethod        = "Zn"
ProbabliltiesGenerationMethod = "Equal"
Acuracy                       = 100

colors = [plt.cm.rainbow(i / (LastGroup-FirstGroup)) for i in range(LastGroup+1-FirstGroup)]

for i,iGroup in enumerate(range(FirstGroup,LastGroup+1)):
    NumberOfMatrices              = iGroup
    MatrixDimension               = iGroup
    Overlap                       = 0
    SuccesProbabilityList         = []
    OverlapsList                  = []

    for iPoint in range (1,Acuracy):
        Overlap += 1/Acuracy
        Conditions = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod, Overlap = Overlap)
        SuccesProbabilityList.append(Conditions.getSRMSuccessProbability())
        OverlapsList.append(Overlap)
    plt.plot(OverlapsList, SuccesProbabilityList, label = f'{MatrixGenerationMethod[0]}{iGroup}', color = colors[i])

plt.xlabel("Overlap")
plt.ylabel("Success Probability")
plt.title(f"Group generated {MatrixGenerationMethod}")
plt.legend()
plt.savefig(f"../Plots/OverlapVSSucessProbability{MatrixGenerationMethod}from{FirstGroup}to{LastGroup}.pdf")

