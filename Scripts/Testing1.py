########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import QSExSetUp

NumberOfMatrices              = 3
MatrixDimension               = 3
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [1,2,3]
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

print(Conditions)

print("\n####################################\n\
####### DISCRIMINATION SDP #########\n\
####################################\n")

print("------------------------\n\
Computing the primal SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPMinimumError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("------------------------\n\
Computing the Dual SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPDualMinimumError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("----------------------------\n\
Computing the Zero Error SDP\n\
----------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPZeroError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("---------------------------------\n\
Computing the Dual SDP Zero Error\n\
---------------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPDualZeroError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

print("\n####################################\n\
########## EXCLUSION SDP ###########\n\
####################################\n")

print("------------------------\n\
Computing the primal SDP\n\
------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPExclusionMinimumError(Conditions)

print("The given problem has been:\n", Solution['SDPSolution'])

print("The success probability is: ", round(Solution['SDPSolution'],4))

for i,iPOVM in enumerate(Solution['POVMs']):
    print(f"POVM {i}:\n",iPOVM)
print("Sum of POVM (Identity):\n")
print(sum(Solution['POVMs']))

# print("------------------------\n\
# Computing the Dual SDP\n\
# ------------------------")
#
# Solution = QSExSetUp.SDPSolver.SolveSDPExclusionDualMinimumError(Conditions)
#
# print("The given problem has been:\n", Solution['SDPSolution'])
#
# print("The success probability is: ", round(Solution['SDPSolution'],4))
