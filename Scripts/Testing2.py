########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import QSExSetUp

NumberOfMatrices              = 4
MatrixDimension               = 4
MatrixGenerationMethod        = "ZnEigenValues"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [1,1,1,1]
Conditions                    = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValues, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

print(Conditions)

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

