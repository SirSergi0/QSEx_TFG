########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     24/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import QSExSetUp
import picos
import numpy as np
from scipy.linalg import sqrtm

NumberOfMatrices              = 4
MatrixDimension               = 3
MatrixGenerationMethod        = "GroupGeneratedStates"
ProbabliltiesGenerationMethod = "Equal"
initialState                  = picos.Constant([1,0,0])
sigmaX                        = picos.Constant([[0,1,0],[1,0,0],[0,0,1]])
print(len(sigmaX.shape))
Conditions                    = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension,MatrixGenerationMethod, ProbabliltiesGenerationMethod, seedState = initialState, involutoryMatrix = sigmaX)
DensityMatrices               = Conditions.getDensityMatrices()

print(Conditions)

for iMatix in range(len(DensityMatrices)):
    print(f"Matrix number {iMatix}:\n",DensityMatrices[iMatix])

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

print("-----------------------\n\
Pretty Good Measurement\n\
-----------------------")

GramMatrix         = Conditions.getGramMatrixWithPriors()
SquareRoot         = picos.Constant(sqrtm(GramMatrix))
SquareRootDiagonal = np.diagonal(SquareRoot.value)

sumSquare = 0
for iElement in SquareRootDiagonal:
    sumSquare += abs(iElement)**2

print("The GramMatrix with priors is         :\n",GramMatrix)

print("The square of the squareRoot Matrix is:\n", SquareRoot)

print("Computing G-S^2=0?                    :\n", GramMatrix-(SquareRoot*SquareRoot))

print("The sum of the diagonalSquare is      :", sumSquare)

print("----------------------------\n\
Computing the Zero Error SDP\n\
----------------------------")

Solution = QSExSetUp.SDPSolver.SolveSDPZeroError(Conditions)

print("The given problem has been  :\n", Solution['SDPSolution'])

print("The POVMs are               :\n", Solution['POVMs'])

print("Uncertainity probability is :\n", round(Solution['SDPSolution'],4))

print("The probability of perfect exclussion is:\n", round(1-float(Solution['SDPSolution']),4))


UnknownPOVM = Solution["UncertainPOVM"]

print("We do not now what is going on if the POVM is:\n",UnknownPOVM)

for iElement, iPOVM in enumerate(Solution['POVMs']):
    print(f"POVM {iElement}:")
    print(iPOVM)
    UnknownPOVM += iPOVM
    print(f"DensityMatrices{iElement}:\n",picos.trace(DensityMatrices[iElement]*iPOVM))

print("Does the sum equal the identity?\n",UnknownPOVM)

