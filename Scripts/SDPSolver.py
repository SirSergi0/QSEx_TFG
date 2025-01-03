########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Store the SDP solver function                                                             #
#                                                                                                      #
########################################################################################################

import picos
import numpy as np
import DensityMatricesAndPriorsClass

def SolveSDP(Conditions):
    if not isinstance(Conditions,DensityMatricesAndPriorsClass.DensityMaticesAndPriors):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors")
    MyDenisityMatrices     = Conditions.getDesityMatrices()
    MyPriorsPropbabilities = Conditions.getPriorProbabilities()
    NumberOfMatrices       = Conditions.getNumberOfMatrices()
    NumberOfDimensions     = Conditions.getMatrixDimension()

    MySDP = picos.Problem()

    POVMlist = []
    for iPOVM in range(NumberOfMatrices):
        POVM = picos.HermitianVariable(f"POVM_{iPOVM}", MyDenisityMatrices[iPOVM].shape)
        POVMlist.append(POVM)

    SuccessProbability = 0

    for iPOVM in range(NumberOfMatrices):
        SuccessProbability += MyPriorsPropbabilities[iPOVM]*picos.trace(POVMlist[iPOVM]*MyPriorsPropbabilities[iPOVM])

    MySDP.set_objective("max", SuccessProbability)

    for iPOVM in range(NumberOfMatrices):
        MySDP.add_constraint(POVMlist[iPOVM]>>0)

    I = np.eye(NumberOfDimensions) 
    MySDP.add_constraint(sum(POVMlist) == picos.Constant(I))

    MySDP.solve(solver="cvxopt")

    return {'SDPSolution' : MySDP, 'POVMs' : POVMlist}
