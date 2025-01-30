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

def SolveSDPMinimumError(Conditions):
    if not isinstance(Conditions,DensityMatricesAndPriorsClass.DensityMaticesAndPriors):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors")
    MyDensityMatrices      = Conditions.getDensityMatrices()
    MyPriorsPropbabilities = Conditions.getPriorProbabilities()
    NumberOfMatrices       = Conditions.getNumberOfMatrices()
    NumberOfDimensions     = Conditions.getMatrixDimension()

    MySDP = picos.Problem()

    POVMlist = []
    for iPOVM in range(NumberOfMatrices):
        POVM = picos.HermitianVariable(f"POVM_{iPOVM}", MyDensityMatrices[iPOVM].shape)
        POVMlist.append(POVM)

    SuccessProbability = 0

    for iPOVM in range(NumberOfMatrices):
        SuccessProbability += MyPriorsPropbabilities[iPOVM]*picos.trace(POVMlist[iPOVM]*MyDensityMatrices[iPOVM])

    MySDP.set_objective("max", SuccessProbability)

    for iPOVM in range(NumberOfMatrices):
        MySDP.add_constraint(POVMlist[iPOVM]>>0)

    I = np.eye(NumberOfDimensions) 
    MySDP.add_constraint(sum(POVMlist) == picos.Constant(I))

    MySDP.solve(solver="cvxopt")

    return {'SDPSolution' : MySDP, 'POVMs' : POVMlist}

def SolveSDPDualMinimumError(Conditions):
    if not isinstance(Conditions,DensityMatricesAndPriorsClass.DensityMaticesAndPriors):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors")
    MyDensityMatrices      = Conditions.getDensityMatrices()
    MyPriorsPropbabilities = Conditions.getPriorProbabilities()
    NumberOfMatrices       = Conditions.getNumberOfMatrices()

    MySDP                  = picos.Problem()
    LagrangeMultiplierY    = picos.HermitianVariable("LagrangeMultiplier", MyDensityMatrices[0].shape)
    ErrorProbability       = picos.trace(LagrangeMultiplierY)

    MySDP.set_objective("min", ErrorProbability)
    
    for iContrain in range(NumberOfMatrices):
        MySDP.add_constraint(LagrangeMultiplierY - MyPriorsPropbabilities[iContrain]*MyDensityMatrices[iContrain]>>0)

    MySDP.solve(solver="cvxopt")

    return {'SDPSolution' : MySDP }

def SolveSDPZeroError(Conditions):
    if not isinstance(Conditions,DensityMatricesAndPriorsClass.DensityMaticesAndPriors):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors")
    MyDensityMatrices      = Conditions.getDensityMatrices()
    MyPriorsPropbabilities = Conditions.getPriorProbabilities()
    NumberOfMatrices       = Conditions.getNumberOfMatrices()
    NumberOfDimensions     = Conditions.getMatrixDimension()

    MySDP                  = picos.Problem()
    
    POVMlist = []
    for iPOVM in range(NumberOfMatrices):
        POVM = picos.HermitianVariable(f"POVM_{iPOVM}", MyDensityMatrices[iPOVM].shape)
        POVMlist.append(POVM)

    UncertainPOVM          = picos.HermitianVariable(f"UncertainPOVM", MyDensityMatrices[0].shape)

    UncertainProbability   = 0

    for iElement in range(NumberOfMatrices):
        UncertainProbability += MyPriorsPropbabilities[iElement]*picos.trace(UncertainPOVM*MyDensityMatrices[iElement])

    MySDP.set_objective("min", UncertainProbability)

    I = np.eye(NumberOfDimensions) 
    MySDP.add_constraint(sum(POVMlist) + UncertainPOVM == picos.Constant(I))
    
    for iPOVM in range(NumberOfMatrices):
        MySDP.add_constraint(picos.trace(MyDensityMatrices[iPOVM] * POVMlist[iPOVM]) == 0)
        MySDP.add_constraint(POVMlist[iPOVM] >> 0)

    MySDP.add_constraint(UncertainPOVM >> 0)
    
    MySDP.solve(solver="cvxopt")
    
    return {'SDPSolution' : MySDP, 'POVMs' : POVMlist, 'UncertainPOVM' : UncertainPOVM}
