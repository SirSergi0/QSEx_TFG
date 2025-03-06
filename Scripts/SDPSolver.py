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
import GramGeneratedStates

def SolveSDPMinimumError(Conditions):
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
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

# def SolveSDPMinimumError(Conditions):
#     if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
#         raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
#     NumberOfMatrices     = Conditions.getNumberOfMatrices()
#     GramMatrixWithPriors = Conditions.getGramMatrixWithPriors()
#     MySDP                = picos.Problem()
#     XMatrix              = picos.ComplexVariable('XMatrix', GramMatrixWithPriors.shape)
#     SuccessProbability   = 0
#
#     for iElement in range(NumberOfMatrices):
#         OrthonormalVectorList = []
#         for jElement in range(NumberOfMatrices):
#             if (iElement == jElement):
#                 OrthonormalVectorList.append(1)
#             else:
#                 OrthonormalVectorList.append(0)
#         OrthonormalVector = picos.Constant(OrthonormalVectorList)
#         SuccessProbability += (OrthonormalVector.H*XMatrix*OrthonormalVector)*(OrthonormalVector.H*XMatrix.H*OrthonormalVector)
#     MySDP.set_objective("max", SuccessProbability)
#     MySDP.add_constraint(XMatrix.H*XMatrix == GramMatrixWithPriors)
#     MySDP.solve(solver = 'cvxopt')
#
#     return {'SDPSolution' : MySDP}


def SolveSDPDualMinimumError(Conditions):
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
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
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")

    NumberOfMatrices     = Conditions.getNumberOfMatrices()
    GramMatrixWithPriors = Conditions.getGramMatrixWithPriors()
    MySDP                = picos.Problem()
    GammaMatrix          = picos.HermitianVariable("GammaMatrix", GramMatrixWithPriors.shape)
    
    MySDP.set_objective("max", picos.trace(GammaMatrix))
    MySDP.add_constraint(GammaMatrix>>0)
    GammaMatrixDiagonalList = []
    for iElement in range(NumberOfMatrices):
        OrthonormalVectorList = []
        for jElement in range(NumberOfMatrices):
            if (iElement == jElement):
                OrthonormalVectorList.append(1)
            else:
                OrthonormalVectorList.append(0)
        OrthonormalVector = picos.Constant(OrthonormalVectorList)
        GammaMatrixDiagonalList.append((OrthonormalVector.H*GammaMatrix*OrthonormalVector)*OrthonormalVector*OrthonormalVector.H)
    MySDP.add_constraint(sum(GammaMatrixDiagonalList)<<GramMatrixWithPriors)
    
    MySDP.solve(solver="cvxopt")
    return {'SDPSolution' : MySDP}

def SolveSDPDualZeroError(Conditions):
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
    GramMatrixWithPriors = Conditions.getGramMatrixWithPriors()
    NumberOfMatrices     = Conditions.getNumberOfMatrices()
    ZMatrix              = picos.HermitianVariable("ZMatrix", GramMatrixWithPriors.shape)
    MySDP                = picos.Problem()
    
    MySDP.set_objective("min",picos.trace(GramMatrixWithPriors*ZMatrix))
    MySDP.add_constraint(ZMatrix >> 0) 
    for iElement in range(NumberOfMatrices):
        OrthonormalVectorList = []
        for jElement in range(NumberOfMatrices):
            if (iElement == jElement):
                OrthonormalVectorList.append(1)
            else:
                OrthonormalVectorList.append(0)
        OrthonormalVector = picos.Constant(OrthonormalVectorList)
        MySDP.add_constraint(OrthonormalVector.H*ZMatrix*OrthonormalVector >> 1)
    
    MySDP.solve(solver="cvxopt")
    
    return {'SDPSolution' : MySDP }

def SolveSDPExclusionMinimumError(Conditions):
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
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

    MySDP.set_objective("min", SuccessProbability)

    for iPOVM in range(NumberOfMatrices):
        MySDP.add_constraint(POVMlist[iPOVM]>>0)

    I = np.eye(NumberOfDimensions) 
    MySDP.add_constraint(sum(POVMlist) == picos.Constant(I))

    MySDP.solve(solver="cvxopt")

    return {'SDPSolution' : MySDP, 'POVMs' : POVMlist}

def SolveSDPExclusionDualMinimumError(Conditions):
    if not isinstance(Conditions,(DensityMatricesAndPriorsClass.DensityMaticesAndPriors,GramGeneratedStates.GramGeneratedStatesClass)):
        raise TypeError("The given conditions must be a DensityMatricesAndPriorsClass.DensityMaticesAndPriors or GramGeneratedStates.GramGeneratedStatesClass")
    MyDensityMatrices      = Conditions.getDensityMatrices()
    MyPriorsPropbabilities = Conditions.getPriorProbabilities()
    NumberOfMatrices       = Conditions.getNumberOfMatrices()

    MySDP                  = picos.Problem()
    LagrangeMultiplierY    = picos.HermitianVariable("LagrangeMultiplier", MyDensityMatrices[0].shape)
    ErrorProbability       = picos.trace(LagrangeMultiplierY)

    MySDP.set_objective("max", ErrorProbability)
    
    for iContrain in range(NumberOfMatrices):
        MySDP.add_constraint(LagrangeMultiplierY - MyPriorsPropbabilities[iContrain]*MyDensityMatrices[iContrain]>>0)

    MySDP.solve(solver="cvxopt")

    return {'SDPSolution' : MySDP}
