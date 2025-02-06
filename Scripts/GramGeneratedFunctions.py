########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     04/02/2025                                                                                #
#  Purpose:  Store functions used in the Gram Generated States Class                                   #
#                                                                                                      #
########################################################################################################

import numpy as np
import picos
from scipy.linalg import sqrtm

def Hermitian(matrix, SilentMode = True):
    if not isinstance(matrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given matrix is not valid for computing the hermicity")
    if np.allclose(matrix.value, matrix.H.value):
        if not SilentMode : print("The matrix is Hermitian")
        return True
    if not SilentMode : print("The matrix is NOT hermitian")
    return False

def SetOfRandomProbabilities (NumberOfMatrices):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    priorProbabilities = np.random.rand(NumberOfMatrices)
    priorProbabilities /= sum(priorProbabilities)
    return priorProbabilities

def SetOfEqualProbabilities (NumberOfMatrices):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    return np.array([1/NumberOfMatrices for iMatrix in range (NumberOfMatrices)])

def SetOfMatrices (MatrixMethod, NumberOfMatrices, MatrixDimension, seedState, involutionalMatrix):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    if not isinstance(MatrixMethod, str):
        raise TypeError("The given MatrixMethod must be a string")
    ValidMethods = ["RandomMixedStates", "RandomPureStates", "GroupGeneratedStates"]
    if not MatrixMethod in ValidMethods:
        raise ValueError(f"The given MatrixMethod is not a valid method, the valid methods are: {ValidMethods}")
    if MatrixMethod == "RandomMixedStates":
        return SetOfRandomDensityMatricesMixed(NumberOfMatrices,MatrixDimension)
    if MatrixMethod == "RandomPureStates":
        return SetOfRandomDensityMatricesPure(NumberOfMatrices, MatrixDimension)
    if MatrixMethod == "GroupGeneratedStates":
        return setOfGropGeneratedDensityMatrices(seedState, NumberOfMatrices, involutionalMatrix)
    raise ValueError("The SetOfMatrices function has not worked")

def SetOfProbabilities (ProbabilitiesMethod, NumberOfMatrices):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    if not isinstance(ProbabilitiesMethod, str):
        raise TypeError("The given Method must be a string")
    if not ProbabilitiesMethod in ["Random","Equal"]:
        raise ValueError("The given ProbabilitiesMethod is not a valid method")
    if ProbabilitiesMethod == "Random":
        return SetOfRandomProbabilities(NumberOfMatrices)
    if ProbabilitiesMethod == "Equal":
        return SetOfEqualProbabilities(NumberOfMatrices)
    raise ValueError("The SetOfMatrices function has not worked")


########################
# Generating functions #
########################

def GenerateGramMatrix(NumberOfMatrices, MatrixMethod, Overlap):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatrices must be an Integer")
    ValidMethods = ["Zn"]
    if not MatrixMethod in ValidMethods:
        raise ValueError(f"The given MatrixMethod is not a valid method, the valid methods are: {ValidMethods}")
    if not isinstance(Overlap, (float, int)):
        raise TypeError("The given Overlap must be an Integer or a float")
    if (Overlap > 1 or Overlap < 0):
        raise ValueError("The given Overlap must be a number between 0 and 1.")
    if MatrixMethod == "Zn":
        return ZnGramMatrix(NumberOfMatrices, Overlap)

def GenerateGramMatrixWithPriors(GramMatrix, priorProbabilities, NumberOfStates):
    GramMatrixWithPriors = np.array(GramMatrix)
    for iState in range(NumberOfStates):
        for jState in range(NumberOfStates):
            GramMatrixWithPriors[iState][jState] *= (np.sqrt(priorProbabilities[iState])*np.sqrt(priorProbabilities[jState]))
    return picos.Constant(GramMatrixWithPriors)

def ZnGramMatrix (NumberOfStates, Overlap):
    GramMatrix = []
    for iState in range(NumberOfStates):
        GramMatrixRow = []
        for jState in range(NumberOfStates):
            if iState == jState: GramMatrixRow.append(1)
            else: GramMatrixRow.append(Overlap)
        GramMatrix.append(GramMatrixRow)
    return picos.Constant(GramMatrix)

def GetGramDensityMatrices(GramMatrix, NumberOfStates):
    SquareRoot = np.array(sqrtm(GramMatrix))
    States     = []
    for iState in range(NumberOfStates):
        state  = picos.Constant(SquareRoot[:,iState])
        States.append(picos.Constant(state*state.H))
    return States

def SquareRootMeasurementSuccessPorbability(SquareRoot):
    SquareRootDiagonal = np.diagonal(SquareRoot)
    sumSquare          = 0
    for iElement in SquareRootDiagonal: sumSquare += abs(iElement)**2
    return float(sumSquare)


