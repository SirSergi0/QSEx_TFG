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
import GramGeneratedStates

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

def GenerateGramMatrix(NumberOfMatrices, MatrixMethod, GenerationConditions):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatrices must be an Integer")
    ValidMethods = ["ZnOverlap", "ZnEigenValues"]
    if not MatrixMethod in ValidMethods:
        raise ValueError(f"The given MatrixMethod is not a valid method, the valid methods are: {ValidMethods}")
    if not isinstance(GenerationConditions, (GramGeneratedStates.ZnGramMatrixConditionsOverlap,GramGeneratedStates.ZnGramMatrixConditionsEigenValues)):
        raise TypeError("The GenerationConditions type is not valid")
    if MatrixMethod == "ZnOverlap":
        return ZnGramMatrixOverlap(NumberOfMatrices, GenerationConditions)
    if MatrixMethod == "ZnEigenValues":
        return ZnGramMatrixEigenValues(GenerationConditions)

def GenerateGramMatrixWithPriors(GramMatrix, priorProbabilities, NumberOfStates):
    GramMatrixWithPriors = np.array(GramMatrix)
    for iState in range(NumberOfStates):
        for jState in range(NumberOfStates):
            GramMatrixWithPriors[iState][jState] *= (np.sqrt(priorProbabilities[iState])*np.sqrt(priorProbabilities[jState]))
    return picos.Constant(GramMatrixWithPriors)

def ZnGramMatrixOverlap(NumberOfStates, GenerationConditions):
    if not isinstance(GenerationConditions, GramGeneratedStates.ZnGramMatrixConditionsOverlap):
        raise TypeError("The GenerationConditions type is not valid")
    GramMatrix  = []
    OverlapList = GenerationConditions.getOverlapList()
    PhaseList   = GenerationConditions.getPhaseList()
    n           = GenerationConditions.getN()
    def ListToMatrixZn (iElement,jElement,n):
        if abs(iElement-jElement)<n/2: return int(abs(iElement-jElement)-1)
        return int(n - 1 - abs(iElement-jElement))
    for iState in range(NumberOfStates):
        GramMatrixRow = []
        for jState in range(NumberOfStates):
            if iState == jState:
                GramMatrixRow.append(1)
                continue
            if abs(iState-jState) == n/2:
                GramMatrixRow.append(OverlapList[ListToMatrixZn(iState,jState,n)])
                continue
            if ((iState > jState) and (abs(iState-jState) < n/2)) or ((iState < jState) and (abs(iState-jState) > n/2)):
                GramMatrixRow.append(OverlapList[ListToMatrixZn(iState,jState,n)]*np.exp(-1j * PhaseList[ListToMatrixZn(iState,jState,n)]))
                continue
            GramMatrixRow.append(OverlapList[ListToMatrixZn(iState,jState,n)]*np.exp(1j * PhaseList[ListToMatrixZn(iState,jState,n)]))
        GramMatrix.append(GramMatrixRow)
    return 0.5*(picos.Constant(GramMatrix)+picos.Constant(GramMatrix).H)

def ZnGramMatrixEigenValues (GenerationConditions):
    if not isinstance(GenerationConditions, GramGeneratedStates.ZnGramMatrixConditionsEigenValues):
        print(type(GenerationConditions))
        raise TypeError("The GenerationConditions type is not valid")
    def fourierBasisMatrix(n):
        def foureirVector(baseElement, n):
            return [np.exp((2j*np.pi*iElement*baseElement)/(n))/(np.sqrt(n)) for iElement in range (n)]
        fourierMatrixList = []
        for iBaseElement in range(1,n+1):
            fourierMatrixList.append(foureirVector(iBaseElement,n))
        fourierMatrix = picos.Constant(fourierMatrixList)
        return fourierMatrix.T
    def diagonalMatrix(diagonal,n):
        diagonalMatrixList = []
        for iElement in range(n):
            diagonalMatrixRow = []
            for jElement in range(n):
                if jElement == iElement:
                    diagonalMatrixRow.append(diagonal[jElement])
                    continue
                diagonalMatrixRow.append(0)
            diagonalMatrixList.append(diagonalMatrixRow)
        return picos.Constant(diagonalMatrixList)
    EigenValues  = GenerationConditions.getEigenValues()
    n            = GenerationConditions.getN()
    fourierBasis = fourierBasisMatrix(n)
    ZnGramMatrix = picos.Constant(fourierBasis*diagonalMatrix(EigenValues,n)*fourierBasis.H)
    return 0.5*(ZnGramMatrix+ZnGramMatrix.H)

def GetGramDensityMatrices(GramMatrix, NumberOfStates):
    SquareRoot = np.array(sqrtm(GramMatrix))
    States     = []
    for iState in range(NumberOfStates):
        state  = picos.Constant(SquareRoot[:,iState])
        state /= abs(state)
        States.append(picos.Constant(state*state.H))
    return States

def SquareRootMeasurementSuccessPorbability(SquareRoot):
    SquareRootDiagonal = np.diagonal(SquareRoot)
    sumSquare          = 0
    for iElement in SquareRootDiagonal: sumSquare += abs(iElement)**2
    return float(sumSquare)

def PerfectExlusion(GramMatrix):
    GramEigenValuesSqrt = []
    for iEigenValue in (np.linalg.eigvals(np.array(GramMatrix.value))):
        GramEigenValuesSqrt.append(float(np.sqrt(abs(iEigenValue))))
    GramEigenValuesSqrt = sorted(GramEigenValuesSqrt, reverse = True)
    if GramEigenValuesSqrt[0] <= sum(GramEigenValuesSqrt[1:]) : return True
    return False

def PerfectExlusionLowerBoundMinimumError(Conditions):
    GramMatrix          = Conditions.getGramMatrixWithPriors()
    GramEigenValuesSqrt = []
    for iEigenValue in (np.linalg.eigvals(np.array(GramMatrix.value))):
        GramEigenValuesSqrt.append(float(np.sqrt(abs(iEigenValue))))
    GramEigenValuesSqrt = sorted(GramEigenValuesSqrt, reverse = True)
    if GramEigenValuesSqrt[0] <= sum(GramEigenValuesSqrt[1:]) : return 1
    return 1-((GramEigenValuesSqrt[0] - sum(GramEigenValuesSqrt[1:]))/(Conditions.getNumberOfMatrices()))**2

def PerfectExlusionLowerBoundZeroError(Conditions):
    GramMatrix           = Conditions.getGramMatrixWithPriors()
    GramEigenValuesSqrt  = []
    for iEigenValue in (np.linalg.eigvals(np.array(GramMatrix.value))):
        GramEigenValuesSqrt.append(float(np.sqrt(abs(iEigenValue))))
    GramEigenValuesSqrt = sorted(GramEigenValuesSqrt, reverse = True)
    if GramEigenValuesSqrt[0] <= sum(GramEigenValuesSqrt[1:]) : return 1
    return 1-(sum(GramEigenValuesSqrt)*(GramEigenValuesSqrt[0] - sum(GramEigenValuesSqrt[1:]))/(Conditions.getNumberOfMatrices()))

