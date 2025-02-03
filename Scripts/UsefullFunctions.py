########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Store rellevant functions and utilities                                                   #
#                                                                                                      #
########################################################################################################

import numpy as np 
import picos
import qiskit
from scipy.linalg import expm

def Hermitian(matrix, SilentMode = True):
    if not isinstance(matrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given matrix is not valid for computing the hermicity")
    if np.allclose(matrix.value, matrix.H.value):
        if not SilentMode : print("The matrix is Hermitian")
        return True
    if not SilentMode : print("The matrix is NOT hermitian")
    return False

def Involutional(matrix, SilentMode = True):
    if not isinstance(matrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given matrix is not valid")
    if not matrix.shape[0] == matrix.shape[1]: 
        raise ValueError(f"The given matrix must be a square matrix")
    if np.allclose(matrix*matrix.value,np.eye(matrix.shape[0])):
        if not SilentMode : print("The matrix is Involutional")
        return True
    if not SilentMode : print("The matrix is NOT Involutional")
    return False

def SetOfRandomDensityMatricesMixed(NumberOfMatrices, MatrixDimension):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    def RandomDensityMatrixMixed(Dimension):
        return picos.Constant(qiskit.quantum_info.random_density_matrix(Dimension).data)
    return [RandomDensityMatrixMixed(MatrixDimension) for iDensityMatrix in range(NumberOfMatrices)]

def SetOfRandomDensityMatricesPure (NumberOfMatrices, MatrixDimension):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    def RandomDensityMatrixPure(Dimension):
        randomState         = qiskit.quantum_info.random_statevector(Dimension)
        randomDensityMatrix = np.outer(randomState.data, np.conj(randomState.data))
        return picos.Constant("randomDensityMatrix",(randomDensityMatrix + randomDensityMatrix.T.conj())/2)
    return [RandomDensityMatrixPure(MatrixDimension) for iDensityMatrix in range(NumberOfMatrices)]

def setOfGropGeneratedDensityMatrices(seedState, groupSize, involutionalMatrix):
    if not isinstance(groupSize, int): raise TypeError("The given value of the groupSize must be an integer.")
    if not isinstance(seedState,(picos.expressions.exp_affine.AffineExpression,
                                 picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given seed state is not valid state type")
    if not isinstance(involutionalMatrix,(picos.expressions.exp_affine.AffineExpression,
                                 picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given involutional matrix is not a valid matrix type")
    if not involutionalMatrix.shape[0]==seedState.shape[0]:
        raise ValueError(f"The given seed state and seed matrix are not of the same size")
    if not Involutional(involutionalMatrix): 
        raise ValueError(f"The given matrix must be an involutional matrix (i.e) its square must be the identity.")
    if not round(float(abs(seedState)**2),7)==1:
        raise ValueError(f"The given seed state mus be normalized.")
    def groupGenerateStateDensityMatrix(seedState, groupSize, elementNumber, involutionalMatrix):
        if not isinstance(elementNumber, int): raise TyperError("The given value of the elementNumber must be an integer.")
        if elementNumber%groupSize == 0: return picos.Constant(seedState*seedState.H)
        state = picos.Constant(expm(1j*involutionalMatrix*np.pi*(elementNumber%groupSize)/groupSize))*seedState
        if not np.allclose(state.value, seedState.value): return picos.Constant(state*state.H)
        raise ValueError(f"The given seed state is an eigen state of the given matrix.")
    return [groupGenerateStateDensityMatrix(seedState,groupSize,iElement, involutionalMatrix) for iElement in range(groupSize)]

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

def GramMatrix(NumberOfMatrices, MyDensityMatrices):
    def Overlap (braDenistyMatrix, ketDensityMatrix):
        braEigenValues, braEigenVectors = np.linalg.eigh(braDenistyMatrix)
        ketEigenValues, ketEigenVectors = np.linalg.eigh(ketDensityMatrix)
        return np.vdot(braEigenVectors[:, np.argmax(braEigenValues)], ketEigenVectors[:, np.argmax(ketEigenValues)])
    GramMatrix = []
    for iMatrix in range(NumberOfMatrices):
        GramRow = []
        for jMatrix in range(NumberOfMatrices):
            GramRow.append(Overlap(MyDensityMatrices[iMatrix],MyDensityMatrices[jMatrix]))
        GramMatrix.append(GramRow)
    return picos.Constant(GramMatrix)

def GramMatrixWithPriors(NumberOfMatrices, MyDensityMatrices, priorProbabilities):
    def Overlap (braDenistyMatrix, ketDensityMatrix):
        braEigenValues, braEigenVectors = np.linalg.eigh(braDenistyMatrix)
        ketEigenValues, ketEigenVectors = np.linalg.eigh(ketDensityMatrix)
        return np.vdot(braEigenVectors[:, np.argmax(braEigenValues)], ketEigenVectors[:, np.argmax(ketEigenValues)])
    GramMatrix = []
    for iMatrix in range(NumberOfMatrices):
        GramRow = []
        for jMatrix in range(NumberOfMatrices):
            GramRow.append(Overlap(MyDensityMatrices[iMatrix],MyDensityMatrices[jMatrix])*np.sqrt(priorProbabilities[iMatrix])*np.sqrt(priorProbabilities[jMatrix]))
        GramMatrix.append(GramRow)
    return picos.Constant(GramMatrix)
