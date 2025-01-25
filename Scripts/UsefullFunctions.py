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

def Hermitian(matrix, SilentMode = True):
    if not isinstance(matrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given matrix is not valid for computing the hermicity")
    if np.allclose(matrix.value, matrix.H.value):
        if not SilentMode : print("The matrix is Hermitian")
        return True
    if not SilentMode : print("The matrix is NOT hermitian")
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

def SetOfMatrices (MatrixMethod, NumberOfMatrices, MatrixDimension):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    if not isinstance(MatrixMethod, str):
        raise TypeError("The given MatrixMethod must be a string")
    if not MatrixMethod in ["RandomMixedStates", "RandomPureStates"]:
        raise ValueError("The given MatrixMethod is not a valid method")
    if MatrixMethod == "RandomMixedStates":
        return SetOfRandomDensityMatricesMixed(NumberOfMatrices,MatrixDimension)
    if MatrixMethod == "RandomPureStates":
        return SetOfRandomDensityMatricesPure(NumberOfMatrices, MatrixDimension)
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
