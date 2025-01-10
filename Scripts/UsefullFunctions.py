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

def SetOfRandomDensityMatrices (NumberOfMatrices, MatrixDimension):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    def RandomDensityMatrix(Dimension):
        return picos.Constant(qiskit.quantum_info.random_density_matrix(Dimension).data)
    return [RandomDensityMatrix(MatrixDimension) for iDensityMatrix in range(NumberOfMatrices)]

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
    if not MatrixMethod in ["Random"]:
        raise ValueError("The given MatrixMethod is not a valid method")
    if MatrixMethod == "Random":
        return SetOfRandomDensityMatrices(NumberOfMatrices,MatrixDimension)
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

