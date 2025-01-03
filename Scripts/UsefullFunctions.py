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

def Hermitian(matrix, SilentMode = True):
    if not isinstance(matrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
        raise TypeError(f"The given matrix is not valid for computing the hermicity")
    if np.allclose(matrix.value, matrix.H.value):
        if not SilentMode : print("The matrix is Hermitian")
        return True
    if not SilentMode : print("The matrix is NOT hermitian")
    return False

def RandomSetOfDensityMatrices (NumberOfMatrices, MatrixDimension):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    def RandomDensityMatrix(Dimension):
        NonHemitianMatix = picos.Constant([[np.random.rand()+1j*np.random.rand() for jDimension in range(Dimension)] 
                                            for iDimension in range(Dimension)])
        HermitianMatix   = NonHemitianMatix.H*NonHemitianMatix
        return HermitianMatix/picos.trace(HermitianMatix)
    return [RandomDensityMatrix(MatrixDimension) for iDensityMatrix in range(NumberOfMatrices)]

def RandomSetOfProbabilities (NumberOfMatrices):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    priorProbabilities = np.random.rand(NumberOfMatrices)
    priorProbabilities /= sum(priorProbabilities)
    return priorProbabilities

def GrupGeneratedDensityMatrices (NumberOfMatrices, MatrixDimension, Overlap):
    if not isinstance(MatrixDimension,int):
        raise TypeError("The given MatrixDimension must be an Integer")
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    def GroupGenerateDensityMatrix(Dimension, Overlap):
        return picos.Constant([[1 if iDimension == jDimension else (Overlap if iDimension < jDimension else Overlap.conjugate()) for jDimension in range(MatrixDimension)]for iDimension in range(MatrixDimension)])
    return [GroupGenerateDensityMatrix(MatrixDimension, Overlap) for iDensityMatrix in range(NumberOfMatrices)]

def SetOfEqualProbabilities (NumberOfMatrices):
    if not isinstance(NumberOfMatrices,int):
        raise TypeError("The given NumberOfMatices must be an Integer")
    return np.array([1/NumberOfMatrices for iMatrix in range (NumberOfMatrices)])
