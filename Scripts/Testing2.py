########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     29/01/2025                                                                                #
#  Purpose:  Testing                                                                                   #
#                                                                                                      #
########################################################################################################

import picos
import numpy as np
from scipy.linalg import expm

groupSize = 7

sigmaX = picos.Constant([[0,1],[1,0]])

initialState = picos.Constant([1/np.sqrt(2),1/np.sqrt(2)])

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

def setOfGropGeneratedStates(seedState, groupSize, involutionalMatrix):
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
    def groupGenerateStateDensityMatrix(seedState, groupSize, elementNumber):
        if not isinstance(elementNumber, int): raise TyperError("The given value of the elementNumber must be an integer.")
        if elementNumber%groupSize == 0: return seedState*seedState.H
        state = picos.Constant(expm(1j*sigmaX*np.pi*(elementNumber%groupSize)/groupSize))*initialState
        if not np.allclose(state.value, seedState.value): return state*state.H
        raise ValueError(f"The given seed state is an eigen state of the given matrix.")
    return [groupGenerateStateDensityMatrix(seedState,groupSize,iElement) for iElement in range(groupSize)]

setOfStates = setOfGropGeneratedStates(initialState,groupSize,sigmaX)

for i, iState in enumerate(setOfStates):
    print(f"State {i}:")
    print(iState)

