########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     03/01/2025                                                                                #
#  Purpose:  Defining the DensityMaticesAndPriors class                                                #
#                                                                                                      #
########################################################################################################

import UsefullFunctions
import numpy as np 
import picos

class DensityMaticesAndPriors:
    def __init__(self, NumberOfMatrices, MatrixDimension, MatrixConstructionMethod, ProbabilitiesContructionMethod = "Equal", seedState = None, involutoryMatrix = None):
        if not isinstance(MatrixDimension,int):
            raise TypeError("The given MatrixDimension must be an Integer")
        if not isinstance(NumberOfMatrices,int):
            raise TypeError("The given NumberOfMatrices must be an Integer")
        if not isinstance(MatrixConstructionMethod,str):
            raise TypeError("The given MatrixConstructionMethod is not a String")
        if not isinstance(ProbabilitiesContructionMethod,str):
            raise TypeError("The given ProbabilitiesContructionMethod is not a String")
        if MatrixConstructionMethod == "GroupGeneratedStates":
            if not isinstance(seedState,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
                raise TypeError("The given seed state is not valid.")
            if not isinstance(involutoryMatrix,(picos.expressions.exp_affine.AffineExpression,
                              picos.expressions.exp_affine.ComplexAffineExpression)):
                raise TypeError("The given involutory matrix is not valid.")
            if not seedState.shape[1] == 1:
                raise ValueError("The given seed state must be a vector")
            if not seedState.shape[0] == MatrixDimension:
                raise ValueError("The given seed state dimention does not coincide with the MatrixDimension")
            if not (len(involutoryMatrix.shape) == 2 and (involutoryMatrix.shape[0] == involutoryMatrix.shape[1])):
                raise TypeError("The given involutory matrix must be a square matrix.")
            if not involutoryMatrix.shape[0] == MatrixDimension:
                raise ValueError("The given involutory matrix dimension does not coincide with the MatrixDimension")
        self.NumberOfMatrices               = NumberOfMatrices
        self.MatrixDimension                = MatrixDimension
        self.MatrixConstructionMethod       = MatrixConstructionMethod
        self.ProbabilitiesContructionMethod = ProbabilitiesContructionMethod
        self.PriorProbabilities             = UsefullFunctions.SetOfProbabilities(ProbabilitiesContructionMethod,NumberOfMatrices)
        self.DensityMatrices                = UsefullFunctions.SetOfMatrices(MatrixConstructionMethod, NumberOfMatrices, MatrixDimension, seedState, involutoryMatrix)

    def to_dict(self):
        return {
            'NumberOfMatrices'               : self.NumberOfMatrices,
            'MatrixDimension'                : self.MatrixDimension,
            'PriorProbabilities'             : self.PriorProbabilities,
            'DensityMatrices'                : self.DensityMatrices,
            'MatrixConstructionMethod'       : self.MatrixConstructionMethod,
            'ProbabilitiesContructionMethod' : self.ProbabilitiesContructionMethod,
            'GramMatrix'                     : UsefullFunctions.GramMatrix(self.NumberOfMatrices, self.DensityMatrices),
            'GramMatrixWithPriors'           : UsefullFunctions.GramMatrixWithPriors(self.NumberOfMatrices, self.DensityMatrices, self.PriorProbabilities),
            'SquareRootMeasurement'          : UsefullFunctions.SquareRootMeasurement(self.NumberOfMatrices, self.DensityMatrices),
            'SquareRootSuccessProbability'   : UsefullFunctions.SquareRootMeasurementSuccessPorbability(self.NumberOfMatrices, self.DensityMatrices, self.PriorProbabilities)
        }
    def __getitem__(self, key):
        return self.to_dict()[key]

    def getDensityMatrices(self):
        return self.to_dict()['DensityMatrices']

    def getPriorProbabilities(self):
        return self.to_dict()['PriorProbabilities']

    def getNumberOfMatrices(self):
        return self.to_dict()['NumberOfMatrices']   

    def getMatrixDimension(self):
        return self.to_dict()['MatrixDimension']   

    def getGramMatrix(self):
        return self.to_dict()['GramMatrix']

    def getGramMatrixWithPriors(self):
        return self.to_dict()['GramMatrixWithPriors']

    def __repr__(self):
        PrintingText = "The working paramaters are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText

