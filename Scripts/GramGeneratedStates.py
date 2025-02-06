########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     04/02/2025                                                                                #
#  Purpose:  Defining the GramGeneratedStatesClass                                                     #
#                                                                                                      #
########################################################################################################

import GramGeneratedFunctions
from scipy.linalg import sqrtm
import picos

class GramGeneratedStatesClass:
    def __init__(self, NumberOfMatrices, MatrixDimension, MatrixConstructionMethod, ProbabilitiesContructionMethod = "Equal", Overlap = 0):
        if not isinstance(MatrixDimension,int):
            raise TypeError("The given MatrixDimension must be an Integer")
        if not isinstance(NumberOfMatrices,int):
            raise TypeError("The given NumberOfMatrices must be an Integer")
        if not isinstance(MatrixConstructionMethod,str):
            raise TypeError("The given MatrixConstructionMethod is not a String")
        if not isinstance(ProbabilitiesContructionMethod,str):
            raise TypeError("The given ProbabilitiesContructionMethod is not a String")
        self.NumberOfMatrices               = NumberOfMatrices
        self.MatrixDimension                = MatrixDimension
        self.MatrixConstructionMethod       = MatrixConstructionMethod
        self.ProbabilitiesContructionMethod = ProbabilitiesContructionMethod
        self.PriorProbabilities             = GramGeneratedFunctions.SetOfProbabilities(ProbabilitiesContructionMethod,NumberOfMatrices)
        self.GramMatrix                     = GramGeneratedFunctions.GenerateGramMatrix(NumberOfMatrices,MatrixConstructionMethod, Overlap)
        self.SquareRoot                     = picos.Constant(sqrtm(self.GramMatrix))
        self.GramMatrixWithPriors           = GramGeneratedFunctions.GenerateGramMatrixWithPriors(self.GramMatrix, self.PriorProbabilities,self.NumberOfMatrices)
        self.SquareRootWithPriors           =picos.Constant(sqrtm(self.GramMatrixWithPriors))
        self.DensityMatrices                = GramGeneratedFunctions.GetGramDensityMatrices(self.GramMatrix, self.NumberOfMatrices)
        self.SquareRootSuccessProbability   = GramGeneratedFunctions.SquareRootMeasurementSuccessPorbability(self.SquareRootWithPriors)
    def to_dict(self):
        return {
            'NumberOfMatrices'               : self.NumberOfMatrices,
            'MatrixDimension'                : self.MatrixDimension,
            'PriorProbabilities'             : self.PriorProbabilities,
            'MatrixConstructionMethod'       : self.MatrixConstructionMethod,
            'ProbabilitiesContructionMethod' : self.ProbabilitiesContructionMethod,
            'GramMatrix'                     : self.GramMatrix,
            'SquareRoot'                     : self.SquareRoot,
            'GramMatrixWithPriors'           : self.GramMatrixWithPriors,
            'DensityMatrices'                : self.DensityMatrices,
            'SquareRootSuccessProbability'   : self.SquareRootSuccessProbability,
        }
    def __getitem__(self, key):
        return self.to_dict()[key]

    def getPriorProbabilities(self):
        return self.to_dict()['PriorProbabilities']

    def getNumberOfMatrices(self):
        return self.to_dict()['NumberOfMatrices']   

    def getMatrixDimension(self):
        return self.to_dict()['MatrixDimension']   

    def getDensityMatrices(self):
        return self.to_dict()['DensityMatrices']   

    def getSRMSuccessProbability(self):
        return self.to_dict()['SquareRootSuccessProbability']

    def getGramMatrix(self):
        return self.to_dict()['GramMatrix']
    
    def __repr__(self):
        PrintingText = "The working paramaters are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText

