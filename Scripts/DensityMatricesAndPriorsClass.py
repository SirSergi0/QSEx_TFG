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
    def __init__(self, NumberOfMatrices, MatrixDimension):
        if not isinstance(MatrixDimension,int):
            raise TypeError("The given MatrixDimension must be an Integer")
        if not isinstance(NumberOfMatrices,int):
            raise TypeError("The given NumberOfMatrices must be an Integer")
        self.NumberOfMatrices   = NumberOfMatrices
        self.MatrixDimension    = MatrixDimension
    def to_dict(self):
        return {
            'NumberOfMatrices'   : self.NumberOfMatrices,
            'MatrixDimension'    : self.MatrixDimension,
            'PriorProbabilities' : UsefullFunctions.RandomSetOfProbabilities(self.NumberOfMatrices),
            'DensityMatrices'    : UsefullFunctions.RandomSetOfDensityMatrices(self.NumberOfMatrices,self.MatrixDimension)
        }
    def __getitem__(self, key):
        return self.to_dict()[key]

    def getDesityMatrices(self):
        return self.to_dict()['DensityMatrices']

    def getPriorProbabilities(self):
        return self.to_dict()['PriorProbabilities']

    def getNumberOfMatrices(self):
        return self.to_dict()['NumberOfMatrices']   

    def getMatrixDimension(self):
        return self.to_dict()['MatrixDimension']   

    def __repr__(self):
        PrintingText = "The working paramaters are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText

