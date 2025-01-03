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
    def __init__(self, NumberOfMatrices, MatrixDimension, ConstructionMethod, Overlap = None):
        if not isinstance(MatrixDimension,int):
            raise TypeError("The given MatrixDimension must be an Integer")
        if not isinstance(NumberOfMatrices,int):
            raise TypeError("The given NumberOfMatrices must be an Integer")
        if not isinstance(ConstructionMethod,str):
            raise TypeError("The given ConstructionMethod is not a String")
        if not (ConstructionMethod in ['Random', 'GroupGenerated']):
            raise ValueError("The given ConstructionMethod must be either: ", ['Random','GroupGenerated'])
        if (ConstructionMethod != 'Random') and (Overlap == None):
            raise ValueError("A Overlap value must be given.")
            if (not isinstance(Overlap,complex)) or (not abs(Overlap)<1):
                raise ValueError("The given Overlap be a complex number with modulus less than 1.")
        self.NumberOfMatrices   = NumberOfMatrices
        self.MatrixDimension    = MatrixDimension
        self.ConstructionMethod = ConstructionMethod
        self.Overlap            = Overlap
    def to_dict(self):
        if self.ConstructionMethod == "Random":
            return {
                'NumberOfMatrices'   : self.NumberOfMatrices,
                'MatrixDimension'    : self.MatrixDimension,
                'PriorProbabilities' : UsefullFunctions.RandomSetOfProbabilities(self.NumberOfMatrices),
                'DensityMatrices'    : UsefullFunctions.RandomSetOfDensityMatrices(self.NumberOfMatrices,self.MatrixDimension),
                'Overlap'            : None,
                'ConstructionMethod' : "Random"
            }
        if self.ConstructionMethod == "GroupGenerated":
            return {
                'NumberOfMatrices'   : self.NumberOfMatrices,
                'MatrixDimension'    : self.MatrixDimension,
                'PriorProbabilities' : UsefullFunctions.RandomSetOfProbabilities(self.NumberOfMatrices),
                'DensityMatrices'    : UsefullFunctions.GrupGeneratedDensityMatrices(self.NumberOfMatrices,self.MatrixDimension, self.Overlap),
                'Overlap'            : self.Overlap,
                'ConstructionMethod' : "GroupGenerated"
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

    def getOverlap(self):
        return self.to_dict()['Overlap']

    def __repr__(self):
        PrintingText = "The working paramaters are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText

