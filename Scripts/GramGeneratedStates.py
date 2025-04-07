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
import numpy as np

class GramGeneratedStatesClass:
    def __init__(self, NumberOfMatrices, MatrixDimension, MatrixConstructionMethod, GenerationConditions, ProbabilitiesContructionMethod = "Equal"):
        if not isinstance(MatrixDimension,int):
            raise TypeError("The given MatrixDimension must be an Integer")
        if not isinstance(NumberOfMatrices,int):
            raise TypeError("The given NumberOfMatrices must be an Integer")
        if not isinstance(MatrixConstructionMethod,str):
            raise TypeError("The given MatrixConstructionMethod is not a String")
        if not isinstance(ProbabilitiesContructionMethod,str):
            raise TypeError("The given ProbabilitiesContructionMethod is not a String")
        if not isinstance(GenerationConditions, (ZnGramMatrixConditionsOverlap,ZnGramMatrixConditionsEigenValues)):
            raise TypeError("The given GenerationConditions type in not valid.")
        self.NumberOfMatrices               = NumberOfMatrices
        self.MatrixDimension                = MatrixDimension
        self.MatrixConstructionMethod       = MatrixConstructionMethod
        self.ProbabilitiesContructionMethod = ProbabilitiesContructionMethod
        self.PriorProbabilities             = GramGeneratedFunctions.SetOfProbabilities(ProbabilitiesContructionMethod,NumberOfMatrices)
        GramMatrixAndUnitaryMatrix          = GramGeneratedFunctions.GenerateGramMatrix(NumberOfMatrices,MatrixConstructionMethod, GenerationConditions)
        self.GramMatrix                     = GramMatrixAndUnitaryMatrix["GramMatrix"]
        self.SquareRoot                     = picos.Constant(sqrtm(self.GramMatrix))
        self.GramMatrixWithPriors           = GramGeneratedFunctions.GenerateGramMatrixWithPriors(self.GramMatrix, self.PriorProbabilities,self.NumberOfMatrices)
        self.SquareRootWithPriors           = picos.Constant(sqrtm(self.GramMatrixWithPriors))
        self.DensityMatrices                = GramGeneratedFunctions.GetGramDensityMatrices(self.GramMatrix, self.NumberOfMatrices)
        self.SquareRootSuccessProbability   = GramGeneratedFunctions.SquareRootMeasurementSuccessPorbability(self.SquareRootWithPriors)
        self.PerfectExlusion                = GramGeneratedFunctions.PerfectExlusion(self.GramMatrixWithPriors)
        self.UnitaryMatrix                  = GramMatrixAndUnitaryMatrix["UnitaryMatrix"]
    def to_dict(self):
        return {
            'NumberOfMatrices'               : self.NumberOfMatrices,
            'MatrixDimension'                : self.MatrixDimension,
            'PriorProbabilities'             : self.PriorProbabilities,
            'MatrixConstructionMethod'       : self.MatrixConstructionMethod,
            'ProbabilitiesContructionMethod' : self.ProbabilitiesContructionMethod,
            'GramMatrix'                     : self.GramMatrix,
            'GramMatrixWithPriors'           : self.GramMatrixWithPriors,
            'SquareRoot'                     : self.SquareRoot,
            'SquareRootWithPriors'           : self.SquareRootWithPriors,
            'DensityMatrices'                : self.DensityMatrices,
            'SquareRootSuccessProbability'   : self.SquareRootSuccessProbability,
            'PerfectExlusion'                : self.PerfectExlusion,
            'UnitaryMatrix'                  : self.UnitaryMatrix
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

    def getGramMatrixWithPriors(self):
        return self.to_dict()['GramMatrixWithPriors']
    
    def getSquareRootMatrix(self):
        return self.to_dict()['SquareRoot']

    def getPerfectExlusion(self):
        return self.to_dict()['PerfectExlusion']

    def getPerfectExlusionLowerBoundMinimumError(self):
        if self.getPerfectExlusion(): return 1 
        return GramGeneratedFunctions.PerfectExlusionLowerBoundMinimumError(self)

    def getPerfectExlusionLowerBoundZeroError(self):
        if self.getPerfectExlusion(): return 1 
        return GramGeneratedFunctions.PerfectExlusionLowerBoundZeroError(self)

    def getOverlapsAndPhases(self):
        GramMatrix   = self.getGramMatrix()
        OverlapsList = []
        PhaseList    = []
        for iElement in GramMatrix:
            Overlap = round(np.abs(iElement.value),10)
            Phase = round(abs(np.angle(iElement.value)),10)
            if not Overlap in OverlapsList : OverlapsList.append(Overlap)
            if not Phase in PhaseList : PhaseList.append(Phase)
        return OverlapsList, PhaseList
    
    def getUnitaryMatrix(self):
        return self.to_dict()['UnitaryMatrix']

    def __repr__(self):
        PrintingText = "The working paramaters are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText

class ZnGramMatrixConditionsOverlap:
    def __init__(self, OverlapsList, PhaseList, n):
        if not isinstance(n,int):
            raise TypeError("The given n value must be an Integer.")
        if not isinstance(OverlapsList, list):
            raise TypeError("The given OverlapList must be a list.")
        if not isinstance(PhaseList, list):
            raise TypeError("The given PhaseList must be a list.")
        for iElement in range(len(OverlapsList)):
            if not isinstance(OverlapsList[iElement], (int,float)):
                raise TypeError(f"The element {OverlapsList[iElement]} of the OverlapsList must be an integer or a float")
            if (OverlapsList[iElement] > 1) or (OverlapsList[iElement] < 0):
                raise ValueError(f"The element {OverlapsList[iElement]} of the OverlapsList must be a number between 0 and 1.")
        for iElement in range(len(PhaseList)):
            if not isinstance(PhaseList[iElement], (int,float)):
                raise TypeError(f"The element {PhaseList[iElement]} of the PhaseList must be an integer or a float")
        if not len(OverlapsList) == int(np.floor(n/2)):
            raise ValueError(f"The given OverlapsList must be a list of {int(np.floor(n/2))} element(s).")
        if not len(PhaseList) == (int(np.floor((n-1)/2))):
            raise ValueError(f"The given PhaseList must be a list of {int(np.floor((n-1)/2))} element(s).")

        self.OverlapsList = OverlapsList
        self.PhaseList    = PhaseList
        self.n            = n
    
    def to_dict(self):
        return {
            'OverlapsList'  : self.OverlapsList,
            'PhaseList'     : self.PhaseList,
            'n'             : self.n
        }

    def getOverlapList(self):
        return self.to_dict()['OverlapsList']

    def getPhaseList(self):
        return self.to_dict()['PhaseList']

    def getN(self):
        return self.to_dict()['n']
    
    def __repr__(self):
        PrintingText = "The Zn Gram Matrix Conditions are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText
        
class ZnGramMatrixConditionsEigenValues:
    def __init__(self, EigenValues, n):
        if not isinstance(n,int):
            raise TypeError(f"The n type is must be an integer.")
        if not isinstance(EigenValues,list):
            raise TypeError(f"The EigenValues type is must be a list.")
        if not len(EigenValues) == n:
            raise ValueError(f"The EigenValues list lenght {len(EigenValues)} must be equal to the 'n'={n} value.")
        for iElement in range(n):
            if not isinstance(EigenValues[iElement],(int,float)):
                raise TypeError(f"The {iElement}-th element of the EigenValues must be eighter an integer or a float.")
            if iElement < 0:
                raise ValueError(f"The gram matrix eigenvalues must be positive.")
        if not round(sum(EigenValues),14) == float(n):
            raise ValueError(f"The sum of the eigenvalues must add up to {n}.")
        self.EigenValues  = EigenValues
        self.n            = n
    
    def to_dict(self):
        return {
            'EigenValues'  : self.EigenValues,
            'n'            : self.n
        }

    def getEigenValues(self):
        return self.to_dict()['EigenValues']

    def getN(self):
        return self.to_dict()['n']
    
    def __repr__(self):
        PrintingText = "The Zn Gram Matrix Conditions are:\n"
        for iVariable in self.to_dict():
            PrintingText+= f'{iVariable :<30}:{self.to_dict()[iVariable] }\n'
        return PrintingText
        
