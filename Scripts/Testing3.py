import QSExSetUp
import numpy as np
import picos
from tqdm import tqdm
from scipy.stats import unitary_group

NumberOfMatrices              = 5
MatrixDimension               = 5
MatrixGenerationMethod        = "Random"
ProbabliltiesGenerationMethod = "Equal"
EigenValues                   = [0.1,0.1,0.1,0.1,100]
EigenValuesNormalized         = [iEigenValue*NumberOfMatrices/sum(EigenValues) for iEigenValue in EigenValues]
Iterations                    = 10000
GroupGeneratedConditions      = QSExSetUp.GramGeneratedStates.GramGeneratedStatesClass(NumberOfMatrices,MatrixDimension, "ZnEigenValues", QSExSetUp.GramGeneratedStates.ZnGramMatrixConditionsEigenValues(EigenValuesNormalized, NumberOfMatrices), ProbabilitiesContructionMethod = ProbabliltiesGenerationMethod)

print(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - GroupGeneratedConditions.getUnitaryMatrix(),1))

DistanceFrobenius             = []
Distance1                     = []
DistanceInfinite              = []
DistanceSpectral              = []
DistanceNuclear               = []
ExlusionSuccessProbability    = []

def random_unitary(n):
    # Create a random complex matrix
    A = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    # QR decomposition
    Q, R = np.linalg.qr(A)
    # Normalize to make Q unitary
    D = np.diag(R) / np.abs(np.diag(R))
    Q = Q @ np.diag(D.conj())
    return Q
from scipy.linalg import expm

def skew_hermitian(n):
    A = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    return A - A.conj().T

def unitary_from_exponential(n):
    A = skew_hermitian(n)
    return expm(A)

for iIteration in tqdm(range(Iterations),"Computing SDP"):
    DistanceFrobenius.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - random_unitary(NumberOfMatrices),2))
    Distance1.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - unitary_group.rvs(NumberOfMatrices),2))
    DistanceNuclear.append(picos.Norm(GroupGeneratedConditions.getUnitaryMatrix() - unitary_from_exponential(NumberOfMatrices),2))


print(min(DistanceFrobenius))
print(min(Distance1))
print(min(DistanceNuclear))
