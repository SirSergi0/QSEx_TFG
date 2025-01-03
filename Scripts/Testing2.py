import QSExSetUp
import picos as pc
import numpy as np

NumberOfMatrices       = 10
MatrixDimension        = 5
RandomConditions       = QSExSetUp.DensityMatricesAndPriorsClass.DensityMaticesAndPriors(NumberOfMatrices,MatrixDimension)
MyDenisityMatrices     = RandomConditions.getDesityMatrices()
MyPriorsPropbabilities = RandomConditions.getPriorProbabilities()
# Assuming RandomConditions is a class you have already defined, which provides the matrices
# Replace this with your actual data:
MyDenisityMatrices = RandomConditions.getDesityMatrices()  # Replace with actual method to get matrices
MyPriorsProbabilities = RandomConditions.getPriorProbabilities()  # Replace with actual method to get priors

# Define the problem
F = pc.Problem()

# Define Z as a complex decision variable (Hermitian matrix)
Z = F.add_variable("Z", MyDenisityMatrices[0].shape, "hermitian")

# Set the objective: maximize trace(Z * MyDenisityMatrices[0])
F.set_objective("max", pc.trace(Z * MyDenisityMatrices[0]))

# Add constraint: trace(Z) == 1
F.add_constraint(pc.trace(Z) == 1)

# Add positive semidefinite constraint: Z must be positive semidefinite (Z >> 0)
F.add_constraint(Z >> 0)

# Solve the problem using cvxopt solver
F.solve(solver="cvxopt")

# Output the results
print("\nOptimal value:", round(F.value, 4))
print("Optimal Z:", Z.value, sep="\n")
