########################################################################################################
#                                                                                                      #
#  Project:  Physics TFG                                                                               #
#  Author:   Sergio Castañeiras Morales                                                                #
#  Date:     02/01/2025                                                                                #
#  Purpose:  Testing the Square Root Measurement                                                       #
#                                                                                                      #
########################################################################################################
import numpy
import picos

# Lets create an hermitian we start by defining a non hermitian matrix '_P'
_P = picos.Constant([
[ 1 -1j, 2 +2j, 1 ],
[ 3j, -2j, -1 -1j],
[ 1 +2j, -0.5+1j, 1.5 ]])

# now we construct P as _P*_P^dagger which ensures that P is hermitian

P = (_P*_P.H).renamed("P")

# We do the same procedure with a Q matrix

_Q = picos.Constant([
[-1 -2j, 2j, 1.5 ],
[ 1 +2j, -2j, 2.0-3j],
[ 1 +2j, -1 +1j, 1 +4j]])
Q = (_Q*_Q.H).renamed("Q")

# Define the problem.
F = picos.Problem()
Z = picos.ComplexVariable("Z", P.shape)
F.set_objective("max", 0.5*picos.trace(Z + Z.H))
F.add_constraint(((P & Z) // (Z.H & Q)) >> 0)
print(F)
# Solve the problem.
F.solve(solver = "cvxopt")
print("\nOptimal value:", round(F, 4))
print("Optimal Z:", Z.value, sep="\n")
# Also compute the fidelity via NumPy for comparison.
PP = numpy.matrix(P.value)
QQ = numpy.matrix(Q.value)
S,U = numpy.linalg.eig(PP)
sqP = U * numpy.diag([s**0.5 for s in S]) * U.H # Square root of P.
S,U = numpy.linalg.eig(QQ)
sqQ = U * numpy.diag([s**0.5 for s in S]) * U.H # Square root of Q.
Fnp = sum(numpy.linalg.svd(sqP * sqQ)[1]) # Trace-norm of sqrt(P)·sqrt(Q).
print("Fidelity F(P,Q) computed by NumPy:", round(Fnp, 4))
