import simplex
import numpy

## MAIN PROGRAM
splx = simplex.Simplex()
v = splx.simplex( 1.5,
    1.3,
    0.8,
    [numpy.array((-1.2,1)) , numpy.array((-3,5)), numpy.array((1.57,2.6))], 
    0.00005,
    splx.F_Rosenbrock)
print("::::::::::::::FINAL RESULTS::::::::::::::")
print("P =",v[0],", computed in %i iterations" % (v[1]))