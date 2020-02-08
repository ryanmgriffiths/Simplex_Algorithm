# Simplex_Algorithm
Python implementation of a Nelder-Mead simplex algorithm adapted from the original paper [1] using numpy functions. Simplex method is contained within the ```simplex.py``` python module.

```python 
import simplex

# Arbitrary function to minimise, input is a 
# vertex of the simplex as a numpy array, e.g. numpy.array((10,20))
def F(vertex):
    return vertex[0] ** 2 - 5 * vertex[1]

# Simplex class object
splx = simplex.Simplex()

# Calculate the minimum of function F with the simplex method
# ALPHA, GAMMA, BETA are the reflection, expansion, contraction coefficients
# THRESHOLD is the minimisation check for the smallest vertex of the simplex.
# F is the function to minimise as described above.
minimum = splx.simplex(ALPHA, GAMMA, BETA, THRESHOLD, F)
print(minimum)
```

```python
# Output
>> numpy.array((MIN_VECTOR)), number_of_iterations
```
simplex_test.py contains an example in which we minimise the Rosenbrock function with a = 1, b = 100 and is able to converge on the minimum at (1,1).

[1] J. A. Nelder, R. Mead, A Simplex Method for Function Minimization, The Computer Journal, Volume 7, Issue 4, January 1965, Pages 308–313, https://doi.org/10.1093/comjnl/7.4.308
