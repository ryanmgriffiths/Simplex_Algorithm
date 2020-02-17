import numpy

class Simplex:
    """ Class containing functions used to construct a Nelder-Mead simplex algorithm, also contains a method to minimise a function of n free parameters with the simplex algorithm. Requires numpy to function."""

    def __init__(self):
        return None

    def F_Rosenbrock(self,input):
        '''Rosenbrock function for testing, input is arry type (x,y)'''
        return (1-input[0])**2 + 100*(input[1]-input[0]**2)**2

    def centroid(self,input):
        """Calculate the centroid of list of simplex points \n
        Simplex.centroid(input) \n
        Input is a python list of numpy arrays of the simplex points [P(0)...P(n+1)] where P(0) is an n dimensional vector reoresenting a vertex of the simplex. \n
        Returns the centroid vector as a numpy array. \n"""
        centroid = numpy.zeros(len(input[0]))
        idx = 0

        for x in input:
            centroid = centroid+x
            idx = idx + 1 
        
        centroid = centroid/idx
        return centroid

    def reflection(self, pbar, ph, alpha):
        """Returns simplex reflection coordinate: \n
        Simplex.reflection(centroid, ph, alpha) \n
        Centroid is a numpy array of centroid vector, ph is the largest simplex point, alpha is the reflection coefficient. \n
        """
        return (1+alpha)*pbar - alpha*ph

    def expansion(self, pbar, pstar, gamma):
        """Returns simplex expansion coordinate: \n
        Simplex.reflection(centroid, p*, gamma) \n
        Centroid is a numpy array of centroid vector, p* see Nelder-Mead paper, gamma is the expansion coefficient. \n
        """
        return gamma*pstar + (1-gamma)*pbar

    def contraction(self, pbar, ph, beta):
        """Returns simplex reflection coordinate: \n
        Simplex.reflection(centroid, ph, beta) \n
        Centroid is a numpy array of centroid vector, ph is the largest simplex point, beta is the contraction coefficient. \n"""
        return beta*ph + (1-beta)*pbar

    def simplex(self, ALPHA, GAMMA, BETA, P0, THRESHOLD, f):
        """ Use the Nelder-Mead simplex algorithm to minimise a function with n free parameters.\n
        Simplex.simplex(ALPHA, GAMMA, BETA, P0, THRESHOLD)\n
        ALPHA, GAMMA, BETA are the reflection, expansion and contraction coefficients, P0 are our initial simplex vertex coordinates\n
        THRESHOLD is the end condition for our simplex target, and F is the function to minimise of the form F([P(0), ...P(n+1)]), i.e a single argument\n
        where P(i) is an n dimensional vector representing a vertex of the simplex, F must also return a single value.\n
        Returns a tuple of (pi, iterations) where pi is the minimised coordinate and iterations are the number of process iterations to reach this.\n
        """
        # Assert inputs
        assert(ALPHA > 0)
        assert(GAMMA > 1)
        assert(BETA > 0 and BETA < 1)
        assert(len(P0) == len(P0[0])+1)
        
        pi = P0
        yi = numpy.array([f(p) for p in pi])
        yh = numpy.amax(yi)
        yl = numpy.amin(yi)
        ph = pi[int(numpy.where(yi==yh)[0])]
        pl = pi[int(numpy.where(yi==yl)[0])]
        itera = 0
        
        while (yl>THRESHOLD):
            try:
                p_not_h = [x for x in pi if (x!=ph).all()]
                y_not_h = [x for x in yi if (x!=yh)]
                pbar = self.centroid(p_not_h)
                pstar = self.reflection(pbar,ph,ALPHA)        
                ystar = f(pstar)
            
                if ystar < yl:
                    
                    pstarstar = self.expansion(pbar, pstar, GAMMA)
                    ystarstar = f(pstarstar)
                    
                    if ystarstar < yl:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstarstar
                        yi[yi==yh] = ystarstar
                    else:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                        yi[yi==yh] = ystar
                else:
                    y_not_h.append(ystar)
                    if ystar == numpy.amax(numpy.array(y_not_h)):
                        if ystar > yh:
                            pass 
                        else:
                            pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                            yi[yi==yh] = ystar
                        
                        yh = numpy.amax(yi)
                        yl = numpy.amin(yi)
                        ph = pi[int(numpy.where(yi==yh)[0])]
                        pl = pi[int(numpy.where(yi==yl)[0])]
                        
                        pbar = self.centroid(p_not_h)
                        pstarstar = self.contraction(pbar, ph ,BETA)
                        
                        ystarstar = f(pstarstar)
                        
                        if ystarstar > yh:            
                            pi = [(x+pl)/2 for x in pi]
                            yi = numpy.array([f(p) for p in pi])

                        else:
                            pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstarstar
                            yi[yi==yh] = ystarstar
                    
                    else:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                        yi[yi==yh] = ystar

                yh = numpy.amax(yi)
                yl = numpy.amin(yi)
                ph = pi[int(numpy.where(yi==yh)[0])]
                pl = pi[int(numpy.where(yi==yl)[0])]
                itera += 1 

            except:
                if itera <=1:
                    raise ValueError("INCORRECT INITIAL VALUES GIVEN TO SIMPLEX")
                else:
                    print("SIMPLEX COMPLETE")
                break
        
        return (pl, itera)
