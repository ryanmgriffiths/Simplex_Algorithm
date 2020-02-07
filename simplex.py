import numpy

class Simplex:

    def __init__(self):
        return None

    def f(self,input):
        '''Rosenbrock function for testing, input is array like (x,y)'''
        return (1-input[0])**2 + 100*(input[1]-input[0]**2)**2

    def centroid(self,array):
        centroid = numpy.zeros(2)
        idx = 0

        for x in array:
            centroid = centroid+x
            idx = idx + 1 
        
        centroid = centroid/idx
        return centroid

    def reflection(self,pbar, ph, alpha):
        return (1+alpha)*pbar - alpha*ph

    def expansion(self, pbar, pstar, gamma):
        return gamma*pstar + (1-gamma)*pbar

    def contraction(self, pbar, phigh, beta):
        return beta*phigh + (1-beta)*pbar

    def simplex(self, ALPHA, GAMMA, BETA, p0, THRESHOLD):
        # Assert inputs
        assert(ALPHA > 0)
        assert(GAMMA > 1)
        assert(BETA > 0 and BETA < 1)
        assert(len(p0) == len(p0[0])+1)
        
        pi = p0
        yi = [f(p) for p in pi]
        ph = pi[int(numpy.where(yi==numpy.amax(yi))[0])]
        yh = f(ph)
        pl = pi[int(numpy.where(yi==numpy.amin(yi))[0])]
        yl = f(pl)
        itera = 0
        
        while (yl>THRESHOLD):
            
            try:
                p_not_h = [x for x in pi if (x!=ph).all()]
                y_not_h = [f(x) for x in p_not_h]
                pbar = centroid(p_not_h)
                pstar = reflection(pbar,ph,ALPHA)        
                ystar = f(pstar)
            
                if ystar < yl:
                    
                    pstarstar = expansion(pbar, pstar, GAMMA)
                    ystarstar = f(pstarstar)
                    
                    if ystarstar < yl:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstarstar
                        yi = [f(p) for p in pi]
                    else:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                        yi = [f(p) for p in pi]
                else:
                    y_not_h.append(ystar)
                    if ystar == numpy.amax(numpy.array(y_not_h)):
                        if ystar > yh:
                            pass 
                        else:
                            pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                            yi = [f(p) for p in pi]
                        
                        ph = pi[int(numpy.where(yi==numpy.amax(yi))[0])]
                        pl = pi[int(numpy.where(yi==numpy.amin(yi))[0])]
                        yh = f(ph)
                        yl = f(pl)
                        pbar = centroid(p_not_h)
                        pstarstar = contraction(pbar, ph ,BETA)
                        
                        ystarstar = f(pstarstar)
                        
                        if ystarstar > yh:            
                            pi = [(x+pl)/2 for x in pi]
                            yi = [f(p) for p in pi]

                        else:
                            pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstarstar
                            yi = [f(p) for p in pi]
                    
                    else:
                        pi[[numpy.array_equal(ph,x) for x in pi].index(True)] = pstar
                        yi = [f(p) for p in pi]

                ph = pi[int(numpy.where(yi==numpy.amax(yi))[0])]
                yh = f(ph)
                pl = pi[int(numpy.where(yi==numpy.amin(yi))[0])]
                yl = f(pl)
                itera += 1 
            
            except:
                if itera <=1:
                    raise ValueError("INCORRECT VALUES FOR SIMPLEX PARAMS")
                else:
                    print("SIMPLEX COMPLETE")
                break
        return (pl, itera)

if __name__ == "__main__:
    ## MAIN PROGRAM
    s = Simplex()
    v = s.simplex( 1.5,
        1.3,
        0.8,
        [numpy.array((-1.2,1)) , numpy.array((-3,5)), numpy.array((1.57,2.6))], 
        0.00005)
    print("::::::::::::::FINAL RESULTS::::::::::::::")
    print("P =",v[0],", computed in %i iterations" % (v[1]))