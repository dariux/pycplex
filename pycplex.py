## Copyright (c) 2006-2007 Darius Braziunas

## Permission is hereby granted, free of charge, to any person obtaining 
## a copy of this software and associated documentation files (the "Software"), 
## to deal in the Software without restriction, including without limitation the 
## rights to use, copy, modify, merge, publish, distribute, sublicense, 
## and/or sell copies of the Software, and to permit persons to whom 
## the Software is furnished to do so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in all 
## copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
## OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
## OTHER DEALINGS IN THE SOFTWARE.



import numpy as N
import CPX
from cplex_const import *

Inf = CPX_INFBOUND
LP = CPXPROB_LP
QP = CPXPROB_QP
MILP = CPXPROB_MILP

OPTIMAL = [CPX_STAT_OPTIMAL, CPXMIP_OPTIMAL, CPXMIP_OPTIMAL_TOL]
UNBOUNDED = [CPX_STAT_UNBOUNDED, CPXMIP_UNBOUNDED]
INFEASIBLE = [CPX_STAT_INFEASIBLE, CPXMIP_INFEASIBLE, CPXMIP_INForUNBD]



class MPProb:
    """Mathematical Programming problem"""
    def __init__(self, numrows, numcols):

        assert numcols > 0
        assert numrows > 0
        self.numcols = numcols
        self.numrows = numrows
        self.probtype = LP # CPXPROB_LP, CPXPROB_QP or CPXPROB_MILP
        self.objsen = CPX_MAX # maximization problem
        
        self.obj = N.ones((numcols,)) # numvars
        self.rhs = N.zeros((numrows,))
        # 'L' (<=), 'E' (=), 'G' (>=), 'R' (range)
        self.sense = N.empty((numrows,),'|S1')
        self.sense[:] = 'E' # set by default
        
        self.matbeg = self.matcnt = self.matind = None # integers
        self.matval = None # doubles
        
        self.lb = -Inf * N.ones((numcols,))
        self.ub = Inf * N.ones((numcols,))
        self.rngval = None # or N.zeros((numrows,))

        
    def makesparse(self,A):
        """Convert constraint matrix A to CPLEX sparse representation"""
        self.A = N.asarray(A, dtype=float)
        self.matval = self.A.transpose().copy()
        self.matval.shape = (self.numrows*self.numcols,)

        # CPLEX sparse index matrices
        self.matbeg = N.empty((self.numcols,), dtype=N.int32)
        self.matcnt = N.empty((self.numcols,), dtype=N.int32)
        self.matind = N.empty((self.numrows*self.numcols,), dtype=N.int32)

        for i in xrange(0,self.numcols):
            self.matbeg[i] = i * self.numrows
            for j in xrange(0, self.numrows):
                self.matind[i*self.numrows + j] = j
                self.matcnt[i] = self.numrows

    def verify(self):
        assert(self.numcols > 0)
        assert(self.numrows > 0)
        assert(self.objsen in (CPX_MIN, CPX_MAX))

        assert(self.obj != None and len(self.obj) == self.numcols)

        assert(self.rhs != None and len(self.rhs) == self.numrows)
        assert(self.sense != None and len(self.sense) == self.numrows)
        assert(self.lb != None and len(self.lb) == self.numcols)
        assert(self.ub != None and len(self.ub) == self.numcols)
        assert(self.rngval is None or len(self.rngval) == self.numrows)
        
        self.obj = N.asarray(self.obj, dtype=float)
        
        self.rhs = N.asarray(self.rhs, dtype=float)
        self.sense = N.asarray(self.sense, '|S1')
        self.lb = N.asarray(self.lb, dtype=float)
        self.ub = N.asarray(self.ub, dtype=float)
        if self.rngval is not None: self.rngval = N.asarray(self.rngval, dtype=float)

        if self.probtype == MILP:
            assert(self.ctype is not None and len(self.ctype) == self.numcols)
            self.ctype = N.asarray(self.ctype, '|S1')
            # types: 'C' (continuous), 'B' (binary), 'I' (integer)
    
def test1():
    """Call CPLEX to solve LP

    Minimize
    obj: x + 4 y + 9 z
    Subject To
    c1: y + x <= 5
    c2: z + x >= 10
    c3: z - y = 7
    Bounds
    x <= 4
    -1 <= y <= 1
    """

    # define MPProb object
    p = MPProb(3,3) # numcols, numrows
    p.objsen = CPX_MIN
    p.obj = [1,4,9]
    A = [[1,1,0], # constraint matrix
         [1,0,1],
         [0,-1,1]]
    p.makesparse(A) # sets matval, matcnt, matind, matbeg
    p.rhs = [5,10,7]
    p.sense = ['L','G','E']
    p.ub[0] = 4
    p.lb[1], p.ub[1] = -1,1

    # verify all input parameters are set up
    p.verify()
    
    # open CPLEX environment, set some parameters
    env = CPX.openCPLEX()
    CPX.setintparam(env, CPX_PARAM_SCRIND, CPX_ON) # print info to screen
    CPX.setintparam(env, CPX_PARAM_DATACHECK, CPX_ON)

    # create CPLEX problem, add objective and constraints to it
    lp = CPX.createprob(env, 'test1')
    CPX.copylp(env, lp, p.numcols, p.numrows, p.objsen, 
               p.obj, p.rhs, p.sense,
               p.matbeg, p.matcnt, p.matind, p.matval, 
               p.lb, p.ub) #rngval optional
    status = CPX.writeprob(env, lp, "test1.lp") # write problem to file                

    # solve the problem
    CPX.lpopt(env,lp)

    # get the solution
    lpstat = CPX.getstat(env,lp)
    if lpstat in OPTIMAL:
        x, objval = CPX.getx(env,lp), CPX.getobjval(env,lp)
        print x, objval
    elif lpstat in UNBOUNDED: print "Solution unbounded"
    elif lpstat in INFEASIBLE: print "Solution infeasible"
    else: print "Solution code: ", lpstat
    print
    
    # free the problem
    CPX.freeprob(env,lp)
    # close CPLEX environment
    CPX.closeCPLEX(env)

def test2():
    """Call CPLEX to solve MILP

    Minimize
    obj: x + 4 y + 9 z
    Subject To
    c1: y + x <= 5
    c2: z + x >= 10
    c3: z - y = 7
    Bounds
    x <= 4
    -1 <= y <= 1
    Types
    x is continuous, y and z boolean
    """

    # define MPProb object
    p = MPProb(3,3) # numcols, numrows
    p.objsen = CPX_MIN
    p.obj = [1,4,9]
    A = [[1,1,0], # constraint matrix
         [1,0,1],
         [0,-1,1]]
    p.makesparse(A) # sets matval, matcnt, matind, matbeg
    p.rhs = [5,10,7]
    p.sense = ['L','G','E']
    p.ub[0] = 4
    p.lb[1], p.ub[1] = -1,1

    # make the problem mixed integer problem
    p.probtype = MILP
    p.ctype = ['I','B','I']

    # verify all input parameters are set up
    p.verify()
    
    # open CPLEX environment, set some parameters
    env = CPX.openCPLEX()
    CPX.setintparam(env, CPX_PARAM_SCRIND, CPX_ON) # print info to screen
    CPX.setintparam(env, CPX_PARAM_DATACHECK, CPX_ON)

    # create CPLEX problem, add objective and constraints to it
    lp = CPX.createprob(env, 'test2')
    
    CPX.copylp(env, lp, p.numcols, p.numrows, p.objsen, 
               p.obj, p.rhs, p.sense,
               p.matbeg, p.matcnt, p.matind, p.matval, 
               p.lb, p.ub) # rngval is optional
    CPX.copyctype(env, lp, p.ctype)
    status = CPX.writeprob(env, lp, "test2.lp") # write problem to file
                
    # solve the problem
    CPX.mipopt(env,lp) # can call lpopt, too

    # get the solution
    lpstat = CPX.getstat(env,lp)
    if lpstat in OPTIMAL:
        # can call getx, getobjval, too
        x, objval = CPX.getmipx(env,lp), CPX.getmipobjval(env,lp)
        print x, objval
    elif lpstat in UNBOUNDED: print "Solution unbounded"
    elif lpstat in INFEASIBLE: print "Solution infeasible"
    else: print "Solution code: ", lpstat
    print
    
    # free the problem
    CPX.freeprob(env,lp)
    # close CPLEX environment
    CPX.closeCPLEX(env)


if __name__ == "__main__":

    # run 2 tests
    test1()
    test2()
