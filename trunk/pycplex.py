## Copyright (c) 2006-2008 Darius Braziunas

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
MIQP = CPXPROB_MIQP

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

        self.qmatbeg = self.qmatcnt = self.qmatind = None # integers
        self.qmatval = None # doubles

        
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

    # for quadratic matrix Q
    def qmakesparse(self,Q):
        """Convert matrix Q to CPLEX sparse representation"""
        self.Q = N.asarray(Q, dtype=float)
        self.qmatval = self.Q.transpose().copy()
        self.qmatval.shape = (self.numcols*self.numcols,)

        # CPLEX sparse index matrices
        self.qmatbeg = N.empty((self.numcols,), dtype=N.int32)
        self.qmatcnt = N.empty((self.numcols,), dtype=N.int32)
        self.qmatind = N.empty((self.numcols*self.numcols,), dtype=N.int32)

        for i in xrange(0,self.numcols):
            self.qmatbeg[i] = i * self.numcols
            for j in xrange(0, self.numcols):
                self.qmatind[i*self.numcols + j] = j
                self.qmatcnt[i] = self.numcols


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

        if self.probtype in [MILP,MIQP]:
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

    # change bounds of x1 and x3 (x1 <= 3, x3=5)
    #indices = N.array([0,2],dtype=N.int32)
    #lu = N.array(['L','B'])
    #bd = N.array([-3.0,5.0])
    #CPX.chgbds(env,lp,2,indices,lu,bd)

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
    
    CPX.freeprob(env,lp) # free the problem
    CPX.closeCPLEX(env) # close CPLEX environment

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
    
    CPX.freeprob(env,lp) # free the problem
    CPX.closeCPLEX(env) # close CPLEX environment


# thanks to Christina Fabritius Eskebaek for adding this test case
def test3():
    """Call CPLEX to solve MIQP

    Maximize
    obj: x + 2y + 3z - 0.5(33*x*x + 22*y*y + 11*z*z - 12*x*y - 23*y*z)
    Subject To
    c1: -x + y + z <= 20
    c2: x - 3y + z <= 30
    Bounds
    0 <= x <= 40
    0 <= y
    0 <= z
    Types
    """

    # define MPProb object
    p = MPProb(2,3) # numrows, numcols
    p.objsen = CPX_MAX
    p.obj = [1,2,3]
    A = [[-1,1,1], # constraint matrix
         [1,-3,1]]
    p.makesparse(A) # sets matval, matcnt, matind, matbeg
    p.rhs = [20,30]
    p.sense = ['L','L']
    p.lb[0], p.ub[0] = 0, 40
    p.lb[1] = 0
    p.lb[2], p.ub[2] = 0, 10

    p.probtype = QP

    # make the problem a quadratic problem
    Q = [[ -33,   6,     0   ],
	 [   6, -33,    11.5 ],
	 [   0,  11.5, -33   ]]
    p.qmakesparse(Q)

    # verify all input parameters are set up
    p.verify()

    # open CPLEX environment, set some parameters
    env = CPX.openCPLEX()
    CPX.setintparam(env, CPX_PARAM_SCRIND, CPX_ON) # print info to screen
    CPX.setintparam(env, CPX_PARAM_DATACHECK, CPX_ON)

    # create CPLEX problem, add objective and constraints to it
    lp = CPX.createprob(env, 'test3')
    
    CPX.copylp(env, lp, p.numcols, p.numrows, p.objsen, 
               p.obj, p.rhs, p.sense,
               p.matbeg, p.matcnt, p.matind, p.matval, 
               p.lb, p.ub) # rngval is optional
    CPX.copyquad(env, lp, p.qmatbeg, p.qmatcnt, p.qmatind, p.qmatval)

    # to make the problem mixed integer quadratic problem, do this:
    # p.probtype = MIQP
    # p.ctype = ['I', 'I', 'I']
    # CPX.copyctype(env, lp, p.ctype)

    status = CPX.writeprob(env, lp, "test3.lp") # write problem to file
                
    # solve the problem
    CPX.qpopt(env,lp) # can call lpopt, to


    # get the solution
    lpstat = CPX.getstat(env,lp)
    if lpstat in OPTIMAL:
        x, objval = CPX.getmipx(env,lp), CPX.getobjval(env,lp)
        print x, objval
    elif lpstat in UNBOUNDED: print "Solution unbounded"
    elif lpstat in INFEASIBLE: print "Solution infeasible"
    else: print "Solution code: ", lpstat
    print

    CPX.freeprob(env,lp) # free the problem
    CPX.closeCPLEX(env) # close CPLEX environment

if __name__ == "__main__":

    # run 3 tests
    test1()
    test2()
    test3()
