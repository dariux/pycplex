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
from pycplex import CPX, cplexcodes as C

Inf = C.CPX_INFBOUND
LP = C.CPXPROB_LP
QP = C.CPXPROB_QP
MILP = C.CPXPROB_MILP
MIQP = C.CPXPROB_MIQP

OPTIMAL = [C.CPX_STAT_OPTIMAL, C.CPXMIP_OPTIMAL, C.CPXMIP_OPTIMAL_TOL]
UNBOUNDED = [C.CPX_STAT_UNBOUNDED, C.CPXMIP_UNBOUNDED]
INFEASIBLE = [C.CPX_STAT_INFEASIBLE, C.CPXMIP_INFEASIBLE, C.CPXMIP_INForUNBD]


class MPProb(object):
    """Mathematical Programming problem
    By default, objsen = max, and obj is all zeros
    """
    def __init__(self, numrows, numcols):
        """
        Objective: all zeros, maximize
        Variable bounds: [-Inf, Inf]
        Probtype: LP
        """

        assert numcols > 0
        assert numrows >= 0
        self.numcols = numcols
        self.numrows = numrows
        self.probtype = LP # CPXPROB_LP, CPXPROB_QP or CPXPROB_MILP
        self.objsen = C.CPX_MAX # maximization problem

        self.obj = N.zeros((numcols,)) # numvars
        self.ctype = N.empty((numcols,), '|S1')
        self.ctype[:] = 'C' # {C,B,I} C = continuous
        self.lb = -Inf * N.ones((numcols,))
        self.ub = Inf * N.ones((numcols,))

        self.A = N.zeros((numrows,numcols))
        self.Q = None # for quadratic objective
        self.rhs = N.zeros((numrows,))
        # 'L' (<=), 'E' (=), 'G' (>=), 'R' (range)
        self.sense = N.empty((numrows,),'|S1')
        self.sense[:] = 'E' # set by default
        
        self.matbeg = self.matcnt = self.matind = None # integers
        self.matval = None # doubles
        self.qmatbeg = self.qmatcnt = self.qmatind = None # integers
        self.qmatval = None # doubles
        self.rngval = None # or N.zeros((numrows,))
        self._sparseA_in_sync = False
        self._sparseQ_in_sync = False

    @staticmethod
    def cplexsparse(A):
        """Convert matrix A to CPLEX sparse representation
        Thanks to Stephen Hartke for the code"""
        A = N.asarray(A, dtype=float)
        numrows, numcols = A.shape
        # n is the number of non-zero entries in A
        n = len(A.nonzero()[0])

        matval = N.empty((n,), dtype=float)
        matind = N.empty((n,), dtype=N.int32)
        matbeg = N.empty((numcols,), dtype=N.int32)
        matcnt = N.empty((numcols,), dtype=N.int32)
       
        i = 0
        for col in xrange(0, numcols):
            #if col % 100 == 0: 
            #    print col, " ", 
            matbeg[col] = i
            cur_row_count = 0
            for row in xrange(0, numrows):
                if A[row][col] != 0:
                    matval[i] = A[row][col]
                    matind[i] = row
                    i += 1
                    cur_row_count += 1
            matcnt[col] = cur_row_count
        assert(i == n)

        return {'matval':matval, 'matind':matind,
                'matbeg':matbeg, 'matcnt':matcnt}


    def setA(self,A):
        """Set constraints matrix A"""
        self.A = N.asarray(A, dtype=float)
        self.numrows, numcols = self.A.shape
        assert(self.numcols == numcols)
        self._sparseA_in_sync = False

    def setQ(self,Q):
        """Set quadratic objective matrix Q"""
        self.Q = N.asarray(Q, dtype=float)


    def setRHS(self,rhs):
        """Set rhs vector rhs"""
        self.rhs = N.asarray(rhs, dtype=float)


    def setSense(self,sense):
        """Set sense vector sense"""
        self.sense = N.asarray(sense, '|S1')


    def makeSparseA(self):
        """Convert constraint matrix A to CPLEX sparse representation"""
        if not self._sparseA_in_sync:
            s = MPProb.cplexsparse(self.A)
            self.matval = s['matval']
            self.matind = s['matind']
            self.matbeg = s['matbeg']
            self.matcnt = s['matcnt']
            self._sparseA_in_sync = True


    def makeSparseQ(self):
        """Convert matrix Q to CPLEX sparse representation"""
        if self.Q is not None and not self._sparseQ_in_sync:
            s = MPProb.cplexsparse(self.Q)
            self.qmatval = s['matval']
            self.qmatind = s['matind']
            self.qmatbeg = s['matbeg']
            self.qmatcnt = s['matcnt']
            self._sparseQ_in_sync = True

    def addConstraint(self, c):
        # c = {'indices', 'coeffs', 'sense', 'rhs'}
        assert(len(c['indices']) == len(c['coeffs']))
        Arow = N.zeros((1,self.numcols))
        for i,a in zip(c['indices'],c['coeffs']):
            Arow[0,i] = a
        self.addConstraintRows((Arow, [c['rhs']], [c['sense']]))
        
    def addComparisonConstraint(self, c):
        # c = {'index1', 'sense', 'index2'}
        if c['index1'] != c['index2']:
            Arow = N.zeros((1,self.numcols))
            Arow[0, c['index1']] = 1
            Arow[0, c['index2']] = -1
            self.addConstraintRows((Arow, [0], [c['sense']]))

    def addBoundConstraint(self, c):
        # c = {'index', 'sense', 'val'}
        Arow = N.zeros((1,self.numcols))
        Arow[0, c['index']] = 1
        self.addConstraintRows((Arow, [c['val']], [c['sense']]))

    def addConstraintRows(self, r):
        # Arows = r[0], rhs = r[1], sense = r[2]
        self.A = N.vstack([self.A, r[0]])
        self._sparseA_in_sync = False
        self.rhs = N.concatenate([self.rhs, r[1]])
        self.sense = N.concatenate([self.sense, r[2]])
        self.numrows += len(r[1])
        
    def removeLastConstraint(self):
        self.A = self.A[:-1:]
        self._sparseA_in_sync = False
        self.rhs = self.rhs[:-1]
        self.sense = self.sense[:-1]
        self.numrows -= 1


    def prepare(self):
        self.makeSparseA()
        self.makeSparseQ()
        
        assert(self.numcols > 0)
        assert(self.numrows >= 0)
        assert(self.objsen in (C.CPX_MIN, C.CPX_MAX))

        assert(self.obj != None and len(self.obj) == self.numcols)
        assert(self.lb != None and len(self.lb) == self.numcols)
        assert(self.ub != None and len(self.ub) == self.numcols)

        assert(self.rhs != None and len(self.rhs) == self.numrows)
        assert(self.sense != None and len(self.sense) == self.numrows)
        assert(self.rngval is None or len(self.rngval) == self.numrows)
        
        self.obj = N.asarray(self.obj, dtype=float)
        self.lb = N.asarray(self.lb, dtype=float)
        self.ub = N.asarray(self.ub, dtype=float)
        
        self.rhs = N.asarray(self.rhs, dtype=float)
        self.sense = N.asarray(self.sense, '|S1')
        if self.rngval is not None: self.rngval = N.asarray(self.rngval, dtype=float)

        if self.probtype in [MILP,MIQP]:
            assert(self.ctype is not None and len(self.ctype) == self.numcols)
            self.ctype = N.asarray(self.ctype, '|S1')
            # types: 'C' (continuous), 'B' (binary), 'I' (integer)

