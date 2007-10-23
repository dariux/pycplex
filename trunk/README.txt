
Version 1.2, October 2007

This software provides a Python interface to the CPLEX C Callable
Library. It implements a subset of the most commonly used
functions. If you need some other function, it should be very easy it
add it yourself (or drop me a line, and I will add it to the next
release).
The software is released under MIT Licence.

Darius Braziunas <darius@cs.toronto.edu>

For the latest version, notes, and acknowledgments, please check
http://www.cs.toronto.edu/~darius/software/pycplex


Installation
------------

You have to have numpy installed. See http://scipy.org

First, you need to compile the CPX.so module. Change the required
paths to CPLEX, Python and numpy in the Makefile, and type "make".

After, you can test the software by typing "python pycplex.py". I get
the following output:
=========
Created problem test1
Default variable names x1, x2 ... being created.
Default row names c1, c2 ... being created.
Tried aggregator 1 time.
LP Presolve eliminated 1 rows and 1 columns.
Aggregator did 2 substitutions.
All rows and columns eliminated.
Presolve time =   -0.00 sec.
[ 4. -1.  6.] 54.0

Created problem test2
Default variable names x1, x2 ... being created.
Default row names c1, c2 ... being created.
Tried aggregator 1 time.
MIP Presolve eliminated 1 rows and 1 columns.
Aggregator did 2 substitutions.
All rows and columns eliminated.
Presolve time =    0.01 sec.
[ 3.  0.  7.] 66.0
=========


Usage
-----

The functions test1() and test2() in pycplex.py provide examples on
how to call CPLEX. The class MPProb in pycplex.py is not strictly
necessary, but provides some convenient methods for preparing the
problem before calling CPLEX. If you have your input matrices and
vectors set up properly, you only need the wrappers in CPX.so.


List of CPLEX Callable library functions
----------------------------------------
Some documentation is obviously lacking. The functions are generally
equivalent to their CPLEX (CPX...) counterparts. 


chgcoef(...)

chgcoeflist(...)
    CPXchgcoeflist (env, lp, numcoefs, rowlist, collist, vallist) 
    changes a list of matrix coefficients 
    The list is prepared as a set of triples (i, j, value).

chgobj(...)
    CPX.chgobj(env,lp,cnt,indices,values) 
    cnt = len(indices), indices is int vector, values is double vector

    Changes the linear objective coefficients.

closeCPLEX(...)

copyctype(...)

copylp(...)

copyquad(...)

copymipstart(...)

createprob(...)

delrows(...)

addrows(...)

freeprob(...)

getmipobjval(...)

getmipx(...)

getobjval(...)

getstat(...)
    lpstat = CPX.getstat(env,lp) 
    Returns integer solution status of the problem after an LP, QP,
QCP, or MIP optimization.

getx(...)

lpopt(...)

qpopt(...)

mipopt(...)

newcols(...)
    CPXnewcols(env, lp, ccnt, obj, lb, ub, ctype) 
    ccnt is column count, obj, lb, ub, ctype are numpy vectors 
    In this function, we don't accept the column names as the last
argument.

newrows(...)
    CPXnewrows (env, lp, rcnt, rhs, sense, [rngval]) 
    Adds empty rows (without coefficients) 
    In this function, we don't accept the row names as the last
argument.

openCPLEX(...)
    env = CPX.openCPLEX()
    Initializes a CPLEX environment.

setintparam(...)

writeprob(...)
    CPX.writeprob(env,lp,filename)
    Format is deduced from file extension
    Formats: SAV, MPS, LP, REW, RMP, RLP.
 
