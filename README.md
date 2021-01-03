This software provides a Python interface to the CPLEX C Callable
Library. It implements a subset of the most commonly used
functions. If you need some other function, it should be easy it
add it yourself (or drop me a line, and I will add it to the next
release).

The software is released under MIT Licence.

CPLEX® is a registered trademark of ILOG.

# Installation
Starting with Version 1.6 (April 2008), pycplex uses distutils. Therefore, to install pycplex:
  * Download and unzip pycplex package
  * Go to pycplex-2.x/ directory
  * Edit setup.py file and follow the instructions to specify correct CPLEX include and library locations as well as the name of the CPLEX library. You only need to change the following three lines according to your system setup:
```
cplexinc = '/pkgs/ilog/cplex111/include'
cplexlib = '/pkgs/ilog/cplex111/lib/x86-64_debian4.0_4.1/static_pic'
cplexlibrary = 'cplex'
```

Run
```
python setup.py install
```
  * If you get no errors, the package will be installed in <python path>/site-packages/pycplex location
  * Go to <python path>/site-packages/pycplex directory and test the installation:
```
python tests.py
```

Your output should be similar to:
```
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

            Created problem test3
            Default variable names x1, x2 ... being created.
            Default row names c1, c2 ... being created.
            Number of nonzeros in lower triangle of Q = 2
            Using Approximate Minimum Degree ordering
            ...
            [ 0.05133332  0.11566658  0.13121714] 0.338158944865
```

# Notes

You have to have numpy installed. See scipy.org

You need to compile the CPX.so module. Change the required paths to CPLEX, Python and numpy in the Makefile, and type "make".

You can test the software by typing ```python pycplex.py```. I get the following output:
```
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

            Created problem test3
            Default variable names x1, x2 ... being created.
            Default row names c1, c2 ... being created.
            Number of nonzeros in lower triangle of Q = 2
            Using Approximate Minimum Degree ordering
            ...
            [ 0.05133332  0.11566658  0.13121714] 0.338158944865
```

## Note 
If you're getting "CPLEX Error 1228: Count entry 0 indicates overlapping entries", the problem is with 64bit machines. numpy/numarray internally uses 64 bit for int arrays on x86-64, but CPLEX expects the C "int" type which is only 32 bits there. Version 1.1 contains a fix in the pycplex.py by allocating the integer arrays as explicit 32-bit versions.
```
        self.matbeg = N.empty((self.numcols,), dtype=N.int32)
	    self.matcnt = N.empty((self.numcols,), dtype=N.int32)
	    self.matind = N.empty((self.numrows*self.numcols,), dtype=N.int32)
```
Version 1.2 checks for the right item size (only integer) numpy arrays in CPX.c code.
The rule is: numpy integer arrays passed to CPLEX have to be of size sizeof(int) on your platform. This usually means that you have to explicitly initialize them as N.int32 on most 64-bit platforms (for future code compatibility, do the same on 32-bit platforms, too).
Thanks to Christoph Lampert for the fix.

## Note
Please make sure you have at least one constraint row, otherwise you will run into strange memory errors.


# Usage
The functions test1() and test2() in pycplex.py provide examples on how to call CPLEX. The class MPProb in pycplex.py (and the whole module pycplex.py) is not strictly necessary, but provides some convenient methods for preparing the problem before calling CPLEX. If you have your input matrices and vectors set up properly, you only need the wrappers in CPX.so.
List of implemented (wrapped) CPLEX Callable library functions
Some documentation is obviously lacking. The functions are generally equivalent to their CPLEX (CPX...) counterparts.

```
addsos(...)

copysos(...)

delsetsos(...)

chgbds(...)

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

copystart(...)

createprob(...)

delrows(...)

addrows()

freeprob(...)

getmipobjval(...)

getmipx(...)

getobjval(...)

getstat(...)
    lpstat = CPX.getstat(env,lp) 
    Returns integer solution status of the problem after an LP, QP,
QCP, or MIP optimization.

getx(...)

getslack(...)

getpi(...)

getdj(...)

getbase(...)

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

setdblparam(...)

writeprob(...)
    CPX.writeprob(env,lp,filename)
    Format is deduced from file extension
    Formats: SAV, MPS, LP, REW, RMP, RLP.
```