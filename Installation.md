Starting with Version 1.6 (April 2008), pycplex uses distutils. Therefore, to install pycplex:
  * Download and unzip pycplex package
  * Go to pycplex-2.x/ directory
  * Edit setup.py file and follow the instructions to specify correct CPLEX include and library locations as well as the name of the CPLEX library. You only need to change the following three lines according to your system setup:
```
cplexinc = '/pkgs/ilog/cplex111/include'
cplexlib = '/pkgs/ilog/cplex111/lib/x86-64_debian4.0_4.1/static_pic'
cplexlibrary = 'cplex'
```
  * Run
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

---

See http://www.cs.toronto.edu/~darius/software/pycplex/ for more information (some of it might be a bit outdated)