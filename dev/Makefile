# Change the following paths according to your installation

# Python
PYTHONINC = $(LOCALDIR)/include/python2.5
NUMPYINC = $(LOCALDIR)/lib/python2.5/site-packages/numpy/core/include/numpy 

# CPLEX
CPLEXINC = /pkgs/ilog/cplex91/include
CPLEXLIB = /pkgs/ilog/cplex91/lib/x86_RHEL3.0_3.2/static_pic

CFLAGS = -I$(PYTHONINC) -I$(NUMPYINC) -I$(CPLEXINC) -fPIC
LDFLAGS = -shared -L$(CPLEXLIB)

default: CPX.so

clean:
	rm -f *.o *.so

CPX.so: CPX.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o CPX.so CPX.o -lcplex
