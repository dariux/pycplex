import os
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

# Linux: change the following lines to reflect correct paths on your system 
cplexinc = '/p/learning/pkgs/ilog/cplex91/include'
cplexlib = '/p/learning/pkgs/ilog/cplex91/lib/x86-64_RHEL3.0_3.2/static_pic'
cplexlibrary = 'cplex'

# For Windows systems, uncomment the following:
#cplexinc = r'C:\ILog\Cplex91\include'
#cplexlib = r'C:\ILog\Cplex91\lib\msvc7\stat_sta'
#cplexlibrary = 'cplex91'

numpyinc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include/numpy')
CPXmodule = Extension('pycplex.CPX',
                      define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '6')],
                      include_dirs = [cplexinc, numpyinc],
                      library_dirs = [cplexlib],
                      libraries = [cplexlibrary],
                      sources = ['pycplex/CPX.c'])

setup(name = 'pycplex',
      version = '1.6',
      description = 'Python interface to CPLEX Callable Library',
      author = 'Darius Braziunas',
      author_email = 'darius@cs.toronto.edu',
      url = 'http://www.cs.toronto.edu/~darius/software/pycplex',
      long_description = '''
      This software provides a Python interface to the ILOG CPLEX Callable Library.
      It implements a subset of the most commonly used functions.''',
      ext_modules = [CPXmodule],
      packages = ['pycplex'],
      package_data = {'pycplex': ['*.txt']},
      )
       
