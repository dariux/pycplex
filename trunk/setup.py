### You have to specify correct CPLEX include and library locations
### and the name of the cplex library (usually cplex or cplex91)
import sys

# Linux: change the following lines to reflect correct paths on your system 
cplexinc = '/p/learning/pkgs/ilog/cplex91/include'
cplexlib = '/p/learning/pkgs/ilog/cplex91/lib/x86-64_RHEL3.0_3.2/static_pic'
cplexlibrary = 'cplex'

# For Windows systems, make sure the following are correct:
if sys.platform == 'win32':
    cplexinc = r'C:\ILog\Cplex91\include'
    cplexlib = r'C:\ILog\Cplex91\lib\msvc7\stat_sta'
    cplexlibrary = 'cplex91'

import os
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

# patch distutils if it can't cope with the "classifiers" or
# "download_url" keywords
if sys.version < '2.2.3':
    from distutils.dist import DistributionMetadata
    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None

numpyinc = os.path.join(get_python_lib(plat_specific=1), 'numpy/core/include/numpy')
CPXmodule = Extension('pycplex.CPX',
                      include_dirs = [cplexinc, numpyinc],
                      library_dirs = [cplexlib],
                      libraries = [cplexlibrary],
                      sources = ['pycplex/CPX.c'])

setup(name = 'pycplex',
      version = '2.0a1',
      description = 'Python interface to CPLEX Callable Library',
      author = 'Darius Braziunas',
      author_email = 'darius@cs.toronto.edu',
      url = 'http://www.cs.toronto.edu/~darius/software/pycplex',
      long_description = '''
      This software provides a Python interface to the ILOG CPLEX Callable Library.
      It implements a subset of the most commonly used functions.''',
      license = 'The MIT Licence (http://www.opensource.org/licenses/mit-license.php)',
      platforms = ['any'],
      ext_modules = [CPXmodule],
      packages = ['pycplex'],
      package_data = {'pycplex': ['*.txt', 'CPX.c']},
      )
       
