/***
Copyright (c) 2006-2007 Darius Braziunas (darius@cs.toronto.edu)
http://www.cs.toronto.edu/~darius/software/pycplex

Permission is hereby granted, free of charge, to any person obtaining 
a copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the 
rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom 
the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.
***/


#include <math.h>
#include "Python.h"
#include "arrayobject.h"
#include "ilcplex/cplex.h"


#define CPXSTATUS if (status) setCPLEXerrorstring(env, status);

static PyObject *openCPLEX(PyObject *self);
static PyObject *closeCPLEX(PyObject *self, PyObject *args);
static PyObject *createprob(PyObject *self, PyObject *args);
static PyObject *writeprob(PyObject *self, PyObject *args);
static PyObject *freeprob(PyObject *self, PyObject *args);
static PyObject *setintparam(PyObject *self, PyObject *args);
static PyObject *copylp(PyObject *self, PyObject *args);
static PyObject *copyquad(PyObject *self, PyObject *args); 
static PyObject *copyctype(PyObject *self, PyObject *args);
static PyObject *copymipstart(PyObject *self, PyObject *args);
static PyObject *chgobj(PyObject *self, PyObject *args);
static PyObject *newcols(PyObject *self, PyObject *args);
static PyObject *newrows(PyObject *self, PyObject *args);
static PyObject *delrows(PyObject *self, PyObject *args);
static PyObject *addrows(PyObject *self, PyObject *args); 
static PyObject *chgcoeflist(PyObject *self, PyObject *args);
static PyObject *chgcoef(PyObject *self, PyObject *args);
static PyObject *lpopt(PyObject *self, PyObject *args);
static PyObject *getx(PyObject *self, PyObject *args);
static PyObject *getobjval(PyObject *self, PyObject *args);
static PyObject *getstat(PyObject *self, PyObject *args);

static void checkintsize(PyArrayObject *intarray);
static void setCPLEXerrorstring(CPXENVptr env, int status);

/* Docs */
static char openCPLEX_doc[] = 
"env = CPX.openCPLEX()\n\
Initializes a CPLEX environment.";

static char writeprob_doc[] = 
"CPX.writeprob(env,lp,filename)\n\
Format is deduced from file extension\n\
Formats: SAV, MPS, LP, REW, RMP, RLP.";

static char newrows_doc[] = 
"CPXnewrows (env, lp, rcnt, rhs, sense, [rngval]) \n\
Adds empty rows (without coefficients) \n\
In this function, we don't accept the row names as the last argument.";

static char newcols_doc[] = 
"CPXnewcols(env, lp, ccnt, obj, lb, ub, ctype) \n\
ccnt is column count, obj, lb, ub, ctype are numpy vectors \n\
In this function, we don't accept the column names as the last argument.";

static char chgcoeflist_doc[] = 
"CPXchgcoeflist (env, lp, numcoefs, rowlist, collist, vallist) \n\
changes a list of matrix coefficients \n\
The list is prepared as a set of triples (i, j, value).";

static char chgobj_doc[] = 
"CPX.chgobj(env,lp,cnt,indices,values) \n\
cnt = len(indices), indices is int vector, values is double vector \n\
Changes the linear objective coefficients.";

static char getstat_doc[] = 
"lpstat = CPX.getstat(env,lp) \n\
Returns integer solution status of the problem after an LP, QP, QCP, or MIP optimization.";


static PyMethodDef methods[] = {
  {"openCPLEX", (PyCFunction)openCPLEX, METH_NOARGS, openCPLEX_doc},
  {"closeCPLEX", closeCPLEX, METH_VARARGS},
  {"createprob", createprob, METH_VARARGS},
  {"writeprob", writeprob, METH_VARARGS, writeprob_doc},
  {"freeprob", freeprob, METH_VARARGS},
  {"setintparam", setintparam, METH_VARARGS},
  {"copylp", copylp, METH_VARARGS},
  {"copyquad", copyquad, METH_VARARGS}, 
  {"copyctype", copyctype, METH_VARARGS},
  {"copymipstart", copymipstart, METH_VARARGS},
  {"chgobj", chgobj, METH_VARARGS, chgobj_doc},
  {"newcols", newcols, METH_VARARGS, newcols_doc},
  {"newrows", newrows, METH_VARARGS, newrows_doc},
  {"delrows", delrows, METH_VARARGS},
  {"addrows", addrows, METH_VARARGS}, 
  {"chgcoeflist", chgcoeflist, METH_VARARGS, chgcoeflist_doc},
  {"chgcoef", chgcoef, METH_VARARGS},
  {"lpopt", lpopt, METH_VARARGS},
  {"mipopt", lpopt, METH_VARARGS}, // same function as lpopt
  {"qpopt", lpopt, METH_VARARGS}, // same function as lpopt
  {"getx", getx, METH_VARARGS},
  {"getmipx", getx, METH_VARARGS}, // same as getx
  {"getobjval", getobjval, METH_VARARGS},
  {"getmipobjval", getobjval, METH_VARARGS}, // same as getobjval
  {"getstat", getstat, METH_VARARGS, getstat_doc},
  {NULL, NULL} /* sentinel */
};

void initCPX(void)
{
  (void) Py_InitModule("CPX", methods);
}

PyObject *openCPLEX(PyObject *self) {
  /* Open CPLEX environment */
  int status;

  CPXENVptr env;
  env = CPXopenCPLEX(&status);
  CPXSTATUS;
  
  // http://www.python.org/doc/2.3/api/cObjects.html
  return PyCObject_FromVoidPtr(env, NULL);
}

PyObject *closeCPLEX(PyObject *self, PyObject *args) {
  /* Close CPLEX environment 
     http://www.cs.toronto.edu/~darius/ref/cplex/
     refcallablelibrary/html/functions/CPXcloseCPLEX.html
  */
  int status;

  CPXENVptr env;
  PyObject *pyenv;
  if (!PyArg_ParseTuple(args, "O", &pyenv))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);

  status = CPXcloseCPLEX(&env); 
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *createprob(PyObject *self, PyObject *args) {
  /* CPXcreateprob */
  int status;
  char *probname;
  CPXLPptr lp;

  CPXENVptr env; 
  PyObject *pyenv;
  if (!PyArg_ParseTuple(args, "Os", &pyenv, &probname))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);

  lp = CPXcreateprob(env, &status, probname);
  CPXSTATUS;
  printf("Created problem %s\n", probname);

  return PyCObject_FromVoidPtr(lp, NULL);
}

PyObject *writeprob(PyObject *self, PyObject *args) {
  /* CPXwriteprob - format is deduced from file extension*/
  int status;
  char *filename_str;
  CPXLPptr lp;
  CPXENVptr env;
  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OOs", &pyenv, &pylp, &filename_str))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXwriteprob(env,lp,filename_str,NULL);
  CPXSTATUS;
  
  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *freeprob(PyObject *self, PyObject *args) {
  /* CPXfreeprob */
  int status;
  CPXLPptr lp;
  CPXENVptr env;
  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OO", &pyenv, &pylp))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  status = CPXfreeprob (env, &lp);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *setintparam(PyObject *self, PyObject *args) {
  /* CPXsetintparam */
  int status;
  int whichparam, newvalue;

  CPXENVptr env; 
  PyObject *pyenv;
  if (!PyArg_ParseTuple(args, "Oii", &pyenv, &whichparam, &newvalue))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);

  status = CPXsetintparam (env, whichparam, newvalue);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *copymipstart(PyObject *self, PyObject *args) {
  /* CPXcopymipstart - copy MIP starting values */

  int status;
  int cnt;
  CPXENVptr env;
  CPXLPptr lp;

  PyObject *pyenv, *pylp;
  PyArrayObject *indices, *value;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiO!O!",
			&pyenv, &pylp, &cnt, 
			&PyArray_Type, &indices,
			&PyArray_Type, &value))
    return NULL;

  checkintsize(indices);

  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXcopymipstart(env, lp, cnt, (int *)indices->data, (double *)value->data);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *copylp(PyObject *self, PyObject *args) {
  /* CPXcopylp - copy LP into CPLEX environment */
  /* rngval is optional */
  /* int CPXcopylp(CPXCENVptr env, CPXLPptr lp, int numcols, int numrows, 
     int objsen, const double * obj, const double * rhs, const char * sense, 
     const int * matbeg, const int * matcnt, const int * matind, 
     const double * matval, const double * lb, 
     const double * ub, const double * rngval) */
  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int numcols, numrows, objsen; // objsen is 1 (CPX_MIN) or -1 (CPX_MAX)
  /* Python input variables */
  PyObject *pyenv, *pylp;
  PyArrayObject *obj, *rhs, *sense;
  PyArrayObject *matbeg, *matcnt, *matind, *matval, *lb, *ub, *rngval=NULL;
  double *rngvaldata = NULL;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiiiO!O!O!O!O!O!O!O!O!|O!", 
			&pyenv, &pylp,
			&numcols, &numrows, &objsen, 
			&PyArray_Type, &obj, 
			&PyArray_Type, &rhs, 
			&PyArray_Type, &sense,
			&PyArray_Type, &matbeg, 
			&PyArray_Type, &matcnt, 
			&PyArray_Type, &matind, 
			&PyArray_Type, &matval, 
			&PyArray_Type, &lb, 
			&PyArray_Type, &ub, 
			&PyArray_Type, &rngval))
    return NULL;

  checkintsize(matbeg);
  checkintsize(matcnt);
  checkintsize(matind);

  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  /*
  numcols = obj->dimensions[0]; // # of variables 
  numrows = rhs->dimensions[0]; // # of rows
  */

  if (rngval != NULL) rngvaldata = (double *)rngval->data;
  
  status = CPXcopylp(env, lp, numcols, numrows, objsen, 
		     (double *)obj->data, 
		     (double *)rhs->data, 
		     (char *)sense->data, 
		     (int *)matbeg->data, 
		     (int *)matcnt->data, 
		     (int *)matind->data,
		     (double *)matval->data,
		     (double *)lb->data, 
		     (double *)ub->data, 
		     rngvaldata);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *copyquad(PyObject *self, PyObject *args) { 
  /* CPXcopyquad - copy qp data into CPLEX environment */ 
  /* int CPXcopyquad(CPXCENVptr env, CPXLPptr lp, const int * qmatbeg, const int * qmatcnt, const int * qmatind, const double * qmatval) */ 
  int status; 
 
  CPXENVptr env; 
  CPXLPptr lp; 
  int numcols, numrows, objsen; // objsen is 1 (CPX_MIN) or -1 (CPX_MAX) 
  /* Python input variables */ 
  PyObject *pyenv, *pylp; 
  PyArrayObject *qmatbeg, *qmatcnt, *qmatind, *qmatval; 
  import_array(); 
 
  if (!PyArg_ParseTuple(args, "OOO!O!O!O!",  
                        &pyenv, &pylp, 
                        &PyArray_Type, &qmatbeg,  
                        &PyArray_Type, &qmatcnt,  
                        &PyArray_Type, &qmatind,  
                        &PyArray_Type, &qmatval)) 
    return NULL; 
 
  checkintsize(qmatbeg);
  checkintsize(qmatcnt);
  checkintsize(qmatind);

  env = PyCObject_AsVoidPtr(pyenv); 
  lp = PyCObject_AsVoidPtr(pylp); 
 
  status = CPXcopyquad(env, lp,  
		       (int *)qmatbeg->data,  
		       (int *)qmatcnt->data,  
		       (int *)qmatind->data, 
		       (double *)qmatval->data); 
  CPXSTATUS; 
 
  Py_INCREF(Py_None); 
  return Py_None; 
} 

PyObject *copyctype(PyObject *self, PyObject *args) {
  /* CPXcopyctype - copy ctype and convert prob to MIP */

  int status;

  CPXENVptr env;
  CPXLPptr lp;

  PyObject *pyenv, *pylp;
  PyArrayObject *ctype;
  import_array();

  if (!PyArg_ParseTuple(args, "OOO!",
			&pyenv, &pylp, 
			&PyArray_Type, &ctype))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXcopyctype(env, lp, (char *)ctype->data);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *newcols(PyObject *self, PyObject *args) {
  /* CPXnewcols adds empty columns */
  /* CPXnewcols(env, lp, ccnt, obj, lb, ub, ctype)
     In this function, we don't accept the column names as the last argument
  */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int ccnt;  

  PyObject *pyenv, *pylp;
  PyArrayObject *obj, *lb, *ub, *ctype;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiO!O!O!O!",
			&pyenv, &pylp, &ccnt,
			&PyArray_Type, &obj,
			&PyArray_Type, &lb,
			&PyArray_Type, &ub,
			&PyArray_Type, &ctype))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXnewcols(env,lp,ccnt,
		      (double *)obj->data,
		      (double *)lb->data,
		      (double *)ub->data,
		      (char *)ctype->data, NULL);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *newrows(PyObject *self, PyObject *args) {
  /* CPXnewrows adds empty rows (without coefficients) */
  /* CPXnewrows (env, lp, rcnt, rhs, sense, [rngval])
     In this function, we don't accept the row names as the last argument
  */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int rcnt;  

  PyObject *pyenv, *pylp;
  PyArrayObject *rhs, *sense, *rngval=NULL;
  double *rngvaldata=NULL;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiO!O!|O!",
			&pyenv, &pylp, &rcnt,
			&PyArray_Type, &rhs,
			&PyArray_Type, &sense,
			&PyArray_Type, &rngval))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  if (rngval != NULL) rngvaldata = (double *)rngval->data;

  status = CPXnewrows(env, lp, rcnt, 
		      (double *)rhs->data, 
		      (char *)sense->data, 
		      rngvaldata, NULL);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *delrows(PyObject *self, PyObject *args) {
  /* CPXdelrows (env, lp, begin, end) */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int begin, end;  

  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OOii", &pyenv, &pylp, &begin, &end))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  status = CPXdelrows(env, lp, begin, end);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *chgcoef(PyObject *self, PyObject *args) {
  /* Change a single coefficient in the constraint matrix, 
     linear objective coefficients, right-hand side or 
     ranges of a CPLEX problem object
  */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int i,j;  
  double newvalue;

  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OOiid",
			&pyenv, &pylp, &i, &j, &newvalue))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXchgcoef(env, lp, i, j, newvalue);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *addrows(PyObject *self, PyObject *args) { 
  /* int CPXPUBLIC CPXaddrows(CPXCENVptr env, CPXLPptr lp, 
   * int ccnt, int rcnt, int nzcnt, const double * rhs,  
   * const char * sense, const int * rmatbeg, const int * rmatind, 
   * const double * rmatval, char ** colname, char ** rowname) 
   */ 
 
  int status; 
 
  CPXENVptr env; 
  CPXLPptr lp; 
  int rcnt;   
 
  PyObject *pyenv, *pylp; 
  int numCols, numRows, nonZeroEntries; 
  PyArrayObject *rhs, *sense; 
  PyArrayObject *rmatbeg, *rmatind, *rmatval; 
  import_array(); 
 
  if (!PyArg_ParseTuple(args, "OOiiiO!O!O!O!O!", 
                        &pyenv, &pylp, &numCols, &numRows, &nonZeroEntries, 
                        &PyArray_Type, &rhs, 
                        &PyArray_Type, &sense, 
                        &PyArray_Type, &rmatbeg, 
                        &PyArray_Type, &rmatind, 
                        &PyArray_Type, &rmatval)) 
    return NULL; 

  checkintsize(rmatbeg);
  checkintsize(rmatind);

  env = PyCObject_AsVoidPtr(pyenv); 
  lp = PyCObject_AsVoidPtr(pylp); 
   
  status = CPXaddrows(env, lp, numCols, numRows,  
		      nonZeroEntries, 
                      (double *)rhs->data,  
                      (char *)sense->data,  
                      (int *)rmatbeg->data,  
                      (int *)rmatind->data, 
                      (double *)rmatval->data, 
		      NULL,NULL); 
  CPXSTATUS; 
 
  return Py_BuildValue("i", status); 
} 

PyObject *chgcoeflist(PyObject *self, PyObject *args) {
  /* CPXchgcoeflist changes a list of matrix coefficients */
  /* The list is prepared as a set of triples (i, j, value) */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int numcoefs;  

  PyObject *pyenv, *pylp;
  PyArrayObject *rowlist, *collist, *vallist;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiO!O!O!",
			&pyenv, &pylp, &numcoefs,
			&PyArray_Type, &rowlist,
			&PyArray_Type, &collist,
			&PyArray_Type, &vallist))
    return NULL;

  checkintsize(rowlist);
  checkintsize(collist);

  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  status = CPXchgcoeflist (env, lp, numcoefs, 
			   (int *)rowlist->data, 
			   (int *)collist->data, 
			   (double *)vallist->data);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *chgobj(PyObject *self, PyObject *args) {
  /* CPXchgobj - change the linear objective coefficients */

  int status;

  CPXENVptr env;
  CPXLPptr lp;
  int cnt;

  PyObject *pyenv, *pylp;
  PyArrayObject *indices, *values;
  import_array();

  if (!PyArg_ParseTuple(args, "OOiO!O!",
			&pyenv, &pylp, &cnt,
			&PyArray_Type, &indices,
			&PyArray_Type, &values))
    return NULL;

  checkintsize(indices);
  //printf("itemsize %d, sizeofint %d, size(int) %d \n", PyArray_ITEMSIZE(indices),SIZEOF_INT, sizeof(int));

  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  status = CPXchgobj(env,lp,cnt,(int *)indices->data,(double *)values->data);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}


PyObject *lpopt(PyObject *self, PyObject *args) {
  /* CPXlpopt or CPXmipopt*/
  int status, probtype;

  CPXENVptr env;
  CPXLPptr lp;
  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OO", &pyenv, &pylp))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  probtype = CPXgetprobtype(env, lp);
  if (probtype == CPXPROB_MILP)
    status = CPXmipopt(env,lp);
  else if (probtype == CPXPROB_QP)
    status = CPXqpopt(env,lp);
  else
    status = CPXlpopt(env,lp);
  CPXSTATUS;

  Py_INCREF(Py_None);
  return Py_None;
}


PyObject *getx(PyObject *self, PyObject *args) {
  /* CPXgetx or CPXgetmipx (works for MIPs too)
     status = CPXgetx (env, lp, x, 0, CPXgetnumcols(env, lp)-1);
     In Python: x = getx(env,lp,begin=0,end=numcols-1)
     If no begin and end specified, return full x
   */
  int status, probtype;

  double *x;
  PyArrayObject *pyx;
  int dimensions[1];

  CPXENVptr env;
  CPXLPptr lp;
  int begin=-1, end=-1;
  PyObject *pyenv, *pylp;

  import_array();
  // |ii means ii is optional
  if (!PyArg_ParseTuple(args, "OO|ii", &pyenv, &pylp, &begin, &end))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);

  if (begin == -1) { // set to default values
    begin = 0;
    end = CPXgetnumcols(env, lp)-1;
  }
  x = malloc(sizeof(double)*(end-begin+1));
  probtype = CPXgetprobtype(env, lp);
  if (probtype == CPXPROB_MILP)
    status = CPXgetmipx (env, lp, x, begin, end);
  else
    status = CPXgetx (env, lp, x, begin, end);
  CPXSTATUS;

  dimensions[0] = end-begin+1;
  pyx = (PyArrayObject *)PyArray_FromDimsAndData(1,dimensions,
						 PyArray_DOUBLE,(char *)x);
  return PyArray_Return(pyx);
}

PyObject *getobjval(PyObject *self, PyObject *args) {
  /* CPXgetobjval or CPXgetmipobjval*/
  int status, probtype;

  CPXENVptr env;
  CPXLPptr lp;
  double objval;
  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OO", &pyenv, &pylp))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  probtype = CPXgetprobtype(env, lp);
  if (probtype == CPXPROB_MILP)
    status = CPXgetmipobjval(env,lp, &objval);
  else
    status = CPXgetobjval(env,lp, &objval);
  CPXSTATUS;

  return Py_BuildValue("d", objval);
}

PyObject *getstat(PyObject *self, PyObject *args) {
  /* CPXgetstat*/
  int lpstat;

  CPXENVptr env;
  CPXLPptr lp;

  PyObject *pyenv, *pylp;
  if (!PyArg_ParseTuple(args, "OO", &pyenv, &pylp))
    return NULL;
  env = PyCObject_AsVoidPtr(pyenv);
  lp = PyCObject_AsVoidPtr(pylp);
  
  lpstat = CPXgetstat(env,lp);
  
  return Py_BuildValue("i", lpstat);
}



void checkintsize(PyArrayObject *intarray) {
  char errorstring[1024];
  if (PyArray_ITEMSIZE(intarray) != sizeof(int)) {
    sprintf(errorstring, "The item size of numpy integer array is %d, \n \
    but CPLEX is expecting %d. Most likely, you can fix the problem by explicitly intializing \n \
    the numpy array as numpy.array(...,dtype=numpy.int32).", PyArray_ITEMSIZE(intarray), sizeof(int));
    PyErr_SetString(PyExc_RuntimeError, errorstring);
  } 
}

/* This helper function is from PULP software written by
   Jean-Sébastien Roy 
   http://www.jeannot.org/~js/code/index.en.html#PuLP
*/
void setCPLEXerrorstring(CPXENVptr env, int status) {
  char errorstring[1024];
  char pythonerrorstring[2048];

  if (CPXgeterrorstring (env, status, errorstring) != NULL) {
    sprintf(pythonerrorstring, "CPLEX error (%d): %s", status, errorstring);
    PyErr_SetString(PyExc_RuntimeError, pythonerrorstring);
  } else {
    PyErr_SetString(PyExc_RuntimeError, "Generic CPLEX error.");
  }
}
