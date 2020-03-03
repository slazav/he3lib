#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <../he3.h>

static PyObject *He3Error;


// functions
static PyObject *
he3_tc(PyObject *self, PyObject *args) {

    double d1;
    if (!PyArg_ParseTuple(args, "d", &d1)) return NULL;

    return PyFloat_FromDouble(he3_tc_(&d1));
}


// list of all functions:
static PyMethodDef He3Methods[] = {
    {"he3_tc",  he3_tc, METH_VARARGS, "test command."},
    {NULL, NULL, 0, NULL}
};

// module structure
static struct PyModuleDef he3module = {
    PyModuleDef_HEAD_INIT,
    "he3",    /* name of module */
    NULL,     /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    He3Methods
};

PyMODINIT_FUNC
PyInit_he3lib(void){
    PyObject *m;

    m = PyModule_Create(&he3module);
    if (m == NULL)
        return NULL;

    He3Error = PyErr_NewException("he3.error", NULL, NULL);
    Py_XINCREF(He3Error);
    if (PyModule_AddObject(m, "error", He3Error) < 0) {
        Py_XDECREF(He3Error);
        Py_CLEAR(He3Error);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}


