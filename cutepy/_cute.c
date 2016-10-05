#include <stdio.h>
#include "Python.h"
#DEFINE NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "cute.h"


static char module_docstring[] =
    "This module is a wrapper around the correlation functions of CUTE, \
        the code for which can be found at github.com/damonge/CUTE. \
        It has been altered to accept numpy arrays instead of requiring \
        a parameter file. \
        It also includes some extra functions for computing cross-\
        correlations between datasets.";
static char xcorr_docstring[] =
    "Computes the cross-correlation function of two data sets";


static PyObject *xcorr(PyObject *dummy, PyObject *args);


static PyMethodDef mymethods[] = {
    { "xcorr", xcorr,
      METH_VARARGS,
      xcorr_docstring},
    {NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC init_cute(void)
{
    PyObject *m = Py_InitModule("_cute", mymethods, module_docstring);

    if (m == NULL)
        return;

    import_array();
}


static PyObject *xcorr(PyObject *dummy, PyObject *args)
{

    double *ra1, *ra2, *dec1, *dec2, *z1, *z2, *w1, *w2;
    double *randra1, *randra2, *randec1, *randec2;
    double *ranz1, *ranz2, *ranw1, *ranw2;

    int ngals1, ngals2, nrands1, nrands2;
    int pi_max, r_p_nbins, ndecades;
    double r_p_max, O_M, O_L;
    double **result;

    PyArray_Descr *dtype;
    PyObject *data1, *data2, *rand1, *rand2;

    if (!PyArg_ParseTuple(args, "OOOOdiiidd&",
                          &data1, &data2, &rand1, &rand2,
                          &r_p_max, &pi_max, &r_p_nbins, &ndecades,
                          &O_M, &O_L))
        // expected input is four 2D numpy arrays with columns RA, Dec, z, w
        // for both data sets, and both randoms
        // along with maximum projected radius of pairs,
        // maximum line of sight separation between pairs,
        // number of projected radial bins, number of bins per decade,
        // and finally cosmological parameters Om and OL
        return NULL;

    // Now to grab all the arrays and get the galaxy numbers
    // and pass the columns to their designated pointers.
    PyArrayObject *dat1arr, *dat2arr, *ran1arr, *ran2arr;
    ngals1 = PyArray_SHAPE(data1)[0];
    ra1 = (double *)
    dat1arr = (PyArrayObject *) PyArray_FROM_OTF(data1,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    ngals2 = PyArray_SHAPE(data2)[0];
    dat2arr = (PyArrayObject *) PyArray_FROM_OTF(data2,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    nrands1 = PyArray_SHAPE(rand1)[0];
    ran1arr = (PyArrayObject *) PyArray_FROM_OTF(rand1,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    nrands2 = PyArray_SHAPE(rand2)[0];
    ran2arr = (PyArrayObject *) PyArray_FROM_OTF(rand2,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);

    if(dat1arr == NULL | dat2arr == NULL | ran1arr == NULL | ran2arr == NULL) {
        Py_XDECREF(dat1arr);
        Py_XDECREF(dat2arr);
        Py_XDECREF(ran1arr);
        Py_XDECREF(ran2arr);
        return NULL;
    }

    ra1 = (double *)malloc(ngals1 * sizeof(double));
    dec1 = (double *)malloc(ngals1 * sizeof(double));
    z1 = (double *)malloc(ngals1 * sizeof(double));
    w1 = (double *)malloc(ngals1 * sizeof(double));
    int k = 0;
    for (i = 0; i < 4 * ngals1; i++) {
        if (i % 4 == 0)
            ra1[k] = dat1arr[i];
        else if (i % 4 == 1)
            dec1[k] = dat1arr[i];
        else if (i % 4 == 2)
            z1[k] = dat1arr[i];
        else {
            w1[k] = dat1arr[i];
            k++;
        }
    }

    ra2 = (double *)malloc(ngals2 * sizeof(double));
    dec2 = (double *)malloc(ngals2 * sizeof(double));
    z2 = (double *)malloc(ngals2 * sizeof(double));
    w2 = (double *)malloc(ngals2 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * ngals2; i++) {
        if (i % 4 == 0)
            ra2[k] = dat2arr[i];
        else if (i % 4 == 1)
            dec2[k] = dat2arr[i];
        else if (i % 4 == 2)
            z2[k] = dat2arr[i];
        else {
            w2[k] = dat2arr[i];
            k++;
        }
    }

    randra1 = (double *)malloc(nrands1 * sizeof(double));
    randec1 = (double *)malloc(nrands1 * sizeof(double));
    ranz1 = (double *)malloc(nrands1 * sizeof(double));
    ranw1 = (double *)malloc(nrands1 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * nrands1; i++) {
        if (i % 4 == 0)
            randra1[k] = ran1arr[i];
        else if (i % 4 == 1)
            randec1[k] = ran1arr[i];
        else if (i % 4 == 2)
            ranz1[k] = ran1arr[i];
        else {
            ranw1[k] = ran1arr[i];
            k++;
        }
    }

    randra2 = (double *)malloc(nrands2 * sizeof(double));
    randec2 = (double *)malloc(nrands2 * sizeof(double));
    ranz2 = (double *)malloc(nrands2 * sizeof(double));
    ranw2 = (double *)malloc(nrands2 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * nrands2; i++) {
        if (i % 4 == 0)
            randra2[k] = ran2arr[i];
        else if (i % 4 == 1)
            randec2[k] = ran2arr[i];
        else if (i % 4 == 2)
            ranz2[k] = ran2arr[i];
        else {
            ranw2[k] = ran2arr[i];
            k++;
        }
    }


    // And let's run this motha
    result = runfull_xcorr(ngals1, ra1, dec1, z1, w1,
                           ngals2, ra2, dec2, z2, w2,
                           nrands1, randra1, randec1, ranz1, ranw1,
                           nrands2, randra2, randec2, ranz2, ranw2,
                           r_p_max, pi_max, r_p_nbins, ndecades,
                           O_M, O_L);

}
