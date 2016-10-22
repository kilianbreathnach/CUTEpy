#include <stdio.h>
#include "Python.h"
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
static char autopairs_docstring[] =
    "Computes the pair counts of a dataset for the given scales";
static char xpairs_docstring[] =
    "Computes the cross pair counts of two  datasets for the given scales";


static PyObject *xcorr(PyObject *dummy, PyObject *args);
static PyObject *autopairs(PyObject *dummy, PyObject *args);
static PyObject *xpairs(PyObject *dummy, PyObject *args);

static PyMethodDef mymethods[] = {
    { "xcorr", xcorr,
      METH_VARARGS,
      xcorr_docstring},
    {"autopairs", autopairs,
      METH_VARARGS,
      autopairs_docstring},
    {"xpairs", xpairs,
      METH_VARARGS,
      xpairs_docstring},
    {NULL, NULL, 0, NULL} /* Sentinel */
};


PyMODINIT_FUNC init_cute(void)
{
    PyObject *m = Py_InitModule("_cute", mymethods);

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

    PyObject *data1, *data2, *rand1, *rand2, *logspace;
    npy_intp *shapeX;

    if (!PyArg_ParseTuple(args, "OOOOdiiiOdd",
                          &data1, &data2, &rand1, &rand2,
                          &r_p_max, &pi_max, &r_p_nbins, &ndecades,
                          &logspace, &O_M, &O_L))
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
    dat1arr = (PyArrayObject *) PyArray_FROM_OTF(data1,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(dat1arr);
    ngals1 = shapeX[0];

    dat2arr = (PyArrayObject *) PyArray_FROM_OTF(data2,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(dat2arr);
    ngals2 = shapeX[0];

    ran1arr = (PyArrayObject *) PyArray_FROM_OTF(rand1,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(ran1arr);
    nrands1 = shapeX[0];

    ran2arr = (PyArrayObject *) PyArray_FROM_OTF(rand2,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(ran2arr);
    nrands2 = shapeX[0];

    if((dat1arr == NULL) | (dat2arr == NULL) | (ran1arr == NULL) | (ran2arr == NULL)) {
        Py_XDECREF(dat1arr);
        Py_XDECREF(dat2arr);
        Py_XDECREF(ran1arr);
        Py_XDECREF(ran2arr);
        return NULL;
    }

    double *dat1_data = (double *) PyArray_DATA(dat1arr);
    ra1 = (double *)malloc(ngals1 * sizeof(double));
    dec1 = (double *)malloc(ngals1 * sizeof(double));
    z1 = (double *)malloc(ngals1 * sizeof(double));
    w1 = (double *)malloc(ngals1 * sizeof(double));
    int i;
    int k = 0;
    for (i = 0; i < 4 * ngals1; i++) {
        if (i % 4 == 0)
            ra1[k] = dat1_data[i];
        else if (i % 4 == 1)
            dec1[k] = dat1_data[i];
        else if (i % 4 == 2)
            z1[k] = dat1_data[i];
        else {
            w1[k] = dat1_data[i];
            k++;
        }
    }

    double *dat2_data = (double *) PyArray_DATA(dat2arr);
    ra2 = (double *)malloc(ngals2 * sizeof(double));
    dec2 = (double *)malloc(ngals2 * sizeof(double));
    z2 = (double *)malloc(ngals2 * sizeof(double));
    w2 = (double *)malloc(ngals2 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * ngals2; i++) {
        if (i % 4 == 0)
            ra2[k] = dat2_data[i];
        else if (i % 4 == 1)
            dec2[k] = dat2_data[i];
        else if (i % 4 == 2)
            z2[k] = dat2_data[i];
        else {
            w2[k] = dat2_data[i];
            k++;
        }
    }

    double *ran1_data = (double *) PyArray_DATA(ran1arr);
    randra1 = (double *)malloc(nrands1 * sizeof(double));
    randec1 = (double *)malloc(nrands1 * sizeof(double));
    ranz1 = (double *)malloc(nrands1 * sizeof(double));
    ranw1 = (double *)malloc(nrands1 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * nrands1; i++) {
        if (i % 4 == 0)
            randra1[k] = ran1_data[i];
        else if (i % 4 == 1)
            randec1[k] = ran1_data[i];
        else if (i % 4 == 2)
            ranz1[k] = ran1_data[i];
        else {
            ranw1[k] = ran1_data[i];
            k++;
        }
    }

    double *ran2_data = (double *) PyArray_DATA(ran2arr);
    randra2 = (double *)malloc(nrands2 * sizeof(double));
    randec2 = (double *)malloc(nrands2 * sizeof(double));
    ranz2 = (double *)malloc(nrands2 * sizeof(double));
    ranw2 = (double *)malloc(nrands2 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * nrands2; i++) {
        if (i % 4 == 0)
            randra2[k] = ran2_data[i];
        else if (i % 4 == 1)
            randec2[k] = ran2_data[i];
        else if (i % 4 == 2)
            ranz2[k] = ran2_data[i];
        else {
            ranw2[k] = ran2_data[i];
            k++;
        }
    }

    double *result = (double *)malloc(r_p_nbins * pi_max * 4 * sizeof(double));

    int logtrue = PyObject_IsTrue(logspace);

    // And let's run this motha
    runfull_xcorr(result, ngals1, ra1, dec1, z1, w1,
                          ngals2, ra2, dec2, z2, w2,
                          nrands1, randra1, randec1, ranz1, ranw1,
                          nrands2, randra2, randec2, ranz2, ranw2,
                          r_p_max, pi_max, r_p_nbins, ndecades,
                          logtrue, O_M, O_L);

    Py_DECREF(dat1arr);
    Py_DECREF(dat2arr);
    Py_DECREF(ran1arr);
    Py_DECREF(ran2arr);

    npy_intp dims[2];
    dims[0] = r_p_nbins * pi_max;
    dims[1] = 4;
    PyObject *retobj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE,
                                                 result);
    return retobj;
}


static PyObject *autopairs(PyObject *dummy, PyObject *args)
{
    double *ra, *dec, *z, *w;

    int ngals;
    int pi_max, r_p_nbins, ndecades;
    double r_p_max, O_M, O_L;

    PyObject *data, *logspace;
    npy_intp *shapeX;

    if (!PyArg_ParseTuple(args, "OdiiiOdd",
                          &data,
                          &r_p_max, &pi_max, &r_p_nbins, &ndecades,
                          &logspace, &O_M, &O_L))
        // expected input is four 2D numpy arrays with columns RA, Dec, z, w
        // for both data sets, and both randoms
        // along with maximum projected radius of pairs,
        // maximum line of sight separation between pairs,
        // number of projected radial bins, number of bins per decade,
        // and finally cosmological parameters Om and OL
        return NULL;

    PyArrayObject *datarr;
    datarr = (PyArrayObject *) PyArray_FROM_OTF(data,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(datarr);
    ngals = shapeX[0];

    if(datarr == NULL) {
        Py_XDECREF(datarr);
        return NULL;
    }

    double *dat_data = (double *) PyArray_DATA(datarr);
    ra = (double *)malloc(ngals * sizeof(double));
    dec = (double *)malloc(ngals * sizeof(double));
    z = (double *)malloc(ngals * sizeof(double));
    w = (double *)malloc(ngals * sizeof(double));
    int i;
    int k = 0;
    for (i = 0; i < 4 * ngals; i++) {
        if (i % 4 == 0)
            ra[k] = dat_data[i];
        else if (i % 4 == 1)
            dec[k] = dat_data[i];
        else if (i % 4 == 2)
            z[k] = dat_data[i];
        else {
            w[k] = dat_data[i];
            k++;
        }
    }

    double *result = (double *)malloc(r_p_nbins * pi_max * 3 * sizeof(double));

    int logtrue = PyObject_IsTrue(logspace);

    // And let's run this motha
    run_auto_3d_ps(result, ngals, ra, dec, z, w,
                           r_p_max, pi_max, r_p_nbins, ndecades,
                           logtrue, O_M, O_L);

    Py_DECREF(datarr);

    npy_intp dims[2];
    dims[0] = r_p_nbins * pi_max;
    dims[1] = 3;
    PyObject *retobj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE,
                                                 result);
    return retobj;
}


static PyObject *xpairs(PyObject *dummy, PyObject *args)
{
    double *ra1, *ra2, *dec1, *dec2, *z1, *z2, *w1, *w2;

    int ngals1, ngals2;
    int pi_max, r_p_nbins, ndecades;
    double r_p_max, O_M, O_L;

    PyObject *data1, *data2, *logspace;
    npy_intp *shapeX;

    if (!PyArg_ParseTuple(args, "OOdiiiOdd",
                          &data1, &data2,
                          &r_p_max, &pi_max, &r_p_nbins, &ndecades,
                          &logspace, &O_M, &O_L))
        // expected input is four 2D numpy arrays with columns RA, Dec, z, w
        // for both data sets, and both randoms
        // along with maximum projected radius of pairs,
        // maximum line of sight separation between pairs,
        // number of projected radial bins, number of bins per decade,
        // and finally cosmological parameters Om and OL
        return NULL;

    // Now to grab all the arrays and get the galaxy numbers
    // and pass the columns to their designated pointers.
    PyArrayObject *dat1arr, *dat2arr;
    dat1arr = (PyArrayObject *) PyArray_FROM_OTF(data1,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(dat1arr);
    ngals1 = shapeX[0];

    dat2arr = (PyArrayObject *) PyArray_FROM_OTF(data2,
                                                 NPY_DOUBLE,
                                                 NPY_ARRAY_IN_ARRAY);
    shapeX = PyArray_SHAPE(dat2arr);
    ngals2 = shapeX[0];

    if((dat1arr == NULL) | (dat2arr == NULL)) {
        Py_XDECREF(dat1arr);
        Py_XDECREF(dat2arr);
        return NULL;
    }

    double *dat1_data = (double *) PyArray_DATA(dat1arr);
    ra1 = (double *)malloc(ngals1 * sizeof(double));
    dec1 = (double *)malloc(ngals1 * sizeof(double));
    z1 = (double *)malloc(ngals1 * sizeof(double));
    w1 = (double *)malloc(ngals1 * sizeof(double));
    int i;
    int k = 0;
    for (i = 0; i < 4 * ngals1; i++) {
        if (i % 4 == 0)
            ra1[k] = dat1_data[i];
        else if (i % 4 == 1)
            dec1[k] = dat1_data[i];
        else if (i % 4 == 2)
            z1[k] = dat1_data[i];
        else {
            w1[k] = dat1_data[i];
            k++;
        }
    }

    double *dat2_data = (double *) PyArray_DATA(dat2arr);
    ra2 = (double *)malloc(ngals2 * sizeof(double));
    dec2 = (double *)malloc(ngals2 * sizeof(double));
    z2 = (double *)malloc(ngals2 * sizeof(double));
    w2 = (double *)malloc(ngals2 * sizeof(double));
    k = 0;
    for (i = 0; i < 4 * ngals2; i++) {
        if (i % 4 == 0)
            ra2[k] = dat2_data[i];
        else if (i % 4 == 1)
            dec2[k] = dat2_data[i];
        else if (i % 4 == 2)
            z2[k] = dat2_data[i];
        else {
            w2[k] = dat2_data[i];
            k++;
        }
    }

    double *result = (double *)malloc(r_p_nbins * pi_max * 3 * sizeof(double));

    int logtrue = PyObject_IsTrue(logspace);
    printf("logtrue = %d\n", logtrue);

    // And let's run this motha
    run_cross_3d_ps(result, ngals1, ra1, dec1, z1, w1,
                            ngals2, ra2, dec2, z2, w2,
                            r_p_max, pi_max, r_p_nbins, ndecades,
                            logtrue, O_M, O_L);

    Py_DECREF(dat1arr);
    Py_DECREF(dat2arr);

    npy_intp dims[2];
    dims[0] = r_p_nbins * pi_max;
    dims[1] = 3;
    PyObject *retobj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE,
                                                 result);
    return retobj;
}
