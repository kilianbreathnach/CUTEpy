from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    ext_modules=[Extension("_cute", ["cutepy/_cute.c",
                                            "cutepy/cute.c",
                                            "cutepy/src/define.c",
                                            "cutepy/src/common.c",
                                            "cutepy/src/cosmo.c",
                                            "cutepy/src/boxes3D.c",
                                            "cutepy/src/correlator.c"],
                            include_dirs=["cutepy",
                                          "cutepy/src",
                                          "/data1/kww231/home/lib/python2.7/site-packages/numpy/core/include",
                                          "/data1/kww231/home/include/gsl"],
                            define_macros=[("_VERBOSE", None),
                                           ("_WITH_WEIGHTS", None),
                                           ("_HISTO_2D_64", None)],
                            library_dirs=["/data1/kww231/home/lib"],
                            libraries=["m", "gsl", "gslcblas", "gomp"],
                            extra_compile_args=["-Wall",
                                                "-O3",
                                                "-fopenmp"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
