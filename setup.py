from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    ext_modules=[Extension("cutepy.xcorr", ["cutepy/_cute.c",
                                            "cutepy/cute.c",
                                            "cutepy/src/define.c",
                                            "cutepy/src/io.c",
                                            "cutepy/src/common.c",
                                            "cutepy/src/cosmo.c",
                                            "cutepy/src/boxes3D.c",
                                            "cutepy/src/correlator.c"],
                            include_dirs=["cutepy",
                                          "cutepy/src",
                                          "/data1/kww231/home/include/gsl"],
                            define_macros=[("_VERBOSE", None),
                                           ("_WITH_WEIGHTS", None),
                                           ("_HISTO_2D_64", None)],
                            library_dirs=["/data1/kww231/home/lib"],
                            libraries=["math", "gsl", "gslcblas"],
                            extra_compile_args=["-Wall",
                                                "-O3",
                                                "-fopenmp"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
