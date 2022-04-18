from setuptools import Extension, setup

ext = Extension(
    name='vmath',
    sources=['./vmath.cpp'],
    define_macros=[('PY_SSIZE_T_CLEAN', None)],
    include_dirs=[],
    library_dirs=[],
    libraries=[],
)

setup(
    name='vmath',
    version='0.2.0',
    ext_modules=[ext],
    data_files=[('.', ['vmath.pyi'])],
)
