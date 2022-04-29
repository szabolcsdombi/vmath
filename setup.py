import sys

from setuptools import Extension, setup

extra_compile_args = []

if sys.platform.startswith('linux'):
    extra_compile_args = ['-fpermissive', '-Wno-write-strings', '-Wno-narrowing']

if sys.platform.startswith('darwin'):
    extra_compile_args = ['-std=c++11', '-Wno-writable-strings', '-Wno-c++11-narrowing']

ext = Extension(
    name='vmath',
    sources=['./vmath.cpp'],
    define_macros=[('PY_SSIZE_T_CLEAN', None)],
    extra_compile_args=extra_compile_args,
)

with open('README.md') as readme:
    long_description = readme.read()

setup(
    name='vmath',
    version='0.3.0',
    ext_modules=[ext],
    data_files=[('.', ['vmath.pyi'])],
    license='MIT',
    python_requires='>=3.6',
    platforms=['any'],
    description='Compact Python OpenGL rendering library',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Szabolcs Dombi',
    author_email='cprogrammer1994@gmail.com',
    url='https://github.com/szabolcsdombi/vmath/',
    project_urls={
        'Documentation': 'https://vmath.readthedocs.io/',
        'Source': 'https://github.com/szabolcsdombi/vmath/',
        'Bug Tracker': 'https://github.com/szabolcsdombi/vmath/issues/',
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Games/Entertainment',
        'Topic :: Multimedia :: Graphics',
        'Topic :: Multimedia :: Graphics :: 3D Rendering',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    keywords=[
        'vector',
        'matrix',
        'quaternion',
    ],
)
