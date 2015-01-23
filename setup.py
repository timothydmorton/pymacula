#!/usr/bin/env python

import os
import sys
from numpy.distutils.core import setup, Extension

# Hackishly inject a constant into builtins to enable importing of the
# package even if numpy isn't installed. Only do this if we're not
# running the tests!
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__PYMACULA_SETUP__ = True
import pymacula
version = pymacula.__version__

# Publish the library to PyPI.
if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()

# Push a new tag to GitHub.
if "tag" in sys.argv:
    os.system("git tag -a {0} -m 'version {0}'".format(version))
    os.system("git push --tags")
    sys.exit()

# Set up the compiled extension.
sources = [os.path.join('pymacula','macula.f90')]
extensions = [Extension("pymacula._macula", sources=sources)]

setup(
    name="pymacula",
    version=version,
    author="Timothy D. Morton",
    author_email="tim.morton@gmail.com",
    url="https://github.com/timothydmorton/pymacula",
    license="MIT",
    packages=["pymacula"],
    ext_modules=extensions,
    description="Starspot modeling",
    long_description=open("README.rst").read(),
    #package_data={"": ["README.rst", "LICENSE"]},
    #include_package_data=True,
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
