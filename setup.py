#!/usr/bin/env python

import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "singlecell",
    version = "0.0.1",
    packages=["singlecell"],
    install_requires=["HTSeq"]
)
