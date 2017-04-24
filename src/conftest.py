from __future__ import print_function
import sys


def pytest_namespace():
    print("Do not run pytest in rif project src dir!")
    print("   location:", __file__)
    sys.exit(-1)  # Do not run pytest in rif project src dir!
