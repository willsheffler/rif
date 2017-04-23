from __future__ import print_function
import pytest
import os
import sys
import rif
from rif.data.resource import resource_path


def test_resource_path():
    with pytest.raises(IOError):
        os.path.exists(resource_path("I do not exist"))
