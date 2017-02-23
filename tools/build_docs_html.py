#!/usr/bin/env python

import os

os.sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from build_utils import make_docs


if __name__ == '__main__':
    make_docs('html')
