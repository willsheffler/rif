import sys
import os


def main():
    if 'CXX' not in os.environ:
        cxx = 'DEFAULT_CXX'
    else:
        cxx = os.environ['CXX']
    print(cxx, sys.executable, sys.version)

if __name__ == '__main__':
    main()
