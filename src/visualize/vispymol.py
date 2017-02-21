import __main__
import pymol


def launch_pymol(doit):
    if __name__ == '__main__':
        __main__.pymol_argv = ['pymol', '-qei']
        pymol.finish_launching()
        doit()
