import __main__

try:
    import pymol
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False


def launch_pymol(doit, args='-qei'):
    # pymol doesn't seem to get these args
    __main__.pymol_argv = ['pymol', args]
    pymol.finish_launching()
    doit()


if __name__ == '__main__':
    launch_pymol(lambda x: x)
