import tempfile
from collections import defaultdict

def is_rosetta_pose(thing):
    try:
        import rosetta
        return isinstance(thing,rosetta.core.pose.Pose)
    except ImportError:
        return False


def pymol_load(things):
    from pymol import cmd
    if not isinstance(things,list):
        things = [things]
    tmpdir = tempfile.mkdtemp()
    seenit = defaultdict(lambda: -1)
    for thing in things:
        name = 'rif_thing'
        seenit[name] += 1
        name += '_%i' % seenit[name]
        if is_rosetta_pose(thing):
            fname = tmpdir + '/' + name + '.pdb'
            thing.dump_pdb(fname)
            cmd.load(fname)
        else:
            raise NotImplementedError("don't know how to show "+ str(type(thing)))


def showme_pymol(what, headless=False):
    import pymol
    if headless:
        pymol.pymol_argv = ['pymol','-c']
    pymol.finish_launching()
    from pymol import cmd
    pymol_load(what)
    if headless:
        cmd.quit()


def showme(*args, how='pymol', **kwargs):
    if how == 'pymol':
        showme_pymol(*args, **kwargs)
    else:
        raise NotImplementedError('showme how="%s" not implemented'%how)

