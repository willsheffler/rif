import tempfile
from collections import defaultdict


def is_rosetta_pose(to_show):
    try:
        import rosetta
        return isinstance(to_show, rosetta.core.pose.Pose)
    except ImportError:
        return False


def pymol_load_pose(pose, name):
    from pymol import cmd
    tmpdir = tempfile.mkdtemp()
    fname = tmpdir + '/' + name + '.pdb'
    pose.dump_pdb(fname)
    cmd.load(fname)


def pymol_xform(name, xform):
    from pymol import cmd
    assert name in cmd.get_object_list()
    cmd.transform_object(name, xform.flatten())


def pymol_load(to_show, state=None):
    if state is None:
        state = dict(seenit=defaultdict(lambda: -1))
    if isinstance(to_show, list):
        for t in to_show:
            state = pymol_load(t, state)
    elif isinstance(to_show, dict):
        assert 'pose' in to_show
        state = pymol_load(to_show['pose'], state)
        pymol_xform(to_show['position'], state['last_obj'])
    elif is_rosetta_pose(to_show):
        name = 'rif_thing'
        state['seenit'][name] += 1
        name += '_%i' % state['seenit'][name]
        pymol_load_pose(to_show, name)
        state['last_obj'] = name
    else:
        raise NotImplementedError(
            "don't know how to show " + str(type(to_show)))
    return state


def showme_pymol(what, headless=False, block=False):
    import pymol
    pymol.pymol_argv = ['pymol']
    if headless:
        pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    from pymol import cmd
    r = pymol_load(what)
    cmd.set('internal_gui_width', '20')
    cmd.do('full')
    cmd.zoom()
    import time
    while block:
        time.sleep(1)
    return r


def showme(*args, how='pymol', **kwargs):
    if how == 'pymol':
        return showme_pymol(*args, **kwargs)
    else:
        raise NotImplementedError('showme how="%s" not implemented' % how)
