from rif.visualize.vispymol import launch_pymol


def test_launch_pymol():
    def foo():
        pass

    # todo why no -qei?????!?!?!?
    launch_pymol(foo)
