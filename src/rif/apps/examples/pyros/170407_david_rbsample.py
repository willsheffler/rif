try:
    import pyrosetta
    from rosetta import *
    import rosetta
except ImportError:
    print('no pyrosetta!')


def main():
    pyrosetta.init('-beta')

    # making xforms
    for iter in range(10):
        x = rosetta.numeric.random.gaussian_random_xform(10.0, 0.5)
        print("should be 10-ish", x.rotation_angle_degrees())
        print('should be 0.5-ish', x.t.length())

    # customizing a score function
    sf = pyrosetta.get_score_function()
    # do print(sf) and look for weights to see what you can modify
    sf.set_weight(core.scoring.fa_dun, 0.0)

    # use xform to move pose around
    chem_manager = rosetta.core.chemical.ChemicalManager
    rts = chem_manager.get_instance().residue_type_set("fa_standard")

    gly = core.pose.Pose()
    asn = core.pose.Pose()
    core.pose.make_pose_from_sequence(gly, 'G', rts)
    core.pose.make_pose_from_sequence(asn, 'N', rts)
    gn = gly
    # anchor to 1st residue, rot round CG of res 2
    gn.append_pose_by_jump(asn, 1, 'CA', 'CG')
    j = gn.jump(1)  # only one jump in this system
    j.set_translation(j.get_translation() +
                      numeric.xyzVector_double_t(5, 0, 0))
    gn.set_jump(1, j)  # move apart a bit

    for i in range(10):
        j = gn.jump(1)
        # this is better than below, probably
        print(j.gaussian_move(1, 0.0, 5.0))
        # x = rosetta.numeric.random.gaussian_random_xform(5.0, 0.25)
        # j.set_rotation(j.get_rotation()x.R)
        # j.set_translation(j.get_translation() + x.t)
        gn.set_jump(1, j)
        gn.dump_pdb('test_%i.pdb' % i)


if __name__ == '__main__':
    main()
