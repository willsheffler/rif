try:
    import pyrosetta
    from rosetta import *
    import rosetta
except ImportError:
    print('no pyrosetta!')


def main():
    pyrosetta.init('-beta -include_patches patches/VirtualBB.txt')
    chem_manager = rosetta.core.chemical.ChemicalManager
    rts = chem_manager.get_instance().residue_type_set("fa_standard")
    ggg = rosetta.core.pose.Pose()
    core.pose.make_pose_from_sequence(ggg, 'AAA', rts)
    virt_bb = rosetta.core.chemical.VIRTUAL_BB
    core.pose.remove_lower_terminus_type_from_pose_residue(ggg, 1)
    core.pose.remove_upper_terminus_type_from_pose_residue(ggg, 3)
    core.pose.add_variant_type_to_pose_residue(ggg, virt_bb, 1)


if __name__ == '__main__':
    main()
