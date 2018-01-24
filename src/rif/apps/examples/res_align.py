try: from pyrosetta import *
except: print('No pyrosetta!')
import numpy as np
from rif import rcl, vis  # 'rosetta compatibility layer'
from rif.data import poselib  # just loads pdbs form src/rif/data/pdb


def main():
    pose1 = poselib.c1
    pose2 = poselib.c2

    stubs1 = rcl.bbstubs(pose1, [len(pose1)])['raw']  # gets 'stub' for reslist
    stubs2 = rcl.bbstubs(pose2, [1])['raw']  # raw field is position matrix
    print(stubs1.shape)  # homo coords numpy array n x 4 x 4
    print(stubs2.shape)
    xalign = stubs1 @ np.linalg.inv(stubs2)  # a @ b is np.matmul(a, b)
    print(xalign.shape)
    rcl.xform_pose(xalign[0], pose2)

    # append residues 2-7 of pose2 to pose1
    # note: O and H positions may be fucked up...
    # not sure easiest way to fix... cart min?
    rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(
        pose1, len(pose1))
    rosetta.core.pose.append_subpose_to_pose(pose1, pose2, 2, 7)
    # vis.showme(pose1)


if __name__ == '__main__':
    try:
        import pyrosetta
        main()
    except ImportError: pass
