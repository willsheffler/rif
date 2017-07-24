from builtins import *

from pyrosetta import *
from rosetta.numeric import xyzVector_double_t as Vec
from rosetta.protocols.sic_dock import trans_pose

nfold = 3

# temporary workaround, fixed in newer pyrosetta


def mat_vec_mul(M, V):
    return rosetta.numeric.xyzVector_double_t(
        M.xx() * V.x + M.xy() * V.y + M.xz() * V.z,
        M.yx() * V.x + M.yy() * V.y + M.yz() * V.z,
        M.zx() * V.x + M.zy() * V.y + M.zz() * V.z
    )


def append_pose_to_pose(pose1, pose2, start_res=1, end_res=None, new_chain=True):
    if end_res is None:
        end_res = pose2.total_residue()
    pose1.append_residue_by_jump(pose2.residue(
        start_res), pose1.total_residue(), "", "", new_chain)
    for i in range(start_res + 1, end_res + 1):
        if pose2.residue(i).is_lower_terminus():
            if i > 1 and pose2.chain(i) == pose2.chain(i - 1):
                pose1.append_residue_by_jump(pose2.residue(
                    i), pose1.total_residue(), "", "", False)
            else:
                pose1.append_residue_by_jump(pose2.residue(
                    i), pose1.total_residue(), "", "", True)
        else:
            pose1.append_residue_by_bond(pose2.residue(i))


def rot_pose(pose, xform):
    for ir in range(1, pose.total_residue() + 1):
        for ia in range(1, pose.residue_type(ir).natoms() + 1):
            aid = AtomID(ia, ir)
            pose.set_xyz(aid, mat_vec_mul(xform, pose.xyz(aid)))


# def trans_pose(pose, trans):
    # for ir in range(1, pose.total_residue()+1):
    #     for ia in range(1, pose.residue_type(ir).natoms()+1):
    #         aid = AtomID(ia,ir)
    #         pose.set_xyz(aid, pose.xyz(aid) + trans)


def make_cyclic(pose, nfold):
    rot = rosetta.numeric.z_rotation_matrix_degrees_double_t(360.0 / nfold)
    cyclic = Pose()
    for i in range(nfold):
        rot_pose(cyclic, rot)
        append_pose_to_pose(cyclic, pose)
    return cyclic


def center_pose_z(pose):
    cenz = rosetta.core.pose.get_center_of_mass(pose).z
    trans_pose(pose, rosetta.numeric.xyzVector_double_t(0, 0, -cenz))


class AxleDock(object):
    def __init__(self, pose1, pose2, nfold, resl=5.0):
        self.pose1 = pose1
        self.pose2 = pose2
        self.nfold = nfold
        self.resl = resl
        self.cyclic2 = make_cyclic(self.pose2, self.nfold)
        center_pose_z(self.pose1)
        center_pose_z(self.cyclic2)
        self.sfxn = ScoreFunction()
        self.sfxn.set_weight(rosetta.core.scoring.fa_rep, 1.0)
        self.pose = Pose(pose1)
        append_pose_to_pose(self.pose, self.cyclic2)
        self.score0 = self.sfxn(self.pose1) + self.sfxn(self.cyclic2)
        self.offset = 0

    def reset(self):
        trans_pose(self.pose, Vec(0, 0, -self.offset),
                   1, self.pose1.total_residue())
        self.offset = 0

    def wheel_vs_wheel_score(self):
        return self.sfxn(self.pose) - self.score0

    def slide_out_of_contact(self, direction):
        self.reset()
        assert direction in 'up down'
        direction = 1.0 if direction is 'up' else -1.0
        delta = direction * 256.0
        while abs(delta) >= 0.1:
            self.offset += delta
            trans_pose(self.pose, Vec(0, 0, delta),
                       1, self.pose1.total_residue())
            score = self.wheel_vs_wheel_score()
            in_contact = abs(score) > 0.001
            moving_away = delta * direction > 0.0
            # print(delta, score, in_contact, moving_away)
            if in_contact != moving_away:
                delta *= -1.0
            delta /= 2.0

if __name__ == '__main__':
    init('-corrections:beta_nov16')
    p1 = pose_from_file("/Users/sheffler/data/C3/C3_1woz_1.pdb.gz")
    p2 = pose_from_file("/Users/sheffler/data/C3/C3_3l8r_1.pdb.gz")
    # cyclic1.dump_pdb("/Users/sheffler/Desktop/cyclic1.pdb")
    # cyclic2.dump_pdb("/Users/sheffler/Desktop/cyclic2.pdb")
    docker = AxleDock(p1, p2, 3)
    assert abs(docker.offset) < 0.001
    docker.slide_out_of_contact('up')
    docker.pose.dump_pdb('/Users/sheffler/Desktop/up1.pdb')
    docker.slide_out_of_contact('down')
    docker.pose.dump_pdb('/Users/sheffler/Desktop/down1.pdb')
    docker.slide_out_of_contact('up')
    docker.pose.dump_pdb('/Users/sheffler/Desktop/up2.pdb')
