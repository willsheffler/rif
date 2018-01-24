try: from pyrosetta import *
except: pass
from concurrent.futures import *
import rif
import time


def make_pose_huge(p, how_huge=6):
    for i in range(how_huge):
        p = rif.rcl.concatenate_pose(p, p)
    return p


def expensive(p):
    # p = rif.data.poselib.c1
    huge = make_pose_huge(p)
    # print('huge pose is huge:', len(huge))
    sfxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(
        'score0')
    # return 0
    return sfxn(huge)


def main():
    rif.rcl.init_check(strict=False)
    expensive(rif.data.poselib.c1)  # initialize whatever

    poses = [rif.data.poselib.c1] * 64
    # poses = [None] * 64

    with ThreadPoolExecutor(1) as pool:
        t = time.time()
        scores = list(pool.map(expensive, poses))
        print('time 1 thread:', time.time() - t)
    print(scores[:4])

    with ThreadPoolExecutor(8) as pool:
        t = time.time()
        scores = list(pool.map(expensive, poses))
        print('time 8 thread:', time.time() - t)
    print(scores[:4])

    with ProcessPoolExecutor(8) as pool:
        t = time.time()
        scores = list(pool.map(expensive, poses))
        print('time 8 procs:', time.time() - t)
    print(scores[:4])

if __name__ == '__main__':
    try:
        import pyrosetta
        main()
    except ImportError:
        print('no pyrosetta')
