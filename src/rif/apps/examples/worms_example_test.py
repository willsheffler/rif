from rif.worm import *
from rif.data import poselib
from rif.vis import showme
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from time import perf_counter
import sys


def doit(nseg, nworker):
    # nsplice = SpliceSite(sele=[':5', ], polarity='N')
    # csplice = SpliceSite(sele=['-5:', ], polarity='C')
    helix = Spliceable(poselib.curved_helix, sites=[(1, 'N'), ('-4:', 'C')])

    # strand = Spliceable(strand_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # loop = Spliceable(loop_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # splicables = [helix, strand, loop]
    segments = ([Segment([helix], exit='C'), ]
                + [Segment([helix], entry='N', exit='C')] * (nseg - 2)
                + [Segment([helix], entry='N')])
    t = perf_counter()
    worms = grow(segments,
                 SegmentSym('C1', lever=20),
                 thresh=10, max_workers=nworker,
                 executor=ProcessPoolExecutor)
    t = perf_counter() - t
    count = np.prod([len(s) for s in segments])
    s = worms.scores

    try: ptile = np.percentile(s, range(0, 100, 20))
    except: ptile = []
    print('quantile', ptile)
    print('best10  ', s[:10])
    print('nseg %2i' % nseg,
          'best %7.3f' % (s[0] if len(s) else 999),
          'tgrow %7.2f' % t,
          'rate %7.3fM/s' % (count / t / 1000000),
          'npass %8i' % len(s))
    sys.stdout.flush()

    for i in range(min(10, len(worms))):
        worms.pose(i, onechain=True).dump_pdb('worm_nseg%i_%i.pdb' % (nseg, i))


def main():
    try: Nseg0 = int(sys.argv[-2])
    except: Nseg0 = 13
    try: Nseg1 = int(sys.argv[-1])
    except: Nseg1 = 13
    for n in range(Nseg0, Nseg1 + 1):
        fail = True
        while fail:
            try:
                doit(n, 0)
                fail = False
            except BrokenProcessPool as e:
                print('failed run on ' + str(n))
                print(e)

if __name__ == '__main__':
    main()
