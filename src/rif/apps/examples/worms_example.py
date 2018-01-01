from rif.worm import *
from rif.data import poselib
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from time import perf_counter


def main():
    try: Nseg = int(sys.argv[-2])
    except: Nseg = 12
    try: Nworker = int(sys.argv[-1])
    except: Nworker = 0

    # nsplice = SpliceSite(sele=[':5', ], polarity='N')
    # csplice = SpliceSite(sele=['-5:', ], polarity='C')
    helix = Spliceable(poselib.curved_helix, sites=[(1, 'N'), ('-4:', 'C')])

    # strand = Spliceable(strand_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # loop = Spliceable(loop_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # splicables = [helix, strand, loop]
    segments = ([Segment([helix], exit='C'), ]
                + [Segment([helix], entry='N', exit='C')] * (Nseg - 2)
                + [Segment([helix], entry='N')])
    t = perf_counter()
    worms = grow(segments,
                 SegmentXform('C1', lever=20),
                 thresh=5, max_workers=Nworker,
                 executor=ProcessPoolExecutor)
    t = perf_counter() - t
    s = worms.scores

    try: ptile = np.percentile(s, range(0, 100, 20))
    except: ptile = []
    print('ptile', ptile)
    print('best10', s[:10])
    print('nseg', Nseg,
          'npass', len(s),
          'best', s[0] if s else 999,
          'tgrow', t)

if __name__ == '__main__':
    main()
