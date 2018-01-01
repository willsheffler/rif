from rif.worm import *
from rif.data import poselib


def main():
    # nsplice = SpliceSite(sele=[':5', ], polarity='N')
    # csplice = SpliceSite(sele=['-5:', ], polarity='C')
    helix = Spliceable(poselib.curved_helix, sites=[(1, 'N'), ('-4:', 'C')])

    # strand = Spliceable(strand_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # loop = Spliceable(loop_pose, sites=[(':3', 'N'), ('-3:', 'C')])
    # splicables = [helix, strand, loop]
    segments = ([Segment([helix], exit='C'), ]
                + [Segment([helix], entry='N', exit='C')] * 13
                + [Segment([helix], entry='N')])
    worms = grow(segments,
                 SegmentXform('C1', lever=20),
                 memlim=float(sys.argv[-1][1:]))
    print(worms.scores)


if __name__ == '__main__':
    main()
