"""
a prototype module for fiber docking. for illustration see:

https://docs.google.com/presentation/d/1qv6yZIy0b2dfmZAzA7bFNJD4pl9h7Nt1mOL4aqiVpOw/pub?start=false&loop=false&delayms=3000#slide=id.g1d6846b840_0_8

(see prev slide for simpler, less-general case)
"""

from rif.io.pdbio import read_pdb
from rif.nest import NestedXforms
from rif.models import ImplicitEnsemble, ImplicitEnsembleScore
from rif.search import TieredBeamSearch

import numpy as np


def rel_prime_flat_tri(Nmax):
    pass


def center_welzl():
    pass


def fiber_dock(pdb_file):
    """turn monomer into helical fiber in "all" possible ways"""

    # setup
    resolutions = (16, 8, 4, 2, 1, 0.5)
    M, N = rel_prime_flat_tri(Nmax=20)  # 1,2 1,3 2,3 1,4 3,4 1,5 2,5 3,5 4,5 ...
    powers = (M / N)[np.newaxis, :]  # shape (1, whatever)
    assert all(powers < 1)  # else Hsearch can't work
    prot, radius = center_welzl(read_pdb(pdb_file))
    bbframes = prot.bbframes(ss='HE', burial='exposed')
    samp_grid = NestedXforms(baseresl=resolutions[0], lever=radius)
    # ImplicitEnsemble represents a structural ensemble in which each member's
    # elements (atoms, rays, stubs, etc) are within well-defined translational and
    # rotational bounds of a reference structure.
    # (rotation mag defined by radius/lever)
    ensemble0 = ImplicitEnsemble(atoms=prot.atoms, rotslots=bbframes,
                                 radius=0.0, lever=radius, bvh_max_radius=1.0,
                                 voxel_convolv_radii=resolutions)
    hsearch = TieredBeamSearch(resls=resolutions, grid=samp_grid, beam_size_m=10,
                               extra_fields=(('M', 'u1'), ('N', 'u1')), )
    score = ImplicitEnsembleScore(
        # todo
    )

    # assuming I can do what I want with numpy dtype/operators...
    def score_samples(stage, samples):
        resl = resolutions[stage]
        indices = samples['index']
        positions = samp_grid.position_of(stage, indices)
        ensembleN = positions * ensemble0.with_radius(resl)
        score0N = score(ensemble0, ensembleN)  # radius 0 + resl
        # better to extract sub-array of score0N < 0?
        if np.any(score0N < 0):
            # in addition to ensemble['position']**(M/N), operator
            # __pow__(ensemble,float) must scale ensemble['radius']
            # linearly by pow... is this well-defined algebra?
            # maybe just use an ensemble_power() func directly?
            ensembleM = np.power(ensembleN, powers)
            # radius 0 + resl * m / n (not constant over array)
            score0M = score(ensemble0, ensembleM)
            # - coupling='linear' means M ensemble motion tied exactly to
            #   N motion, so instead of radius = Arad + Brad,
            #   radius = min(Arad - Brad, Brad - Arad)
            # - this check may not be worth the cost here, skip or
            #   maybe filter on ensembleM good enough first
            scoreMN = score(ensembleM, ensembleN, coupling='linear')
            score0MN = score0N + score0M + scoreMN  # score0N bcast
            idxmin = score0MN.argmin(axis=2)  # best score forall M / N
            samples['M'] = M[idxmin]
            samples['N'] = N[idxmin]
            samples['score'] = score0MN[:, idxmin]
            # todo: eventually check all good M/N? will probably be rare
            #       to have more than one hit
            # todo: also score ensemble0 vs all 1..N) for thoroughness
        return samples

    # hsearch interface option 1
    # simple, but perhaps not as flexible as below
    hsearch.search(score_samples)
    print(hsearch.results())

    # hsearch interface option 2
    # not sure if using iteration as an interface to the search is a
    # good idea or not... note that samples['score'] gets communicated
    # back to the search, which is not a pattern I've seen used much
    # but it is clearer in the sense that you don't actually have to
    # define a function (ie score_samples)
    hsearch_accum = hsearch.accumulator()
    for iresl, samples in hsearch_accum:
        score_samples(iresl, samples)
    print(hsearch_accum.results())

    # hsearch interface option 3
    # more of an async / await co-routine?


if __name__ == '__main__':
    fiber_dock()
