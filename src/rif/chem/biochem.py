
from builtins import *

aa_name1s = list('ACDEFGHIKLMNPQRSTVWY')
aa_name3s = ('ala cys asp glu phe gly his ile lys leu ' +
             'met asn pro gln arg ser thr val trp tyr').upper().split()
rif_atype_names = ('CNH2 COO CH1 CH2 CH3 aroC Ntrp Nhis NH2O Nlys ' +
                   'Narg Npro OH ONH2 OOC S Nbb CAbb CObb OCbb Hpol ' +
                   'Hapo Haro HNbb VIRT CH0 HS NtrR SH1').split()

mm_atype_name = {}


rif_atype_index = dict()
for i, v in enumerate(rif_atype_names):
    rif_atype_index[v] = i


def compute_chi_levers():
    import rif.rcl as rcl
    rcl.pyrosetta.init('-beta -include_patches patches/VirtualBB.txt')
    chem_manager = rcl.rosetta.core.chemical.ChemicalManager
    rts = chem_manager.get_instance().residue_type_set("fa_standard")
    for aa in aa_name3s:
        print(aa)
        rots = rcl.generate_canonical_rotamer_residues(aa)
        for rot in rots:
            for ichi in range(rot.nchi()):
                print(aa, ichi, rot)
