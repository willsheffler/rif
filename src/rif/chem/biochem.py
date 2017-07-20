from builtins import *
from rif import rcl


aa_name1s = list('ACDEFGHIKLMNPQRSTVWY')

aa_name3s = ('ala cys asp glu phe gly his ile lys leu ' +
             'met asn pro gln arg ser thr val trp tyr').upper().split()

rif_atype_names = (
    """
    Nbb  CAbb CObb OCbb CH2
    CH3  HNbb Hapo Hpol Haro
    CNH2 COO  aroC Ntrp Nhis
    NH2O Nlys Narg Npro OH
    ONH2 OOC  S    HOH
    """
).split()
assert len(rif_atype_names) == len(set(rif_atype_names))

mm_atype_name = {}


rif_atype_index = dict()
for i, v in enumerate(rif_atype_names):
    rif_atype_index[v] = i


def check_atom_types(ary):
    ary = ary[ary['atype'] == 255]
    for atom in ary:
        rname3 = aa_name3s[atom['rtype']]
        aname = 'UNKNOWN'
        if rcl.HAVE_PYROSETTA:
            aname = rcl.rts().name_map(rname3).atom_name(atom['anum'] + 1)
        print('unk atom type', rname3, aname, atom)
    assert len(ary) == 0


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
