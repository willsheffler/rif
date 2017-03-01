
from builtins import *


aa_name3s = 'ala cys asp glu phe gly his ile lys leu met asn pro gln arg ser thr val trp tyr'.upper().split()


rif_atype_names = ('CNH2 COO CH1 CH2 CH3 aroC Ntrp Nhis NH2O Nlys Narg Npro OH ONH2 OOC S'.split() +
                   'Nbb CAbb CObb OCbb Hpol Hapo Haro HNbb VIRT CH0 HS NtrR SH1'.split())

mm_atype_name = {}


rif_atype_index = dict()
for i, v in enumerate(rif_atype_names):
    rif_atype_index[v] = i
