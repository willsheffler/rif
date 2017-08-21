from rif import rcl
if rcl.HAVE_PYROSETTA:
    from rif.rcl import atoms, stubs, rays
    from rif.index import stripe_index_3d


def main():
    # shift NH index down by one
    pept0 = rcl.atoms(peptpose, 'N ca c o')  # 344 ires,iatm,crd

    # 2224 ires, conh, orgdirn,crd
    rays0 = rcl.rays(pept, 'c->o n->h', shift=(0, -1))

    raycen0 = rays0[..., 'orig'].sum(axis='conh') / 2.0  # 24 ires,crd

    bbcoords = rcl. atoms(scafpose, 'N CA c o cb (ca)')
    bbstubs = rcl.stubs(bbcoords, 'N ca c')
    clasher = stripe_index_3d(bbcoords, 3.5)
    contacter = stripe_index_3d(bbcoords, stubs, 5.0)

    # haxis = helix axis in pept cpept = atoms(peptpose, 'N ca c o')
    rays = rays(pept, 'c->o n->h')

    samppos0 = xforms_around_axis(
        axis=haxis, cen=hcen, name='sampcyl',
        radius=(4, 8), resl=0.5, bounds=(-10, 10),

    )
    pept = samppos0 * pept0  # outer n344

    sampclash = clasher.clashes(pept).any(axis=('ires', 'iatm'))
    samppos = samppos[~sampclash]
    pept = samppos * pept0  # m344

    raycen = sampos * raycen0  # m24 sampcyl, ires, crd
    jagged_stubs = contacter.contacts(raycen)  # jagged m*2*x,44

    rays = samppos * rays0  # m2224 sampcyl, ires, conh, origdirn, crd
    stublocalrays = jagged_stubs * rays  # m*2*x,224
    raykeys = ray_hasher10d.get_index(stublocalrays)  # m*2*x
    rayraykeys.dictmap(whateverhash)  # ? How to merge?


class Jagged:
    # DataArray source
    # n44 array data
    # nx array index, x is prefix of source dims?
    def __mul__(self, dataarray_da):
        # Check da marches req. Dims
        # Remap dims to linear on index
        return self * daflat[index]


if __name__ == '__main__':
    main()
