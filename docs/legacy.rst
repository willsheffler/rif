--------
Legacy
--------


This section documents the original RIF docking method and the legacy C++ implementation. Docs are here because (1) this is a convenient place for me to put them, and (2) the current library has grown out of the C++ implementation.

Prototype one-sided design implementation
===========================================

Overview
-------------

Rotamer Interaction Field (RIF) Docking for de novo design of binders to small molecules and fixed protein targets works in two phases: RIF generation and scaffold docking.

RIF Generation
~~~~~~~~~~~~~~~~~~~
First, a RIF
tailored to the small molecule target is generated. The RIF contains billions of
rotamers that interact productively with the target in both polar and apolar
fashion. Polar interactions with the ligand are generated based on hydrogen bond
geometry and are tagged with the polar group(s) on the small molecule with which
they hydrogen bond (Figure XXX, groups with numeric labels). All possible polar
interactions are included in the RIF. Apolar interacting rotamers are generated
with a docking process that generates all possible packing interactions above an
configurable amino acid specific rosetta energy threshold (Figure XXX label
APO). The polar and apolar RIF rotamers are stored ~0.5Å gridded representation
(REF 1) of the six dimensional rigid body space, based on the backbone position
of the rotamers generated. To facilitate search, described next, lower
resolution RIF data structures are produced from the primary, gridded at 1.0Å,
2.0Å, 4.0Å, 8.0Å and 16.0Å resolutions.`

Implemented in binary rifgen_

RIF Docking
~~~~~~~~~~~~~~~
Second, after RIF generation, a set of protein backbone scaffolds selected by
the user is docked into the RIF
using a hierarchical branch and bound style search (REF 2). Starting with the
coarsest 16Å resolution RIF, we perform an enumerative search of possible 6D
positions of each scaffold. Foreach scaffold placement, the designable scaffold
backbone positions, as specified by the user, are checked against the RIF to
determine a set of productively interacting rotamers which can be placed on the
scaffold. If the rotamers found satisfy the all polar groups on the small
molecule target (or a user-specified acceptable fraction), the scaffold position
is assigned a score based on the quality of the polar and apolar interactions
that can be formed, else it is rejected. All acceptable scaffold positions up to
a configurable limit (typically 10 million) are sorted by score and promoted to
the next search stage. Each promoted scaffold is split into 64 child positions
by moving the scaffold in the 6D rigid body space, providing a finer sampling
appropriate for the next RIF resolution. Search is done at 16Å, 8Å, 4Å, 2Å, 1Å,
and 0.5Å resolutions. After the final 0.5Å search stage, a monte-carlo based
combinatorial rotamer packing step is performed on a configurable fraction of
the best scaffold placements (typically 1 million) to find internally consistent
rotamer placements. The docking process is highly optimized and takes roughly
one cpu hour per scaffold to generate thousands of possible binding modes.

Implemented in binary rif_dock_test_

Building
---------
The current RIF docking protocol is not implemented in python and is not part of this repository.

Building Rosetta
~~~~~~~~~~~~~

There is only one dependency for the legacy rif implementation, but it's a whopper: Rosetta_. You will need a relatively current clone of rosetta main branch master from http://github.com/rosettacommons/main. Building with a rosetta release version should be possible, but has not been tested.

You must use rosetta's cmake build system as follows:

.. code-block:: bash

    cd <path_to_rosetta>/source/cmake
    python make_project.py all
    cd build_cxx11_omp
    cmake .
    make -j <num_processors>

The rosetta build will take a few cpu hours.

.. _Rosetta: http://www.rosettacommons.org

Obtain legacy scheme source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, obtain the legacy rif source you want to build. There are a few choices. The 'standard' version most people in the bakerlab are using is from willsheffler@gmail.com. It has been largely unchanged for a year and should be relatively stable. There are also versions which have been enhanced with some additional features by other members of the baker lab, most notably one available from Brian Coventry https://github.com/bcov77/scheme/tree/brian.

Building RIF
~~~~~~~~~~~~~

The following should produce a build:

.. code-block:: bash

    git clone git@github.com:willsheffler/scheme
    cd scheme
    git checkout sheffler/devel

The above will get you the 'stable' version from willsheffler@gmail.com. If you with to use another, substitute the repository and branch as appropriate. For example, for Brian Coventry's version:

.. code-block:: bash

    git clone git@github.com:bcov77/scheme
    cd scheme
    git checkout brian

Do the following to produce build config with cmake and build (ninja can be used instead of make). Note that the CMAKE_ROSETTA_PATH environment variable must be specified as this code links against rosetta.

.. code-block:: bash

    export CMAKE_ROSETTA_PATH='<path_to_rosetta>'
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j <num_processors> rifgen rif_dock_test



Usage
------

The build will product two binaries: rifgen_ for `Rif generation`_, and rif_dock_test_ for `RIF Docking`_. Both binaries are multi-threaded via OpenMP. You may control the number of threads used with the OMP_NUM_THREADS environment variable. The default is to use all available hardware threads (which is probably best). Running rifgen and rif_dock_test will in most require at least 64gb of ram per process, hence, it is best to use more threads per process instead of multiple processes. For some cases, particularly if you are not generating a full de novo RIF but are instead manually inputting hotspots, you may be able to get away with less ram.

You really *really* want to contact somebody who is using rif for a similar problem as you are and use their flags files as a starting point. Broadly speaking, people in the Baker lab are using rif for one-sided protein interface design (particularly denovo) and for small molecule binder design. If you contact willsheffler@gmail.com I can put you in touch with somebody.

rifgen
~~~~~~~~~~
like most rosetta-based programs, usage is: rifgen @flags_file

rifgen makes use of a fair amount of NON-TARGET-SPECIFIC cached data. In most cases, this data will be generated if it is not available. This means the first time you run rifgen, it could take a *long* time. This data is stored and accessed in the location passed to '-rifgen:data_cache_dir'. As this data is not target specific, it can be shared between runs, projects, and users. Most people in the baker lab point to one such directory, for example.

**flags you will likely need to modify:**

::

    -rifgen:data_cache_dir <location_for_cached_data>

    -database <path to rosetta database, MUST be the same database as the rosetta you built>

    -rifgen::rif_type RotScore  # use RotScoreSat for small molecules

    -rifgen:target     <pdb file>
    -rifgen:target_res <residue_numbers.txt> OPTIONAL
    -rifgen:outdir     <output_dir_location>
    -rifgen:outfile    <output_file_name.rif.gz>

    -rifgen:beam_size_M 10.0

**flags you may want need to tweak:**

::

    # apolar residues to generate in de novo RIF
    -rifgen:apores VAL ILE LEU MET PHE TRP

    # polar residues to generate in de novo RIF
    -rifgen:donres SER THR TYR     GLN     ASN HIS HIS_D TRP
    -rifgen:accres SER THR TYR GLU GLN ASP ASN HIS HIS_D

    -rifgen:score_cut_adjust 1.0   # scale factor for hydrophobic score cutoffs

    -rifgen::rosetta_field_resl 0.25          # grid spacing for vdw/lksol scoring
    -rifgen::search_resolutions 3.0 1.5 0.75  # search resl for denovo hydrophobic interactions

    -hash_cart_resl    0.7  # memory use will scale as (1/resl)**6 so be careful!
    -hash_angle_resl  14.0  # should stay roughly 20x hash_cart_resl

    -rif_accum_scratch_size_M 24000  # megabytes of memory to use as a scratch space

    -rifgen:score_threshold 0          # score cut for rotamer inclusion in rif
    -rifgen:hbond_weight 1.0           # max score per-hbond
    -rifgen:upweight_multi_hbond 1.0   # extra score factor for bidentate hbonds

    # for sanity-checking only, you can dump a (small!) fraction RIF rotamers:
    -rifgen:rif_hbond_dump_fraction  0.000001
    -rifgen:rif_apo_dump_fraction    0.000001


**flags you must have but probably shouldn't change:**

::

    -rifgen:extra_rotamers false
    -rifgen:extra_rif_rotamers true

    # rosetta stuff
    -renumber_pdb
    -add_orbitals

    -hbond_cart_sample_hack_range 0.33
    -hbond_cart_sample_hack_resl  0.33
    -rifgen:tip_tol_deg        60.0 # for now, do either 60 or 36
    -rifgen:rot_samp_resl       6.0




    -rifgen:hash_preallocate_mult 0.125
    -rifgen:max_rf_bounding_ratio 4.0

    # geometry of bounding grids
    -rifgen:hash_cart_resls   16.0   8.0   4.0   2.0   1.0
    -rifgen:hash_cart_bounds   512   512   512   512   512
    -rifgen:lever_bounds      16.0   8.0   4.0   2.0   1.0
    -rifgen:hash_ang_resls     38.8  24.4  17.2  13.6  11.8 # yes worky worky
    -rifgen:lever_radii        23.6 18.785501 13.324600  8.425850  4.855575


**Output**

All products of rifgen will be in the directory specified by the '-output_dir' flag. In addition, rifgen will dump output like the following at the end of the run. You can copy and paste this output into your rif_dock_test flags file to avoid some tedious and error-prone editing.

::

    ########################################### what you need for docking ###########################################
    -rif_dock:target_pdb            rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_target.pdb.gz
    -rif_dock:target_res            test_input/fz7/Fz7_model.res
    -rif_dock:target_rf_resl        0.125
    -rif_dock:target_rf_cache       rifgen_fz7_v4/__RF_Fz7_model.pdb_CEN_trhash12511471_resl0.125_osamp2_replonlybdry
    -rif_dock:target_bounding_xmaps rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_BOUNDING_RIF_16.xmap.gz
    -rif_dock:target_bounding_xmaps rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_BOUNDING_RIF_08.xmap.gz
    -rif_dock:target_bounding_xmaps rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_BOUNDING_RIF_04.xmap.gz
    -rif_dock:target_bounding_xmaps rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_BOUNDING_RIF_02.xmap.gz
    -rif_dock:target_bounding_xmaps rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif_BOUNDING_RIF_01.xmap.gz
    -rif_dock:target_rif            rifgen_fz7_v4/rif_64_fz7_v4_sca0.7_noKR.rif.gz
    -rif_dock:extra_rotamers        0
    -rif_dock:extra_rif_rotamers    1
    #################################################################################################################




rif_dock_test
~~~~~~~~~~~~~~


Slide Decks
===================



IPD Presentation 2016/08
-------------------------

    more "polish" than 2016/04, benchmark data here are accurate, to best of my knowledge

.. raw:: html

    <iframe src="https://docs.google.com/presentation/d/1We5liWBFqhYPqFNIoRJbSjzR3LFGCFslLg6hshDZNJE/embed?start=false&loop=false&delayms=3000" frameborder="0" width="1365" height="1053" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe><br><br><br>



IPD Presentation 2016/04
------------------------

    more technical info vs 2016/08, protein docking DDG benchmark numbers **inaccurate** here!

.. raw:: html

    <iframe src="https://docs.google.com/presentation/d/1zEMJ3KQwS8i7S8bdOSKbTcCzBZKrR2oz8fj8FFDUyQE/embed?start=false&loop=false&delayms=3000" frameborder="0" width="1365" height="1053" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe><br><br><br>


