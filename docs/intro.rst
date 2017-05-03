-------------------
RIF Introduction
-------------------

Prototype one-sided design implementation
===========================================

Overview
-------------

Rotamer Interaction Field (RIF) Docking for small molecule binder de novo design
works in two phases: RIF generation and scaffold docking.

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
2.0Å, 4.0Å, 8.0Å and 16.0Å resolutions.

RIF Docking
~~~~~~~~~~~~~~~
Second, after RIF generation, a set of protein backbone scaffolds selected by the user is docked into the RIF
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
