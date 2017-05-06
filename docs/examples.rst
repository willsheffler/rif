========
Examples
========

---------------------------
high level
---------------------------

fiber docking
--------------

    This is prototype code to dock a monomer into helical protein fiber. Xform from :math:`B_0` to :math:`B_N` is sampled via hierarchical search. Possible placements a third body :math:`B_M` are checked based on possible :math:`M/N` roots of the transform :math:`X_N`:

    .. math::
        X_M = X_N^{M/N}

    .. image:: img/fiber_dock_B0_BN_Bm.png
       :width: 50%

    .. literalinclude:: ../src/rif/apps/fiber_dock.py

