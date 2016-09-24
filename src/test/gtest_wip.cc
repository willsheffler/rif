#include "gtest_util.hh"

// include desired tests here:

#include "test/hana/hana.gtest.cc"

// #include "scheme/util/hash.gtest.cc"
// #include "scheme/util/hash_thread.gtest.cc" // TODO ARRAY!!
// #include "scheme/util/bloom.gtest.cc"
// #include "scheme/util/StoragePolicy.gtest.cc"
// #include "scheme/util/template_loop.gtest.cc" // 3.5 vs 2.2 wo Eigen
// #include "scheme/util/SimpleArray.gtest.cc"

// #include "scheme/numeric/euler_angles.gtest.cc"
// #include "scheme/numeric/FixedPoint.gtest.cc"
// #include "scheme/numeric/xform_perf.gtest.cc"
// #include "scheme/numeric/rand_xform.gtest.cc"
// #include "scheme/dock/fftdock.gtest.cc"

// #include "scheme/nest/NEST.gtest.cc" // 4.479 vs. 3.0 wo Eigen
// #include "scheme/nest/NEST_neighbor.gtest.cc" // 4.567 // 3.26 wo Eigen
// #include "scheme/nest/MultiNest.gtest.cc"

// #include "scheme/nest/pmap/parameter_maps.gtest.cc" // 7.546 // 4.0 wo Eigen
// #include "scheme/nest/pmap/parameter_maps_test_nbrcell.gtest.cc" // 4.472 // 3.1 wo eigen
// #include "scheme/nest/pmap/SphereDodec.gtest.cc"
// #include "scheme/nest/pmap/SphereQuad.gtest.cc"
// #include "scheme/nest/pmap/HecatonicosachoronMap.gtest.cc"
// #include "scheme/nest/pmap/QuaternionMap.gtest.cc"
// #include "scheme/nest/pmap/EulerAnglesMap.gtest.cc"
// #include "scheme/nest/pmap/TetracontoctachoronMap.gtest.cc"
// #include "scheme/nest/pmap/Rotation1DMap.gtest.cc"
// #include "scheme/nest/pmap/OriTransMap.gtest.cc"

// #include "scheme/numeric/geom_4d.gtest.cc"
// #include "scheme/numeric/bcc_lattice.gtest.cc"
// #include "scheme/numeric/bcc_lattice_orientation.gtest.cc"

// #include "scheme/util/meta/util.gtest.cc"
// #include "scheme/util/container/ContainerInteractions.gtest.cc"
// #include "scheme/util/meta/InstanceMap.gtest.cc"
// #include "scheme/util/meta/InstanceMap_container.gtest.cc"
// #include "scheme/util/meta/InstanceMap_numeric.gtest.cc"

// #include "scheme/objective/ObjectiveFunction.gtest.cc"
// #include "scheme/objective/voxel/VoxelArray.gtest.cc"
// #include "scheme/objective/voxel/FieldCache.gtest.cc"

// #include "scheme/kinematics/Scene.gtest.cc"
// #include "scheme/kinematics/Scene_test_eigen.gtest.cc"
// #include "scheme/kinematics/SceneIterator.gtest.cc"
// #include "scheme/kinematics/Director.gtest.cc"


// #include "scheme/io/dump_pdb_atom.gtest.cc"

// #include "scheme/kinematics/Scene_test_objective.gtest.cc"

// #include "scheme/actor/BBStub.gtest.cc"

// #include "scheme/rosetta/score/AnalyticEvaluation.gtest.cc"
// #include "scheme/rosetta/score/RosettaField.gtest.cc"

// #include "scheme/chemical/ligand_factory.gtest.cc"

// #include "scheme/objective/methods/hbond_5dof.gtest.cc"

// #include "scheme/objective/hash/XformHash.gtest.cc"
// #include "scheme/objective/hash/XformHashFromNest.gtest.cc"
// #include "scheme/objective/hash/XformHashNeighbors.gtest.cc"
// #include "scheme/objective/hash/XformMap.gtest.cc"

// #include "scheme/objective/storage/RotamerScores.gtest.cc"

// #include "scheme/chemical/RotamerIndex.gtest.cc"

// #include "scheme/objective/storage/TwoBodyTable.gtest.cc"

// #include "scheme/actor/BackboneActor.gtest.cc"

// #include "scheme/search/SpatialBandB.gtest.cc"

// #include "scheme/chemical/stub.gtest.cc"

// #include "scheme/numeric/cube_to_sphere.gtest.cc"

// #include "scheme/numeric/rand_xform.gtest.cc"

int main(int argc, char **argv)
{

	std::cout << int( 1.6 * 1 ) << std::endl;

	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
