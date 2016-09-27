#include "gtest_util.hpp"

// include desired tests here:

#include "test/hana/hana.gtest.cc"

// #include "riflib/util/hash.gtest.cc"
// #include "riflib/util/hash_thread.gtest.cc" // TODO ARRAY!!
// #include "riflib/util/bloom.gtest.cc"
// #include "riflib/util/StoragePolicy.gtest.cc"
// #include "riflib/util/template_loop.gtest.cc" // 3.5 vs 2.2 wo Eigen
// #include "riflib/util/SimpleArray.gtest.cc"

// #include "riflib/numeric/euler_angles.gtest.cc"
// #include "riflib/numeric/FixedPoint.gtest.cc"
// #include "riflib/numeric/xform_perf.gtest.cc"
// #include "riflib/numeric/rand_xform.gtest.cc"
// #include "riflib/dock/fftdock.gtest.cc"

// #include "riflib/nest/NEST.gtest.cc" // 4.479 vs. 3.0 wo Eigen
// #include "riflib/nest/NEST_neighbor.gtest.cc" // 4.567 // 3.26 wo Eigen
// #include "riflib/nest/MultiNest.gtest.cc"

// #include "riflib/nest/pmap/parameter_maps.gtest.cc" // 7.546 // 4.0 wo Eigen
// #include "riflib/nest/pmap/parameter_maps_test_nbrcell.gtest.cc" // 4.472 // 3.1 wo eigen
// #include "riflib/nest/pmap/SphereDodec.gtest.cc"
// #include "riflib/nest/pmap/SphereQuad.gtest.cc"
// #include "riflib/nest/pmap/HecatonicosachoronMap.gtest.cc"
// #include "riflib/nest/pmap/QuaternionMap.gtest.cc"
// #include "riflib/nest/pmap/EulerAnglesMap.gtest.cc"
// #include "riflib/nest/pmap/TetracontoctachoronMap.gtest.cc"
// #include "riflib/nest/pmap/Rotation1DMap.gtest.cc"
// #include "riflib/nest/pmap/OriTransMap.gtest.cc"

// #include "riflib/numeric/geom_4d.gtest.cc"
// #include "riflib/numeric/bcc_lattice.gtest.cc"
// #include "riflib/numeric/bcc_lattice_orientation.gtest.cc"

// #include "riflib/util/meta/util.gtest.cc"
// #include "riflib/util/container/ContainerInteractions.gtest.cc"
// #include "riflib/util/meta/InstanceMap.gtest.cc"
// #include "riflib/util/meta/InstanceMap_container.gtest.cc"
// #include "riflib/util/meta/InstanceMap_numeric.gtest.cc"

// #include "riflib/objective/ObjectiveFunction.gtest.cc"
// #include "riflib/objective/voxel/VoxelArray.gtest.cc"
// #include "riflib/objective/voxel/FieldCache.gtest.cc"

// #include "riflib/kinematics/Scene.gtest.cc"
// #include "riflib/kinematics/Scene_test_eigen.gtest.cc"
// #include "riflib/kinematics/SceneIterator.gtest.cc"
// #include "riflib/kinematics/Director.gtest.cc"


// #include "riflib/io/dump_pdb_atom.gtest.cc"

// #include "riflib/kinematics/Scene_test_objective.gtest.cc"

// #include "riflib/actor/BBStub.gtest.cc"

// #include "riflib/rosetta/score/AnalyticEvaluation.gtest.cc"
// #include "riflib/rosetta/score/RosettaField.gtest.cc"

// #include "riflib/chemical/ligand_factory.gtest.cc"

// #include "riflib/objective/methods/hbond_5dof.gtest.cc"

// #include "riflib/objective/hash/XformHash.gtest.cc"
// #include "riflib/objective/hash/XformHashFromNest.gtest.cc"
// #include "riflib/objective/hash/XformHashNeighbors.gtest.cc"
// #include "riflib/objective/hash/XformMap.gtest.cc"

// #include "riflib/objective/storage/RotamerScores.gtest.cc"

// #include "riflib/chemical/RotamerIndex.gtest.cc"

// #include "riflib/objective/storage/TwoBodyTable.gtest.cc"

// #include "riflib/actor/BackboneActor.gtest.cc"

// #include "riflib/search/SpatialBandB.gtest.cc"

// #include "riflib/chemical/stub.gtest.cc"

// #include "riflib/numeric/cube_to_sphere.gtest.cc"

// #include "riflib/numeric/rand_xform.gtest.cc"

int main(int argc, char **argv)
{

	std::cout << int( 1.6 * 1 ) << std::endl;

	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
