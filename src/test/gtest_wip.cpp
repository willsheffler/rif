#include "gtest_util.hpp"

#define SCHEME_BENCHMARK

// include desired tests here:

// #include "riflib/numeric/eigen_xform_perf.gtest.cpp"


// #include "test/hana/hana.gtest.cpp"

// #include "riflib/util/hash.gtest.cpp"
// #include "riflib/util/hash_thread.gtest.cpp" // TODO ARRAY!!
// #include "riflib/util/bloom.gtest.cpp"
// #include "riflib/util/StoragePolicy.gtest.cpp"
// #include "riflib/util/template_loop.gtest.cpp" // 3.5 vs 2.2 wo Eigen
// #include "riflib/util/SimpleArray.gtest.cpp"

// #include "riflib/numeric/euler_angles.gtest.cpp"
// #include "riflib/numeric/FixedPoint.gtest.cpp"
// #include "riflib/numeric/rand_xform.gtest.cpp"
// #include "riflib/dock/fftdock.gtest.cpp"

// #include "riflib/nest/NEST.gtest.cpp" // 4.479 vs. 3.0 wo Eigen
// #include "riflib/nest/NEST_neighbor.gtest.cpp" // 4.567 // 3.26 wo Eigen
// #include "riflib/nest/MultiNest.gtest.cpp"

// #include "riflib/nest/pmap/parameter_maps.gtest.cpp" // 7.546 // 4.0 wo Eigen
// #include "riflib/nest/pmap/parameter_maps_test_nbrcell.gtest.cpp" // 4.472 // 3.1 wo eigen
// #include "riflib/nest/pmap/SphereDodec.gtest.cpp"
// #include "riflib/nest/pmap/SphereQuad.gtest.cpp"
// #include "riflib/nest/pmap/HecatonicosachoronMap.gtest.cpp"
// #include "riflib/nest/pmap/QuaternionMap.gtest.cpp"
// #include "riflib/nest/pmap/EulerAnglesMap.gtest.cpp"
// #include "riflib/nest/pmap/TetracontoctachoronMap.gtest.cpp"
// #include "riflib/nest/pmap/Rotation1DMap.gtest.cpp"
// #include "riflib/nest/pmap/OriTransMap.gtest.cpp"

// #include "riflib/numeric/geom_4d.gtest.cpp"
// #include "riflib/numeric/bcc_lattice.gtest.cpp"
// #include "riflib/numeric/bcc_lattice_orientation.gtest.cpp"

// #include "riflib/util/meta/util.gtest.cpp"
// #include "riflib/util/container/ContainerInteractions.gtest.cpp"
// #include "riflib/util/meta/InstanceMap.gtest.cpp"
// #include "riflib/util/meta/InstanceMap_container.gtest.cpp"
// #include "riflib/util/meta/InstanceMap_numeric.gtest.cpp"

// #include "riflib/objective/ObjectiveFunction.gtest.cpp"
// #include "riflib/objective/voxel/VoxelArray.gtest.cpp"
// #include "riflib/objective/voxel/FieldCache.gtest.cpp"

// #include "riflib/kinematics/Scene.gtest.cpp"
// #include "riflib/kinematics/Scene_test_eigen.gtest.cpp"
// #include "riflib/kinematics/SceneIterator.gtest.cpp"
// #include "riflib/kinematics/Director.gtest.cpp"


// #include "riflib/io/dump_pdb_atom.gtest.cpp"

// #include "riflib/kinematics/Scene_test_objective.gtest.cpp"

// #include "riflib/actor/BBStub.gtest.cpp"

// #include "riflib/rosetta/score/AnalyticEvaluation.gtest.cpp"
// #include "riflib/rosetta/score/RosettaField.gtest.cpp"

// #include "riflib/chemical/ligand_factory.gtest.cpp"

// #include "riflib/objective/methods/hbond_5dof.gtest.cpp"

// #include "riflib/objective/hash/XformHash.gtest.cpp"
// #include "riflib/objective/hash/XformHashFromNest.gtest.cpp"
// #include "riflib/objective/hash/XformHashNeighbors.gtest.cpp"
// #include "riflib/objective/hash/XformMap.gtest.cpp"

// #include "riflib/objective/storage/RotamerScores.gtest.cpp"

// #include "riflib/chemical/RotamerIndex.gtest.cpp"

// #include "riflib/objective/storage/TwoBodyTable.gtest.cpp"

// #include "riflib/actor/BackboneActor.gtest.cpp"

// #include "riflib/search/SpatialBandB.gtest.cpp"

// #include "riflib/chemical/stub.gtest.cpp"

#include "riflib/numeric/cube_to_sphere.gtest.cpp"
#include "sampling/Orientations.gtest.cpp"

// #include "riflib/numeric/rand_xform.gtest.cpp"

int main(int argc, char **argv)
{

	std::cout << int( 1.6 * 1 ) << std::endl;

	std::vector<std::string> args;
	for(int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
	std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
	init_gtest_tests(args);
	return run_gtest_tests();
}
