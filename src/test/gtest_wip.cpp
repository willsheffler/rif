#include "gtest_util.hpp"

#define SCHEME_BENCHMARK

// include desired tests here:

// #include "numeric/eigen_xform_perf.gtest.cpp"

// #include "test/hana/hana.gtest.cpp"

// #include "util/hash.gtest.cpp"
// #include "util/hash_thread.gtest.cpp" // TODO ARRAY!!
// #include "util/bloom.gtest.cpp"
// #include "util/StoragePolicy.gtest.cpp"
// #include "util/template_loop.gtest.cpp" // 3.5 vs 2.2 wo Eigen
// #include "util/SimpleArray.gtest.cpp"

// #include "numeric/euler_angles.gtest.cpp"
// #include "numeric/FixedPoint.gtest.cpp"
// #include "geom/rand_geom.gtest.cpp"
// #include "dock/fftdock.gtest.cpp"

// #include "nest/NEST.gtest.cpp" // 4.479 vs. 3.0 wo Eigen
// #include "nest/NEST_neighbor.gtest.cpp" // 4.567 // 3.26 wo Eigen
// #include "nest/MultiNest.gtest.cpp"

// #include "nest/pmap/parameter_maps.gtest.cpp" // 7.546 // 4.0 wo Eigen
// #include "nest/pmap/parameter_maps_test_nbrcell.gtest.cpp" // 4.472 // 3.1 wo
// eigen
// #include "nest/pmap/SphereDodec.gtest.cpp"
// #include "nest/pmap/SphereQuad.gtest.cpp"
// #include "nest/pmap/HecatonicosachoronMap.gtest.cpp"
// #include "nest/pmap/QuaternionMap.gtest.cpp"
// #include "nest/pmap/EulerAnglesMap.gtest.cpp"
// #include "nest/pmap/TetracontoctachoronMap.gtest.cpp"
// #include "nest/pmap/Rotation1DMap.gtest.cpp"
// #include "nest/pmap/OriTransMap.gtest.cpp"

// #include "numeric/geom_4d.gtest.cpp"
// #include "numeric/lattice.gtest.cpp"
// #include "numeric/lattice_orientation.gtest.cpp"

// #include "util/meta/util.gtest.cpp"
// #include "util/container/ContainerInteractions.gtest.cpp"
// #include "util/meta/InstanceMap.gtest.cpp"
// #include "util/meta/InstanceMap_container.gtest.cpp"
// #include "util/meta/InstanceMap_numeric.gtest.cpp"

// #include "objective/ObjectiveFunction.gtest.cpp"
// #include "objective/voxel/VoxelArray.gtest.cpp"
// #include "objective/voxel/FieldCache.gtest.cpp"

// #include "kinematics/Scene.gtest.cpp"
// #include "kinematics/Scene_test_eigen.gtest.cpp"
// #include "kinematics/SceneIterator.gtest.cpp"
// #include "kinematics/Director.gtest.cpp"

// #include "io/dump_pdb_atom.gtest.cpp"

// #include "kinematics/Scene_test_objective.gtest.cpp"

// #include "actor/BBStub.gtest.cpp"

// #include "rosetta/score/AnalyticEvaluation.gtest.cpp"
// #include "rosetta/score/RosettaField.gtest.cpp"

// #include "chem/ligand_factory.gtest.cpp"

// #include "objective/methods/hbond_5dof.gtest.cpp"

// #include "objective/hash/XformHash.gtest.cpp"
// #include "objective/hash/XformHashFromNest.gtest.cpp"
// #include "objective/hash/XformHashNeighbors.gtest.cpp"
// #include "objective/hash/XformMap.gtest.cpp"

// #include "objective/storage/RotamerScores.gtest.cpp"

// #include "chem/RotamerIndex.gtest.cpp"

// #include "objective/storage/TwoBodyTable.gtest.cpp"

// #include "actor/BackboneActor.gtest.cpp"

// #include "search/SpatialBandB.gtest.cpp"

// #include "chem/stub.gtest.cpp"

#include "geom/Ray.gtest.cpp"
#include "geom/cube_to_sphere.gtest.cpp"

// #include "sampling/orientations.gtest.cpp"

// #include "geom/rand_geom.gtest.cpp"

int main(int argc, char **argv) {
  std::cout << int(1.6 * 1) << std::endl;

  std::vector<std::string> args;
  for (int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
  std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
  init_gtest_tests(args);
  return run_gtest_tests();
}
