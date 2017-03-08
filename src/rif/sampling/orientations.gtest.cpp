#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <iostream>
#include <sampling/orientations.hpp>

#include <sys/param.h>
#include <unistd.h>

std::string get_working_path() {
  char temp[MAXPATHLEN];
  return (getcwd(temp, MAXPATHLEN) ? std::string(temp) : std::string(""));
}

// fails on travis-ci with llvm
// shouldn't touch the fs from the c++ layer
TEST(Orientation, DISABLED_read_karney_datasets) {
  // todo: unzip data files
  // fill in data structure instead of stream
  std::cout << "TEST read_karney_datasets" << std::endl;

  // std::string s;
  // while(in >> s) std::cout << s << endl;
  // todo: assuming run from project root dir
  auto datfile =
      get_working_path() + "/" + "data/orientations/karney/c48u1.grid.gz";
  std::cerr << "TEST read_karney_datasets: " << datfile << std::endl;
  auto tuple = read_karney_orientation_file(datfile);
  Eigen::MatrixXd quats = std::get<0>(tuple);
  Eigen::VectorXd cover = std::get<1>(tuple);

  ASSERT_EQ(cover.size(), 24);
  // for(int i = 0; i < quats.rows(); ++i){
  // std::cout << quats(i,0) << " " << quats(i,1) << " " << quats(i,2) << " " <<
  // quats(i,3) << " " << cover(i) << std::endl;
  // }
  ASSERT_NEAR(quats(0, 0), 1, 0.0001);
  ASSERT_NEAR(quats(0, 1), 0, 0.0001);
  ASSERT_NEAR(quats(0, 2), 0, 0.0001);
  ASSERT_NEAR(quats(0, 3), 0, 0.0001);

  ASSERT_NEAR(quats(1, 0), 0, 0.0001);
  ASSERT_NEAR(quats(1, 1), 1, 0.0001);
  ASSERT_NEAR(quats(1, 2), 0, 0.0001);
  ASSERT_NEAR(quats(1, 3), 0, 0.0001);

  ASSERT_NEAR(quats(2, 0), 0, 0.0001);
  ASSERT_NEAR(quats(2, 1), 0, 0.0001);
  ASSERT_NEAR(quats(2, 2), 1, 0.0001);
  ASSERT_NEAR(quats(2, 3), 0, 0.0001);

  ASSERT_NEAR(quats(3, 0), 0, 0.0001);
  ASSERT_NEAR(quats(3, 1), 0, 0.0001);
  ASSERT_NEAR(quats(3, 2), 0, 0.0001);
  ASSERT_NEAR(quats(3, 3), 1, 0.0001);

  ASSERT_NEAR(quats(23, 2), 0.707107, 0.0001);

  // ASSERT_TRUE(false);
}
