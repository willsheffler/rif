#include <gtest/gtest.h>
#include <iostream>
#include <sampling/orientations.hpp>
#include <Eigen/Dense>

TEST(Orientation, read_karney_datasets){
    // todo: unzip data files
    // fill in data structure instead of stream
    std::cout << "TEST read_karney_datasets" << std::endl;

    // std::string s;
    // while(in >> s) std::cout << s << endl;     
    // assuming run from cmake tmp dir???
    auto tuple = read_karney_orientation_file(
        "data/orientations/karney/c48u1.grid.gz");
    Eigen::MatrixXd quats = std::get<0>(tuple);
    Eigen::VectorXd cover = std::get<1>(tuple);

    ASSERT_EQ(cover.size(),24);
    // for(int i = 0; i < quats.rows(); ++i){
        // std::cout << quats(i,0) << " " << quats(i,1) << " " << quats(i,2) << " " << quats(i,3) << " " << cover(i) << std::endl;
    // }
    ASSERT_NEAR( quats(0,0), 1, 0.0001 );
    ASSERT_NEAR( quats(0,1), 0, 0.0001 );
    ASSERT_NEAR( quats(0,2), 0, 0.0001 );
    ASSERT_NEAR( quats(0,3), 0, 0.0001 );

    ASSERT_NEAR( quats(1,0), 0, 0.0001 );
    ASSERT_NEAR( quats(1,1), 1, 0.0001 );
    ASSERT_NEAR( quats(1,2), 0, 0.0001 );
    ASSERT_NEAR( quats(1,3), 0, 0.0001 );

    ASSERT_NEAR( quats(2,0), 0, 0.0001 );
    ASSERT_NEAR( quats(2,1), 0, 0.0001 );
    ASSERT_NEAR( quats(2,2), 1, 0.0001 );
    ASSERT_NEAR( quats(2,3), 0, 0.0001 );

    ASSERT_NEAR( quats(3,0), 0, 0.0001 );
    ASSERT_NEAR( quats(3,1), 0, 0.0001 );
    ASSERT_NEAR( quats(3,2), 0, 0.0001 );
    ASSERT_NEAR( quats(3,3), 1, 0.0001 );

    ASSERT_NEAR( quats(23,2), 0.707107, 0.0001);

    // ASSERT_TRUE(false);

}   
