#include <gtest/gtest.h>

#include <random>

#include "riflib/numeric/cube_to_sphere.hpp"

namespace scheme{ namespace numeric { namespace cube_to_sphere_test {


TEST( cube_to_sphere, test_permute ){
    std::mt19937 rng;
    std::normal_distribution<> rnorm;
    std::uniform_int_distribution<> randint(0,5);
    for(int i = 0; i < 100000; ++i){
        int facenum = randint(rng);
        float x, x0, y, y0, z, z0;
        x0 = x = rnorm(rng);
        y0 = y = rnorm(rng);
        z0 = z = rnorm(rng);
        permute_cube_face_xyz(facenum, x, y, z);
        inverse_permute_cube_face_xyz(facenum, x, y, z);
        ASSERT_FLOAT_EQ(x0, x);
        ASSERT_FLOAT_EQ(y0, y);
        ASSERT_FLOAT_EQ(z0, z);
    }
}

}}}