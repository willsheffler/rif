#include <gtest/gtest.h>

#include "objective/storage/RotamerScores.hpp"

#include <random>

namespace scheme {
namespace objective {
namespace storage {
namespace rstest {

using std::cout;
using std::endl;

TEST(RotamerScore, test_RotamerScore) {
  int bits = RotamerScore<>::RotamerBits;
  int mask = RotamerScore<>::RotamerMask;
  float divisor = RotamerScore<>::divisor();

  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform(127.0 / divisor, 0.0);
  std::uniform_int_distribution<> rand_rot(0, 511);

  ASSERT_EQ(9, bits);
  ASSERT_EQ(511, mask);

  for (int i = 0; i < 10000; ++i) {
    float score = uniform(rng);
    int rot = rand_rot(rng);
    RotamerScore<> rotscore(rot, score);
    ASSERT_EQ(rotscore.rotamer(), rot);
    // cout << rotscore.score() << " " << score << endl;
    ASSERT_LE(fabs(rotscore.score() - score), -1.0 / divisor);
  }
}

TEST(RotamerScores, test_store_1) {
  ASSERT_EQ(sizeof(RotamerScores<1>), 2);

  RotamerScores<1> rs;
  rs.add_rotamer(0, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), -1.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  rs.add_rotamer(0, -0.9);
  ASSERT_FLOAT_EQ(rs.score(0), -1.0);  // not lower
  ASSERT_EQ(rs.rotamer(0), 0);
  rs.add_rotamer(0, -1.125);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);  // is lower
  ASSERT_EQ(rs.rotamer(0), 0);

  rs.add_rotamer(1, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  // ASSERT_FLOAT_EQ( rs.score( 1 ),  0.0 );
  rs.add_rotamer(1, -1.25);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.25 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 1);
  rs.add_rotamer(0, -1.50);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.5 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);

  rs.add_rotamer(7, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-2.5 * -13.0) / -13.0);
}

TEST(RotamerScores, test_store_2) {
  ASSERT_EQ(sizeof(RotamerScores<2>), 4);
  RotamerScores<2> rs;
  rs.add_rotamer(0, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);
  rs.add_rotamer(0, -0.9);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);  // not lower
  rs.add_rotamer(0, -1.125);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);  // is lower

  rs.add_rotamer(1, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.0 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 1);

  rs.add_rotamer(1, -1.25);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 1);

  rs.add_rotamer(0, -1.50);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 1);

  rs.add_rotamer(7, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-2.50 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 7);
}

TEST(RotamerScores, test_store_3) {
  ASSERT_EQ(sizeof(RotamerScores<3>), 6);

  RotamerScores<3> rs;
  rs.add_rotamer(0, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);
  rs.add_rotamer(0, -0.9);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);  // not lower
  ASSERT_FLOAT_EQ(rs.score(1), int(0.0 * -13.0) / -13.0);   // not lower
  rs.add_rotamer(0, -1.125);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);  // is lower
  ASSERT_FLOAT_EQ(rs.score(1), int(0.0 * -13.0) / -13.0);     // not lower

  rs.add_rotamer(1, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.0 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);
  rs.add_rotamer(1, -1.25);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);
  rs.add_rotamer(0, -1.50);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);

  rs.add_rotamer(7, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(-2.5 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 1);
  ASSERT_EQ(rs.rotamer(2), 7);

  rs.add_rotamer(4, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-2.5 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(-2.5 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 4);
  ASSERT_EQ(rs.rotamer(2), 7);
}

TEST(RotamerScores, test_sort) {
  ASSERT_EQ(sizeof(RotamerScores<12>), 24);
  RotamerScores<12> rs;
  rs.add_rotamer(0, -0.1);
  rs.add_rotamer(1, -3.0);
  rs.add_rotamer(2, -1.0);
  rs.add_rotamer(3, -2.0);
  rs.add_rotamer(4, -4.0);

  // std::cout << "nosort: " << rs << std::endl;
  rs.sort_rotamers();
  // std::cout << "sorted: " << rs << std::endl;
  ASSERT_EQ(rs.rotamer(0), 4);
  ASSERT_EQ(rs.rotamer(1), 1);
  ASSERT_EQ(rs.rotamer(2), 3);
  ASSERT_EQ(rs.rotamer(3), 2);
  ASSERT_EQ(rs.rotamer(4), 0);
  ASSERT_TRUE(rs.empty(5));
}

TEST(RotamerScores, test_sat_data) {
  ASSERT_EQ(sizeof(RotamerScores<3, RotamerScoreSat<>>), 12);
  RotamerScores<3, RotamerScoreSat<>> rs;

  rs.add_rotamer(0, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);
  rs.add_rotamer(0, -0.9);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.0 * -13.0) / -13.0);  // not lower
  ASSERT_FLOAT_EQ(rs.score(1), int(0.0 * -13.0) / -13.0);   // not lower
  rs.add_rotamer(0, -1.125);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);  // is lower
  ASSERT_FLOAT_EQ(rs.score(1), int(0.0 * -13.0) / -13.0);     // not lower

  rs.add_rotamer(1, -1.0);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.0 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);
  rs.add_rotamer(1, -1.25);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.125 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);
  rs.add_rotamer(0, -1.50);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(0.0 * -13.0) / -13.0);

  rs.add_rotamer(7, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-1.25 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(-2.5 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 1);
  ASSERT_EQ(rs.rotamer(2), 7);

  rs.add_rotamer(4, -2.5);
  ASSERT_FLOAT_EQ(rs.score(0), int(-1.50 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(1), int(-2.5 * -13.0) / -13.0);
  ASSERT_FLOAT_EQ(rs.score(2), int(-2.5 * -13.0) / -13.0);
  ASSERT_EQ(rs.rotamer(0), 0);
  ASSERT_EQ(rs.rotamer(1), 4);
  ASSERT_EQ(rs.rotamer(2), 7);
}

TEST(RotamerScores, test_sort_sat) {
  typedef RotamerScores<14, RotamerScoreSat<>> RS;
  ASSERT_EQ(sizeof(RS), 56);
  RS rs;
  rs.add_rotamer(0, -0.1, 3, 7);
  rs.add_rotamer(1, -3.0, -1, -1);
  rs.add_rotamer(2, -1.0, 0, -1);
  rs.add_rotamer(3, -2.0, 1, 2);
  rs.add_rotamer(4, -4.0, 9, 5);

  // std::cout << "nosort: " << rs << std::endl;
  rs.sort_rotamers();
  std::cout << "sorted: " << rs << std::endl;
  ASSERT_EQ(rs.rotamer(0), 4);
  ASSERT_EQ(rs.rotamer(1), 3);  // was 1
  ASSERT_EQ(rs.rotamer(2), 0);  // was 3
  ASSERT_EQ(rs.rotamer(3), 2);
  ASSERT_EQ(rs.rotamer(4), 1);  // was 0
  ASSERT_TRUE(rs.empty(5));

  RS rs2(rs);
  rs2.add_rotamer(2, -1.0, 0, -1);
  ASSERT_EQ(rs, rs2);

  ASSERT_EQ(rs2.rotscores_[3].sat_data_[1].data_, 255);
  rs2.add_rotamer(2, -1.0, 1, -1);
  ASSERT_NE(rs, rs2);
  ASSERT_EQ(rs2.rotscores_[3].sat_data_[1].data_, 1);

  ASSERT_EQ(rs2.rotscores_[4].sat_data_[0].data_, 255);
  ASSERT_EQ(rs2.rotscores_[4].sat_data_[1].data_, 255);
  rs2.add_rotamer(1, -3.1, 2, 3);
  // cout << rs2 << endl;
  ASSERT_LE(rs2.score(4), -3.0);
  ASSERT_EQ(rs2.rotscores_[4].sat_data_[0].data_, 2);
  ASSERT_EQ(rs2.rotscores_[4].sat_data_[1].data_, 3);

  rs2.add_rotamer(5, -9, 3);
  RS rs3(rs2);
  rs3.add_rotamer(5, -8);
  ASSERT_EQ(rs3, rs2);
}
// TEST( RotamerScores, test_store_4 ){

// 	RotamerScores<4> rs;

// 	ASSERT_EQ( rs.name(), "RotamerScores< 4, 2, 1024 >" );

// 	rs.add_rotamer( 0, -1.0 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.0 );
// 	rs.add_rotamer( 0, -0.9 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.0 ); // not lower
// 	rs.add_rotamer( 0, -1.125 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.125 ); // is lower

// 	rs.add_rotamer( 1, -1.0 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.125 );
// 	ASSERT_FLOAT_EQ( rs.score( 1 ), -1.0 );
// 	rs.add_rotamer( 1, -1.25 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.125 );
// 	ASSERT_FLOAT_EQ( rs.score( 1 ), -1.25 );
// 	rs.add_rotamer( 0, -1.50 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.50 );
// 	ASSERT_FLOAT_EQ( rs.score( 1 ), -1.25 );

// 	rs.add_rotamer( 7, -2.5 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.50 );
// 	ASSERT_FLOAT_EQ( rs.score( 1 ), -1.25 );
// 	ASSERT_FLOAT_EQ( rs.score( 7 ), -2.5 );

// 	rs.add_rotamer( 4, -2.5 );
// 	ASSERT_FLOAT_EQ( rs.score( 0 ), -1.50 );
// 	ASSERT_FLOAT_EQ( rs.score( 1 ), -1.25 );
// 	ASSERT_FLOAT_EQ( rs.score( 7 ), -2.5 );
// 	ASSERT_FLOAT_EQ( rs.score( 4 ), -2.5 );

// }

// TEST( RotamerScores, sort ){
// 	std::mt19937 rng((unsigned int)time(0));
// 	std::uniform_real_distribution<> uniform;

// 	for( int iter = 0; iter < 100; ++iter ){

// 		RotamerScores<10> rs;
// 		for( int i = 0; i < rs.size(); ++i ){
// 			rs.add_rotamer( i, uniform(rng) );
// 		}
// 		rs.sort_rotamers();

// 		for( int i = 1; i < rs.size(); ++i ){
// 			ASSERT_LE( rs.scores_[i-1], rs.scores_[i] );
// 		}

// 	}

// }
}
}
}
}
