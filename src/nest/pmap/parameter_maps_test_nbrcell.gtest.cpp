#include <gtest/gtest.h>
#include <boost/assign/std/vector.hpp>  // for 'operator+=()'
#include <boost/foreach.hpp>
#include <random>
#include "nest/NEST.hpp"
#include "nest/pmap/ScaleMap.hpp"
#include "nest/pmap/UnitMap.hpp"

namespace rif {
namespace nest {
namespace pmap {

using std::cout;
using std::endl;

TEST(ParamMap, UnitMap_get_neighboring_cells) {
  {
    typedef UnitMap<2> MapType;
    MapType umap(1);
    std::vector<size_t> cnb;
    std::back_insert_iterator<std::vector<size_t>> biter(cnb);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5), 0, 0.4, biter);
    EXPECT_EQ(cnb.size(), 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(-0.5, 0.5), 0, 0.4, biter);
    EXPECT_EQ(cnb.size(), 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(-0.39, 0.5), 0, 0.4, biter);
    EXPECT_EQ(cnb.size(), 1);
    EXPECT_EQ(cnb[0], 0);

    // cout << "size: "<< cnb.size() << endl;
    // BOOST_FOREACH(size_t i,cnb) cout << "val " << i << std::endl;
  }
  {
    typedef UnitMap<2> MapType;
    MapType umap(100);
    std::vector<size_t> cnb;
    std::back_insert_iterator<std::vector<size_t>> biter(cnb);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 1);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(10.5, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 10);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(10.3, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 2);
    ASSERT_EQ(cnb[0], 9);
    ASSERT_EQ(cnb[1], 10);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(0.3, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(0.3, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(101.0, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(100.3, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 99);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(-0.5, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 0);

    cnb.clear();
    umap.get_neighboring_cells(MapType::ValueType(-0.3, 0.5), 0, 0.4, biter);
    ASSERT_EQ(cnb.size(), 1);
    ASSERT_EQ(cnb[0], 0);
  }
}

TEST(ScaleMap, ScaleMap_cellindices) {
  util::SimpleArray<1, size_t> bs(3);
  ScaleMap<1> m1(bs);
  for (size_t i = 0; i < 3; ++i)
    ASSERT_EQ(i, m1.indices_to_cellindex(m1.cellindex_to_indices(i)));
  ScaleMap<2> m2(util::SimpleArray<2, size_t>(3, 4));
  for (size_t i = 0; i < 3 * 4; ++i)
    ASSERT_EQ(i, m2.indices_to_cellindex(m2.cellindex_to_indices(i)));
  ScaleMap<3> m3(util::SimpleArray<3, size_t>(3, 4, 8));
  for (size_t i = 0; i < 3 * 4 * 8; ++i)
    ASSERT_EQ(i, m3.indices_to_cellindex(m3.cellindex_to_indices(i)));
  ScaleMap<4> m4(util::SimpleArray<4, size_t>(3, 4, 8, 17));
  for (size_t i = 0; i < 3 * 4 * 8 * 17; ++i)
    ASSERT_EQ(i, m4.indices_to_cellindex(m4.cellindex_to_indices(i)));
}

TEST(ScaleMap, ScaleMap_get_neighboring_cells_244) {
  typedef ScaleMap<2> MapType;
  MapType umap(MapType::Params(0, 0), MapType::Params(4, 4),
               MapType::Indices(4, 4));
  std::vector<size_t> cnb;
  std::back_insert_iterator<std::vector<size_t>> biter(cnb);
  int i;

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5), 0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5), 0, 0.6, biter);
  ASSERT_EQ(cnb.size(), 6);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 2);
  ASSERT_EQ(cnb[i++], 4);
  ASSERT_EQ(cnb[i++], 5);
  ASSERT_EQ(cnb[i++], 6);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 1.5), 0, 0.6, biter);
  ASSERT_EQ(cnb.size(), 9);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 2);
  ASSERT_EQ(cnb[i++], 4);
  ASSERT_EQ(cnb[i++], 5);
  ASSERT_EQ(cnb[i++], 6);
  ASSERT_EQ(cnb[i++], 8);
  ASSERT_EQ(cnb[i++], 9);
  ASSERT_EQ(cnb[i++], 10);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.5, 4.5), 0, 0.6, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 15);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.5, 1.5), 0, 0.2, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 6);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.1, 1.5), 0, 0.2, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 5);
  ASSERT_EQ(cnb[i++], 6);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.1, 1.9), 0, 0.2, biter);
  ASSERT_EQ(cnb.size(), 4);
  ASSERT_EQ(cnb[i++], 5);
  ASSERT_EQ(cnb[i++], 6);
  ASSERT_EQ(cnb[i++], 9);
  ASSERT_EQ(cnb[i++], 10);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.4, 4.4), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 0);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.4, 4.2), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 0);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.2, 4.4), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 0);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.2, 4.2), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 15);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(4.2, 3.2), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 11);
  ASSERT_EQ(cnb[i++], 15);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(3.2, 4.2), 0, 0.3, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 14);
  ASSERT_EQ(cnb[i++], 15);
}

TEST(ScaleMap, ScaleMap_get_neighboring_cells_433) {
  typedef ScaleMap<4> MapType;
  MapType umap(MapType::Params(0, 0, 0, 0), MapType::Params(3, 3, 3, 3),
               MapType::Indices(3, 3, 3, 3));
  std::vector<size_t> cnb;
  std::back_insert_iterator<std::vector<size_t>> biter(cnb);
  int i;

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5, 0.5, 0.5), 0, 0.4,
                             biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.5, 0.5, 0.5), 0, 0.4,
                             biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 0.5), 0, 0.4,
                             biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 1.3), 0, 0.4,
                             biter);
  ASSERT_EQ(cnb.size(), 4);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 27);
  ASSERT_EQ(cnb[i++], 28);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.3, 1.8, 2.1, 1.8), 0, 0.4,
                             biter);
  ASSERT_EQ(cnb.size(), 16);
  ASSERT_EQ(cnb[i++], 40);
  ASSERT_EQ(cnb[i++], 41);
  ASSERT_EQ(cnb[i++], 43);
  ASSERT_EQ(cnb[i++], 44);
  ASSERT_EQ(cnb[i++], 49);
  ASSERT_EQ(cnb[i++], 50);
  ASSERT_EQ(cnb[i++], 52);
  ASSERT_EQ(cnb[i++], 53);
  ASSERT_EQ(cnb[i++], 67);
  ASSERT_EQ(cnb[i++], 68);
  ASSERT_EQ(cnb[i++], 70);
  ASSERT_EQ(cnb[i++], 71);
  ASSERT_EQ(cnb[i++], 76);
  ASSERT_EQ(cnb[i++], 77);
  ASSERT_EQ(cnb[i++], 79);
  ASSERT_EQ(cnb[i++], 80);

  // cout << "size: "<< cnb.size() << endl;
  // BOOST_FOREACH(size_t i,cnb) cout << "val " << i << " : " <<
  // umap.cellindex_to_indices(i).transpose() << std::endl;
}

TEST(ScaleMap, ScaleMap_get_neighboring_cells_433_m1) {
  typedef ScaleMap<4> MapType;
  MapType umap(MapType::Params(0, 0, 0, 0) - 1.0,
               MapType::Params(3, 3, 3, 3) - 1.0, MapType::Indices(3, 3, 3, 3));
  std::vector<size_t> cnb;
  std::back_insert_iterator<std::vector<size_t>> biter(cnb);
  int i;

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5, 0.5, 0.5) - 1.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.5, 0.5, 0.5) - 1.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 0.5) - 1.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 1.3) - 1.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 4);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 27);
  ASSERT_EQ(cnb[i++], 28);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.3, 1.8, 2.1, 1.8) - 1.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 16);
  ASSERT_EQ(cnb[i++], 40);
  ASSERT_EQ(cnb[i++], 41);
  ASSERT_EQ(cnb[i++], 43);
  ASSERT_EQ(cnb[i++], 44);
  ASSERT_EQ(cnb[i++], 49);
  ASSERT_EQ(cnb[i++], 50);
  ASSERT_EQ(cnb[i++], 52);
  ASSERT_EQ(cnb[i++], 53);
  ASSERT_EQ(cnb[i++], 67);
  ASSERT_EQ(cnb[i++], 68);
  ASSERT_EQ(cnb[i++], 70);
  ASSERT_EQ(cnb[i++], 71);
  ASSERT_EQ(cnb[i++], 76);
  ASSERT_EQ(cnb[i++], 77);
  ASSERT_EQ(cnb[i++], 79);
  ASSERT_EQ(cnb[i++], 80);

  // cout << "size: "<< cnb.size() << endl;
  // BOOST_FOREACH(size_t i,cnb) cout << "val " << i << " : " <<
  // umap.cellindex_to_indices(i).transpose() << std::endl;
}

TEST(ScaleMap, ScaleMap_get_neighboring_cells_433_o2) {
  typedef ScaleMap<4> MapType;
  MapType umap(MapType::Params(0, 0, 0, 0) / 2.0,
               MapType::Params(3, 3, 3, 3) / 2.0, MapType::Indices(3, 3, 3, 3));
  std::vector<size_t> cnb;
  std::back_insert_iterator<std::vector<size_t>> biter(cnb);
  int i;

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.5, 0.5, 0.5, 0.5) / 2.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.5, 0.5, 0.5) / 2.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 0.5) / 2.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 1.3) / 2.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 4);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 27);
  ASSERT_EQ(cnb[i++], 28);

  cnb.clear();
  i = 0;
  umap.get_neighboring_cells(MapType::ValueType(2.3, 1.8, 2.1, 1.8) / 2.0, 0,
                             0.4, biter);
  ASSERT_EQ(cnb.size(), 16);
  ASSERT_EQ(cnb[i++], 40);
  ASSERT_EQ(cnb[i++], 41);
  ASSERT_EQ(cnb[i++], 43);
  ASSERT_EQ(cnb[i++], 44);
  ASSERT_EQ(cnb[i++], 49);
  ASSERT_EQ(cnb[i++], 50);
  ASSERT_EQ(cnb[i++], 52);
  ASSERT_EQ(cnb[i++], 53);
  ASSERT_EQ(cnb[i++], 67);
  ASSERT_EQ(cnb[i++], 68);
  ASSERT_EQ(cnb[i++], 70);
  ASSERT_EQ(cnb[i++], 71);
  ASSERT_EQ(cnb[i++], 76);
  ASSERT_EQ(cnb[i++], 77);
  ASSERT_EQ(cnb[i++], 79);
  ASSERT_EQ(cnb[i++], 80);

  // cout << "size: "<< cnb.size() << endl;
  // BOOST_FOREACH(size_t i,cnb) cout << "val " << i << " : " <<
  // umap.cellindex_to_indices(i).transpose() << std::endl;
}

TEST(ScaleMap, ScaleMap_get_neighboring_cells_433_o5m1) {
  typedef ScaleMap<4> MapType;
  MapType smap(MapType::Params(0, 0, 0, 0) / 5.0 - 1.0,
               MapType::Params(3, 3, 3, 3) / 5.0 - 1.0,
               MapType::Indices(3, 3, 3, 3));
  std::vector<size_t> cnb;
  std::back_insert_iterator<std::vector<size_t>> biter(cnb);
  int i;

  cnb.clear();
  i = 0;
  smap.get_neighboring_cells(MapType::ValueType(1.5, 0.5, 0.5, 0.5) / 5.0 - 1.0,
                             0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 1);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  smap.get_neighboring_cells(MapType::ValueType(1.3, 0.5, 0.5, 0.5) / 5.0 - 1.0,
                             0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  smap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 0.5) / 5.0 - 1.0,
                             0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 2);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);

  cnb.clear();
  i = 0;
  smap.get_neighboring_cells(MapType::ValueType(1.3, 0.1, 0.1, 1.3) / 5.0 - 1.0,
                             0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 4);
  ASSERT_EQ(cnb[i++], 0);
  ASSERT_EQ(cnb[i++], 1);
  ASSERT_EQ(cnb[i++], 27);
  ASSERT_EQ(cnb[i++], 28);

  cnb.clear();
  i = 0;
  smap.get_neighboring_cells(MapType::ValueType(2.3, 1.8, 2.1, 1.8) / 5.0 - 1.0,
                             0, 0.4, biter);
  ASSERT_EQ(cnb.size(), 16);
  ASSERT_EQ(cnb[i++], 40);
  ASSERT_EQ(cnb[i++], 41);
  ASSERT_EQ(cnb[i++], 43);
  ASSERT_EQ(cnb[i++], 44);
  ASSERT_EQ(cnb[i++], 49);
  ASSERT_EQ(cnb[i++], 50);
  ASSERT_EQ(cnb[i++], 52);
  ASSERT_EQ(cnb[i++], 53);
  ASSERT_EQ(cnb[i++], 67);
  ASSERT_EQ(cnb[i++], 68);
  ASSERT_EQ(cnb[i++], 70);
  ASSERT_EQ(cnb[i++], 71);
  ASSERT_EQ(cnb[i++], 76);
  ASSERT_EQ(cnb[i++], 77);
  ASSERT_EQ(cnb[i++], 79);
  ASSERT_EQ(cnb[i++], 80);

  // cout << "size: "<< cnb.size() << endl;
  // BOOST_FOREACH(size_t i,cnb) cout << "val " << i << " : " <<
  // smap.cellindex_to_indices(i).transpose() << std::endl;
}
}
}
}
