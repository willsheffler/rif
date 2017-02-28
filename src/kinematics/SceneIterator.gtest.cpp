#include <gtest/gtest.h>

#include "actor/ActorConcept_io.hpp"
#include "kinematics/Scene_io.hpp"
#include "numeric/X1dim.hpp"

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>

namespace rif {
namespace kinematics {
namespace test {

using std::cout;
using std::endl;
using boost::tie;
using numeric::X1dim;

typedef actor::ActorConcept<X1dim, int> ADI;
typedef actor::ActorConcept<X1dim, char> ADC;
typedef m::vector<ADI, ADC> Actors;

typedef Scene<Actors, X1dim, size_t> Scene;
typedef size_t Index;
typedef std::pair<size_t, size_t> Index2;
typedef std::pair<Index2, Index2> Index4;

TEST(SceneIterator, interactions_emptyish_onebody) {
  Scene s;

  Scene::iter_type<ADI>::type adibeg, adiend;
  Scene::iter_type<ADC>::type adcbeg, adcend;

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  // cout << *adibeg << " " << *adiend << endl;
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 1));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(*adcbeg++, Index2(1, 0));
  ASSERT_EQ(adcbeg, adcend);

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  // ADI adi = s.get_interaction<ADI>(*adibeg);
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(0, 0));
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(0, 1));
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(s.get_interaction<ADC>(*adcbeg++), ADC(1, '0'));
  ASSERT_EQ(adcbeg, adcend);

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}

TEST(SceneIterator,
     interactions_emptyish_onebody_sym_interaction) {  // one-body should be all
                                                       // asym, so no effect
  Scene s;
  s.add_symframe(10.0);  // testing this has no effect on one-body interactions
  s.add_symframe(20.0);  // testing this has no effect on one-body interactions
  Scene::iter_type<ADI>::type adibeg, adiend;
  Scene::iter_type<ADC>::type adcbeg, adcend;

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 1));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(*adcbeg++, Index2(1, 0));
  ASSERT_EQ(adcbeg, adcend);

  tie(adibeg, adiend) = s.get_interactions<ADI>();
  tie(adcbeg, adcend) = s.get_interactions<ADC>();
  // ADI adi = s.get_interaction<ADI>(*adibeg);
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(0, 0));
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(0, 1));
  ASSERT_EQ(s.get_interaction<ADI>(*adibeg++), ADI(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(s.get_interaction<ADC>(*adcbeg++), ADC(1, '0'));
  ASSERT_EQ(adcbeg, adcend);

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}
TEST(SceneIterator, interactions_emptyish_onebody_asym_getactor) {
  Scene s;
  s.add_symframe(10.0);
  s.add_symframe(20.0);
  Scene::iter_type<ADI>::type adibeg, adiend;
  Scene::iter_type<ADC>::type adcbeg, adcend;

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 1));

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(adibeg, adiend) = s.get_actors_placeholder_asym<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder_asym<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(*adcbeg++, Index2(1, 0));
  ASSERT_EQ(adcbeg, adcend);

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}
TEST(SceneIterator, interactions_emptyish_onebody_sym_getactor_placeholder) {
  Scene s;
  s.add_symframe(10.0);
  s.add_symframe(20.0);
  SceneIter1B<Scene, ADI, Symmetric, PlaceHolder> adibeg, adiend;
  SceneIter1B<Scene, ADC, Symmetric, PlaceHolder> adcbeg, adcend;

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 1));

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(2, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 1));
  ASSERT_EQ(*adibeg++, Index2(4, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 1));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 1));
  ASSERT_EQ(*adibeg++, Index2(3, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 1));
  ASSERT_EQ(*adibeg++, Index2(5, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(adibeg, adiend) = s.get_actors_placeholder<ADI>();
  tie(adcbeg, adcend) = s.get_actors_placeholder<ADC>();
  ASSERT_EQ(*adibeg++, Index2(0, 0));
  ASSERT_EQ(*adibeg++, Index2(0, 1));
  ASSERT_EQ(*adibeg++, Index2(1, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 0));
  ASSERT_EQ(*adibeg++, Index2(2, 1));
  ASSERT_EQ(*adibeg++, Index2(3, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 0));
  ASSERT_EQ(*adibeg++, Index2(4, 1));
  ASSERT_EQ(*adibeg++, Index2(5, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(*adcbeg++, Index2(1, 0));
  ASSERT_EQ(*adcbeg++, Index2(3, 0));
  ASSERT_EQ(*adcbeg++, Index2(5, 0));
  ASSERT_EQ(adcbeg, adcend);

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}

TEST(SceneIterator, interactions_emptyish_onebody_sym_getactor) {
  Scene s;
  s.add_symframe(10.0);
  s.add_symframe(20.0);
  SceneIter1B<Scene, ADI, Symmetric, ActorCopy> adibeg, adiend;
  SceneIter1B<Scene, ADC, Symmetric, ActorCopy> adcbeg, adcend;

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();
  ASSERT_EQ(s.bodies_.size(), 1);
  ASSERT_EQ(s.positions_.size(), 1);

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.add_body();

  ASSERT_EQ(s.bodies_.size(), 2);
  ASSERT_EQ(s.positions_.size(), 2);

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  EXPECT_EQ(*adibeg++, ADI(00, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 1));

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  EXPECT_EQ(*adibeg++, ADI(00, 0));
  EXPECT_EQ(*adibeg++, ADI(00, 1));
  EXPECT_EQ(*adibeg++, ADI(10, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 1));
  EXPECT_EQ(*adibeg++, ADI(20, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 1));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  EXPECT_EQ(*adibeg++, ADI(00, 0));
  EXPECT_EQ(*adibeg++, ADI(00, 1));
  EXPECT_EQ(*adibeg++, ADI(01, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 1));
  EXPECT_EQ(*adibeg++, ADI(11, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 1));
  EXPECT_EQ(*adibeg++, ADI(21, 0));
  ASSERT_EQ(adibeg, adiend);
  ASSERT_EQ(adcbeg, adcend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(adibeg, adiend) = s.get_actors<ADI>();
  tie(adcbeg, adcend) = s.get_actors<ADC>();
  EXPECT_EQ(*adibeg++, ADI(00, 0));
  EXPECT_EQ(*adibeg++, ADI(00, 1));
  EXPECT_EQ(*adibeg++, ADI(01, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 0));
  EXPECT_EQ(*adibeg++, ADI(10, 1));
  EXPECT_EQ(*adibeg++, ADI(11, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 0));
  EXPECT_EQ(*adibeg++, ADI(20, 1));
  EXPECT_EQ(*adibeg++, ADI(21, 0));
  ASSERT_EQ(adibeg, adiend);
  EXPECT_EQ(*adcbeg++, ADC(01, '0'));
  EXPECT_EQ(*adcbeg++, ADC(11, '0'));
  EXPECT_EQ(*adcbeg++, ADC(21, '0'));
  ASSERT_EQ(adcbeg, adcend);

  std::ostringstream out;
  BOOST_FOREACH (ADI a, s.get_actors<ADI>())
    out << a << endl;
  BOOST_FOREACH (ADC a, s.get_actors<ADC>())
    out << a << endl;
  out << s << endl;
  // cout << out.str();
}

typedef std::pair<ADI, ADC> IIC;
typedef std::pair<ADC, ADI> ICI;
typedef std::pair<ADI, ADI> III;
typedef std::pair<ADC, ADC> ICC;

TEST(SceneIterator, interactions_emptyish_twobody) {
  using std::make_pair;

  Scene s;
  Scene::iter_type<IIC>::type iicbeg, iicend;
  Scene::iter_type<ICI>::type icibeg, iciend;
  Scene::iter_type<III>::type iiibeg, iiiend;
  Scene::iter_type<ICC>::type iccbeg, iccend;

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.add_body();

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.add_body();

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(0).add_actor(
      ADC(0, '0'));  // still empty, because intra-body

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_EQ(iicbeg, iicend);
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_EQ(icibeg, iciend);
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  ASSERT_EQ(iiibeg, iiiend);
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '1'));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_EQ(iicbeg, iicend);
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(1, 0)));
  ASSERT_EQ(icibeg, iciend);
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  ASSERT_EQ(iiibeg, iiiend);
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_EQ(iicbeg, iicend);
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(1, 0)));
  ASSERT_EQ(icibeg, iciend);
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(*iiibeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_EQ(iiibeg, iiiend);
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 1));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(*iicbeg++, make_pair(Index2(1, 0), Index2(1, 0)));
  ASSERT_EQ(iicbeg, iicend);
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(*icibeg++, make_pair(Index2(1, 0), Index2(1, 0)));
  ASSERT_EQ(icibeg, iciend);
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(*iiibeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(*iiibeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_EQ(iiibeg, iiiend);
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 0)));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(*iccbeg++, make_pair(Index2(0, 1), Index2(0, 1)));
  ASSERT_EQ(iccbeg, iccend);

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction<IIC>(*iicbeg++),
            make_pair(ADI(0, 0), ADC(1, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction<IIC>(*iicbeg++),
            make_pair(ADI(0, 0), ADC(1, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction<IIC>(*iicbeg++),
            make_pair(ADI(1, 0), ADC(0, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction<IIC>(*iicbeg++),
            make_pair(ADI(1, 1), ADC(0, '0')));
  ASSERT_EQ(iicbeg, iicend);
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction<ICI>(*icibeg++),
            make_pair(ADC(0, '0'), ADI(1, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction<ICI>(*icibeg++),
            make_pair(ADC(0, '0'), ADI(1, 1)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction<ICI>(*icibeg++),
            make_pair(ADC(1, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction<ICI>(*icibeg++),
            make_pair(ADC(1, '1'), ADI(0, 0)));
  ASSERT_EQ(icibeg, iciend);
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction<III>(*iiibeg++), make_pair(ADI(0, 0), ADI(1, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction<III>(*iiibeg++), make_pair(ADI(0, 0), ADI(1, 1)));
  ASSERT_EQ(iiibeg, iiiend);
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction<ICC>(*iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction<ICC>(*iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '1')));
  ASSERT_EQ(iccbeg, iccend);

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}

TEST(SceneIterator, interactions_emptyish_twobody_sym) {
  using std::make_pair;

  Scene s;
  s.add_symframe(10.0);
  Scene::iter_type<IIC>::type iicbeg, iicend;
  Scene::iter_type<ICI>::type icibeg, iciend;
  Scene::iter_type<III>::type iiibeg, iiiend;
  Scene::iter_type<ICC>::type iccbeg, iccend;

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.add_body();

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.add_body();

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ((*iiibeg++).first, Index2(0, 2));
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(0).add_actor(
      ADC(0, '0'));  // still empty, because intra-body
  // BOOST_FOREACH( ADI a, s.get_actors<ADI>() ) cout << a << endl;
  // BOOST_FOREACH( ADC a, s.get_actors<ADC>() ) cout << a << endl;
  // ADI(  0, 0 ) ADI(  0, 0 )
  // ADC(  0, 0 ) ADI( 10, 0 )
  //              ADC(  0, 0 )
  //              ADC( 10, 0 )

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(10, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(0, '0')));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(10, '0'), ADI(0, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(10, 0)));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(10, '0')));
  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();

  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(01, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(10, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(11, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(00, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(01, '0')));

  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(10, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '0'), ADI(0, 0)));

  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(10, 0)));

  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(11, '0')));

  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADC(1, '1'));

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();

  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(01, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(01, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(10, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(11, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(11, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(00, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(01, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(01, '1')));

  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '1'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '1'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(10, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '1'), ADI(0, 0)));

  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(10, 0)));

  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(11, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(11, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(11, '1')));

  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));
  // BOOST_FOREACH( ADI a, s.get_actors<ADI>() ) cout << a << endl;
  // BOOST_FOREACH( ADC a, s.get_actors<ADC>() ) cout << a << endl;

  // ADI( 0, 0 )
  // ADI( 1, 0 )
  // ADI( 10, 0 )
  // ADI( 11, 0 )

  // ADC( 0, 0 )
  // ADC( 1, 0 )
  // ADC( 1, 1 )
  // ADC( 10, 0 )
  // ADC( 11, 0 )
  // ADC( 11, 1 )

  tie(iicbeg, iicend) = s.get_interactions<IIC>();
  tie(icibeg, iciend) = s.get_interactions<ICI>();
  tie(iiibeg, iiiend) = s.get_interactions<III>();
  tie(iccbeg, iccend) = s.get_interactions<ICC>();

  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(1, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(1, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(10, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(11, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(0, 0), ADC(11, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(1, 0), ADC(0, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(1, 0), ADC(10, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(1, 0), ADC(11, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(1, 0), ADC(11, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(0, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(1, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(10, 0), ADC(1, '1')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(11, 0), ADC(0, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(11, 0), ADC(1, '0')));
  ASSERT_NE(iicbeg, iicend);
  ASSERT_EQ(s.get_interaction_absolute(iicbeg++),
            make_pair(ADI(11, 0), ADC(1, '1')));

  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(1, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(0, '0'), ADI(11, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '1'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '1'), ADI(10, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '0'), ADI(11, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(1, '1'), ADI(11, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(10, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(10, '0'), ADI(1, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '0'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '1'), ADI(0, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '0'), ADI(1, 0)));
  ASSERT_NE(icibeg, iciend);
  ASSERT_EQ(s.get_interaction_absolute(icibeg++),
            make_pair(ADC(11, '1'), ADI(1, 0)));

  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(1, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(10, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(0, 0), ADI(11, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(1, 0), ADI(10, 0)));
  ASSERT_NE(iiibeg, iiiend);
  ASSERT_EQ(s.get_interaction_absolute(iiibeg++),
            make_pair(ADI(1, 0), ADI(11, 0)));

  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(1, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(0, '0'), ADC(11, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(10, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '0'), ADC(11, '1')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(11, '0')));
  ASSERT_NE(iccbeg, iccend);
  ASSERT_EQ(s.get_interaction_absolute(iccbeg++),
            make_pair(ADC(1, '1'), ADC(11, '1')));

  ASSERT_EQ(iicbeg, iicend);
  ASSERT_EQ(icibeg, iciend);
  ASSERT_EQ(iiibeg, iiiend);
  ASSERT_EQ(iccbeg, iccend);

  // s.mutable_conformation_asym(1).add_actor(ADI(1,1));

  std::ostringstream out;
  out << s << endl;
  // cout << out.str();
}

TEST(SceneIterator, interactions_twobody_homo) {
  using std::make_pair;

  Scene s(5);
  Scene::iter_type<III>::type iter, end;

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));
  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));
  s.mutable_conformation_asym(2).add_actor(ADI(2, 0));
  s.mutable_conformation_asym(3).add_actor(ADI(3, 0));
  s.mutable_conformation_asym(4).add_actor(ADI(4, 0));

  tie(iter, end) = s.get_interactions<III>();
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 4));
  ASSERT_EQ(iter, end);
}

TEST(SceneIterator, interactions_twobody_hetero) {
  using std::make_pair;

  Scene s(4);
  Scene::iter_type<IIC>::type iter, end;

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));
  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));
  s.mutable_conformation_asym(2).add_actor(ADI(2, 0));
  s.mutable_conformation_asym(3).add_actor(ADI(3, 0));

  tie(iter, end) = s.get_interactions<IIC>();
  ASSERT_EQ(iter, end);

  s.mutable_conformation_asym(0).add_actor(ADC(0, '0'));
  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));
  s.mutable_conformation_asym(2).add_actor(ADC(2, '0'));
  s.mutable_conformation_asym(3).add_actor(ADC(3, '0'));

  // BOOST_FOREACH( Index4 i4, s.get_interactions<IIC>() ) cout << i4 << endl;

  tie(iter, end) = s.get_interactions<IIC>();
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 2));
  ASSERT_EQ(iter, end);
}

TEST(SceneIterator, interactions_twobody_homo_sym) {
  using std::make_pair;

  Scene s(3);
  s.add_symframe(10.0);
  ASSERT_EQ(s.nbodies_asym(), 3);
  ASSERT_EQ(s.nbodies(), 6);
  Scene::iter_type<III>::type iter, end;

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));
  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));
  s.mutable_conformation_asym(2).add_actor(ADI(2, 0));

  // BOOST_FOREACH( ADI a, s.get_actors<ADI>() ) cout << a << endl;
  // BOOST_FOREACH( ADC a, s.get_actors<ADC>() ) cout << a << endl;
  // ADI( 0, 0 ) ADI(  0, 0 )
  // ADI( 1, 0 ) ADI(  1, 0 )
  // ADI( 2, 0 ) ADI(  2, 0 )
  //             ADI( 10, 0 )
  //             ADI( 11, 0 )
  //             ADI( 12, 0 )

  tie(iter, end) = s.get_interactions<III>();
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 5));
  ASSERT_EQ(iter, end);
}

TEST(SceneIterator, interactions_twobody_hetero_sym) {
  using std::make_pair;

  Scene s(4);
  s.add_symframe(10.0);
  Scene::iter_type<IIC>::type iter, end;

  s.mutable_conformation_asym(0).add_actor(ADI(0, 0));
  s.mutable_conformation_asym(1).add_actor(ADI(1, 0));
  s.mutable_conformation_asym(2).add_actor(ADI(2, 0));
  s.mutable_conformation_asym(3).add_actor(ADI(3, 0));

  tie(iter, end) = s.get_interactions<IIC>();
  ASSERT_EQ(iter, end);

  s.mutable_conformation_asym(0).add_actor(ADC(0, '0'));
  s.mutable_conformation_asym(1).add_actor(ADC(1, '0'));
  s.mutable_conformation_asym(2).add_actor(ADC(2, '0'));
  s.mutable_conformation_asym(3).add_actor(ADC(3, '0'));

  // BOOST_FOREACH( Index4 i4, s.get_interactions<IIC>() ) cout << i4 << endl;

  tie(iter, end) = s.get_interactions<IIC>();
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 6));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(0, 7));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 6));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(1, 7));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 6));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(2, 7));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 4));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 5));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 6));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(3, 7));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(4, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(4, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(4, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(4, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(5, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(5, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(5, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(5, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(6, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(6, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(6, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(6, 3));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(7, 0));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(7, 1));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(7, 2));
  ASSERT_NE(iter, end);
  ASSERT_EQ((*iter++).first, Index2(7, 3));
  ASSERT_EQ(iter, end);
}

// TEST(SceneIterator,bodies_hold_actors)
// {
//  using std::make_pair;

//  typedef actor::ActorConcept<X1dim,int> ADI;
//  typedef actor::ActorConcept<X1dim,char> ADC;
//  typedef m::vector< ADI, ADC > Actors;
//  typedef m::vector< ADI, std::pair<ADI,ADC> > Interactions;

//  typedef Scene<Actors,X1dim,size_t> Scene;
//  typedef size_t Index;
//  typedef std::pair<size_t,size_t> Index2;
//  typedef std::pair<Index2,Index2> Index4;

//  typedef std::pair<ADI,ADC> Interaction;
//  typedef SceneIter<Scene,Interaction > Iter2;
//  typedef get_placeholder_type<Interaction,Index>::type PH;

//  shared_ptr<Conformation> conf = make_shared<Conformation>();
//  conf->add_actor( ADI(0.0,0) );
//  // conf->add_actor( ADI(0.0,1) );
//  // conf->add_actor( ADI(0.0,2) );
//  // conf->add_actor( ADC(0.0,'0') );
//  shared_ptr<Conformation> conf2 = make_shared<Conformation>();
//  conf2->add_actor( ADI(1.0,0) );
//  // conf2->add_actor( ADI(1.0,1) );
//  conf2->add_actor( ADC(1.0,'0') );
//  // conf2->add_actor( ADC(1.0,'1') );

//  Scene s;
//  s.add_body(conf);
//  s.add_body(conf2);

//  SceneIter<Scene,ADI> beg,end;
//  tie(beg,end) = s.get_interactions<ADI>();
//  cout << *beg++ << " / " << *end << endl;
//  cout << *beg++ << " / " << *end << endl;
//  cout << *beg++ << " / " << *end << endl;
//  cout << *beg++ << " / " << *end << endl;
//  // BOOST_FOREACH( Index2 ip, iters1 ) cout << ip << endl;

//  // {
//  //  std::pair<Iter2,Iter2> iters2 =
// s.get_interactions<Interaction>();
//  //  BOOST_FOREACH( Index4 ip, iters2 ){
//  //    Interaction::first_type a1;
//  //    Interaction::second_type a2;
//  //    boost::tie(a1,a2) = s.get_interaction<Interaction>(ip);
//  //    cout << ip.first.first << " " << ip.second.first << " "
// <<
// a1 << "    ";
//  //    cout << ip.first.second << " " << ip.second.second << "
// "
// << a2 << endl;
//  //  }
//  // }

//  // {
//  //  SceneIter<Scene,Interaction > beg,end;
//  //  boost::tie(beg,end) = s.get_interactions<Interaction>();
//  //  ASSERT_EQ( *beg, PH(make_pair(make_pair(0,0),make_pair(0,0))) );
// ++beg;
//  // }

//  // ASSERT_EQ( s.body(0).get_actor<ADI>(0), ADI(0.0,0)    );
//  // ASSERT_EQ( s.body(0).get_actor<ADI>(1), ADI(1.0,0)    );
//  // ASSERT_EQ( s.body(0).get_actor<ADC>(0), ADC(2.0,'0')  );
//  // ASSERT_EQ( s.body(1).get_actor<ADI>(0), ADI(0.0,1)    );
//  // ASSERT_EQ( s.body(1).get_actor<ADI>(1), ADI(1.0,1)    );
//  // ASSERT_EQ( s.body(1).get_actor<ADC>(0), ADC(2.0,'1')  );

//  // ASSERT_THROW( s.body(2), std::out_of_range );
//  // ASSERT_THROW( s.body(1).get_actor<ADI>(2), std::out_of_range );
//  // // body.set_position(1.0);
//  // // cout << body.get_actor<ADI>(0) << endl;

//  // cout << s << endl;

// }
}
}
}
