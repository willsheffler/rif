#include <gtest/gtest.h>

#include "kinematics/Scene_io.hpp"
#include "actor/ActorConcept_io.hpp"
#include "objective/ObjectiveFunction.hpp"
#include "objective/ObjectiveVisitor.hpp"
#include "numeric/X1dim.hpp"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/io.hpp>

#include <stdint.h>

namespace scheme { namespace kinematics { namespace objective_test {

using std::cout;
using std::endl;
using boost::tie;
using numeric::X1dim;

typedef actor::ActorConcept<X1dim,int> ADI;
typedef actor::ActorConcept<X1dim,char> ADC;	
typedef size_t Index;
typedef std::pair<size_t,size_t> Index2;
typedef std::pair<Index2,Index2> Index4;

////////////// test scores ////////////////////////

	struct ScoreADI {
		typedef double Result;
		typedef ADI Interaction;
		static std::string name(){ return "ScoreADI"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			return a.data_;
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADI const& s){ return out << s.name(); }

	struct ScoreADIADI {
		static size_t ncalls;
		typedef double Result;
		typedef ADI Actor1;
		typedef ADI Actor2;
		typedef std::pair<ADI,ADI> Interaction;
		static std::string name(){ return "ScoreADIADI"; }
		template<class Config>
		Result operator()(Actor1 const & a1, Actor2 const & a2, Config const& ) const {
			++ncalls;
			// cout << a.first << " " << a.second << endl;
			return distance(a1.position(),a2.position());
		}
		template<class Config>
		Result operator()(Interaction const & i, Config const& c) const {
			return this->template operator()<Config>(i.first,i.second,c);
		}
	};
	size_t ScoreADIADI::ncalls = 0;
	std::ostream & operator<<(std::ostream & out,ScoreADIADI const& s){ return out << s.name(); }

	struct ScoreADC {
		typedef double Result;
		typedef ADC Interaction;
		static std::string name(){ return "ScoreADC"; }
		template<class Config>
		Result operator()(Interaction const & a, Config const& ) const {
			switch(a.data_){
				case '1': return 1; break;
				case '2': return 2; break;
				case '3': return 3; break;
				case '4': return 4; break;
			};
			return 0;
		}
	};
	std::ostream & operator<<(std::ostream & out,ScoreADC const& s){ return out << s.name(); }

	struct ScoreADCADI {
		typedef double Result;
		typedef std::pair<ADC,ADI> Interaction;
		static std::string name(){ return "ScoreADCADI"; }
		template<class Config>
		Result operator()(ADC const & adc, ADI const & adi, Config const& ) const {
			return 10*distance(adc.position(),adi.position());
		}
		template<class Config>
		Result operator()(Interaction const & i, Config const& c ) const {
			return this->operator()( i.first, i.second, c);
		}

	};
	std::ostream & operator<<(std::ostream & out,ScoreADCADI const& s){ return out << s.name(); }

	struct Config {};

/////////////////// tests //////////////////////////

TEST(SceneObjective,basic){
	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreADI,
			ScoreADC,
			ScoreADIADI,
			ScoreADCADI
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;
	ObjFun score;

	typedef m::vector< ADI, ADC > Actors;
	typedef Scene<Actors,X1dim,size_t> Scene;

	Scene scene(3);
	ASSERT_EQ( score(scene).sum(), 0 );

	scene.mutable_conformation_asym(0).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,2) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,3) );
	scene.mutable_conformation_asym(0).add_actor( ADI(0,4) );
	ASSERT_EQ( score(scene), Results(10,0,0,0) );

	scene.mutable_conformation_asym(1).add_actor( ADI(1,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,3) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,4) );
	ASSERT_EQ( score(scene), Results(20,0,16,0) );

	scene.mutable_conformation_asym(0).add_actor( ADC(0,'1') );
	scene.mutable_conformation_asym(0).add_actor( ADC(0,'2') );
	ASSERT_EQ( score(scene), Results(20,3,16,80) );

	scene.mutable_conformation_asym(2).add_actor( ADI(2,100) );
	ASSERT_EQ( score(scene), Results(120,3,28,120) );

	scene.mutable_conformation_asym(2).add_actor( ADI(2,1000) );
	ASSERT_EQ( score(scene), Results(1120,3,40,160) );

}

TEST(SceneObjective,symmetry){
	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreADI,
			ScoreADC,
			ScoreADIADI,
			ScoreADCADI
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;

	ObjFun score;

	typedef m::vector< ADI, ADC > Actors;
	typedef Scene<Actors,X1dim,size_t> Scene;

	Scene scene(2);
	scene.add_symframe(10);

	ASSERT_EQ( score(scene).sum(), 0 );

	scene.mutable_conformation_asym(0).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	// check that symmetric interactiosn are downweighted by 0.5: 40/2 = 20
	ASSERT_EQ( score(scene), Results(3,0,1+20,0) ); 

	scene.mutable_conformation_asym(1).add_actor( ADI(0,1) );
	scene.mutable_conformation_asym(1).add_actor( ADI(1,2) );
	// check that symmetric interactiosn are downweighted by 0.5: 160/2 = 80
	ASSERT_EQ( score(scene), Results(6,0,2+80,0) );

	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'1') );
	// scene.mutable_conformation_asym(1).add_actor( ADC(0,'2') );
	// ASSERT_EQ( score(scene), Results(20,3,16,80) );

}





template<class Scene,class Visitor>
void performance_test_helper(Scene const & scene, Visitor & visitor){
	typedef ADI Actor1;
	typedef ADI Actor2;
	typedef typename Scene::Conformation Conformation;
	typedef typename Scene::Position Position;
	typedef typename Scene::Index Index;

		Index const NBOD = (Index)scene.bodies_.size();
		Index const NSYM = (Index)scene.symframes_.size()+1;
		for(Index i1 = 0; i1 < NBOD*NSYM; ++i1){
			Conformation const & c1 = scene.conformation(i1);
			Position     const & p1 =     scene.position(i1);
			Index const NACT1 = (Index)c1.template get<Actor1>().size();
			for(Index i2 = 0; i2 < NBOD*NSYM; ++i2){
				if( i1 >= NBOD && i2 >= NBOD ) continue;
				if( i2 <= i1 ) continue;
				Conformation const & c2 = scene.conformation(i2);
				Position     const & p2 =     scene.position(i2);
				Index const NACT2 = (Index)c2.template get<Actor2>().size();
				for(Index j1 = 0; j1 < NACT1; ++j1){
					Actor1 a1( c1.template get<Actor1>()[j1], p1 );
					for(Index j2 = 0; j2 < NACT2; ++j2){
						Actor2 a2( c2.template get<Actor2>()[j2], p2 );
						visitor( std::make_pair(a1,a2), i1<NBOD&&i2<NBOD?1.0:0.5 );
					}
				}
			}
		}
	cout << visitor.result_ << " " << (double)ScoreADIADI::ncalls/1000000.0 << endl;
	ScoreADIADI::ncalls = 0;
}



#ifdef SCHEME_BENCHMARK

TEST( SceneObjective, big_scene_performance ){
// TEST(SceneObjective,performance){
	// TODO: speed up SceneIter iteration 
	//       iteration seems to take about 100 cycles per score call overhead
	//       much of this is probably all the conditions for symmetry checks
	//       could template out these and have both sym and asym scenes?
	// FIX with visitation pattern, seems at least 10x faster

	std::cout << "This test performs 301.934M score calls, should take about a second when compiled with optimizations." << std::endl;

	typedef	objective::ObjectiveFunction<
		m::vector<
			ScoreADI,
			ScoreADC,
			ScoreADIADI,
			ScoreADCADI
		>,
		Config
	> ObjFun;
	typedef ObjFun::Results Results;

	ObjFun score;

	typedef m::vector< ADI, ADC > Actors;
	typedef Scene<Actors,X1dim,uint32_t> Scene;

	ScoreADIADI obj;
	Config c;
	objective::ObjectiveVisitor<ScoreADIADI,Config> visitor(obj,c);

	Scene scene; {
		Scene::Index const NBOD = 10;
		Scene::Index const NSYM = 20;
		Scene::Index const NACT = 400;
		for(Scene::Index i = 0; i < NBOD; ++i) scene.add_body();		
		for(Scene::Index i = 0; i < NSYM-1; ++i) scene.add_symframe(i+1);
		for(Scene::Index i = 0; i < NBOD; ++i){
			for(Scene::Index j = 0; j < NACT; ++j){
				scene.mutable_conformation_asym(i).add_actor( ADI(j,i) );
			}
		}
	}

	cout << score(scene).get<ScoreADIADI>() << " " << (double)ScoreADIADI::ncalls/1000000.0 << "M" << endl;
	ScoreADIADI::ncalls = 0;
	return;


	if(false)
	{
			// typedef ADI Actor1;
			// typedef ADI Actor2;
			// typedef Scene::Position Position;
			// Scene::Index const NBOD = scene.bodies_.size();
			// Scene::Index const NSYM = scene.symframes_.size()+1;
			// for(Scene::Index i1 = 0; i1 < NBOD*NSYM; ++i1){
			// 	Conformation const & c1 = scene.conformation(i1);
			// 	Position     const & p1 =     scene.position(i1);
			// 	Scene::Index const NACT1 = c1.get<Actor1>().size();
			// 	for(Scene::Index i2 = 0; i2 < NBOD*NSYM; ++i2){
			// 		if( i1 >= NBOD && i2 >= NBOD ) continue;
			// 		if( i2 <= i1 ) continue;
			// 		Conformation const & c2 = scene.conformation(i2);
			// 		Position     const & p2 =     scene.position(i2);
			// 		Scene::Index const NACT2 = c2.get<Actor2>().size();					
			// 		for(Scene::Index j1 = 0; j1 < NACT1; ++j1){
			// 			Actor1 a1( c1.get<Actor1>()[j1], p1 );
			// 			for(Scene::Index j2 = 0; j2 < NACT2; ++j2){
			// 				Actor2 a2( c2.get<Actor2>()[j2], p2 );
			// 				visitor( a1, a2, i1<NBOD&&i2<NBOD?1.0:0.5 );
			// 			}
			// 		}
			// 	}
			// }
			// cout << visitor.result_ << " " << (double)ScoreADIADI::ncalls/1000000.0 << endl;
			// ScoreADIADI::ncalls = 0;

	}

	if(false)
		performance_test_helper(scene,visitor);
	ScoreADIADI::ncalls = 0;


	scene.visit(visitor);
	cout << visitor.result_ << " " << (double)ScoreADIADI::ncalls/1000000.0 << "M" << endl;
	ScoreADIADI::ncalls = 0;

}

#endif 

}
}
}
