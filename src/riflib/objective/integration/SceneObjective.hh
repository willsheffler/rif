#ifndef INCLUDED_objective_integration_SceneObjective_HH
#define INCLUDED_objective_integration_SceneObjective_HH

#include "riflib/kinematics/SceneBase.hh"

namespace scheme {
namespace objective {
namespace integration {

typedef std::vector<std::pair<int32_t,int32_t> > Rotamers;

// todo: move this into the proper libraries
template< class _Position, class _Index = uint64_t >
struct SceneOjbective {
	typedef _Position Position;
	typedef _Index Index;
	typedef ::scheme::kinematics::SceneBase< Position, Index > SceneBase;
	typedef shared_ptr<SceneBase> SceneP;
	virtual float score( SceneBase const & s ) const = 0;
	virtual void  score( SceneBase const & s , std::vector<float> & vec ) const = 0;
	virtual float score_with_rotamers( SceneBase const & s, Rotamers & rots ) const = 0;
	virtual void  score_with_rotamers( SceneBase const & s , std::vector<float> & vec, Rotamers & rots ) const = 0;
	virtual bool is_compatible( SceneBase const & s ) const = 0;
	virtual bool provides_rotamers() const = 0;
	std::vector<float> scores( SceneBase const & s ) const { std::vector<float> tmp; score(s,tmp); return tmp; }
};

template< class _Scene, class _Objective, class _RotamerMethod = void >
struct SceneObjectiveParametric :
	public SceneOjbective< typename _Scene::Position, typename _Scene::Index >
{
	typedef _Scene Scene;
	typedef _Objective Objective;
	typedef SceneOjbective< typename _Scene::Position, typename _Scene::Index > Super;
	typedef typename Super::SceneBase SceneBase;
	typedef typename Super::SceneP SceneP;
	typedef typename Objective::Config Config;

	Objective objective;
	Config config;

	virtual
	float
	score( SceneBase const & s ) const
	{
		Scene const & scene = static_cast<Scene const &>( s );
		return objective( scene, config ).sum();
	}
	virtual
	float
	score_with_rotamers( SceneBase const & s, Rotamers & rots ) const
	{
		Scene const & scene = static_cast<Scene const &>( s );
		typename Objective::Results results = objective( scene, config );
		fill_rotamers<_RotamerMethod>( results, rots );
		return results.sum();
	}
	virtual
	void
	score( SceneBase const & s, std::vector<float> & vec ) const
	{
		// ::devel::riflib::print_eigenxform( s.position(0) );
		// ::devel::riflib::print_eigenxform( s.position(1) );
		Scene const & scene = static_cast<Scene const &>( s );
		objective( scene, config ).vector( vec );
	}

	virtual
	void
	score_with_rotamers(
		SceneBase const & s ,
		std::vector<float> & vec,
		Rotamers & rots
	) const {
		assert( provides_rotamers() );
		Scene const & scene = static_cast<Scene const &>( s );
		typename Objective::Results results = objective( scene, config );
		results.vector(vec);
		fill_rotamers<_RotamerMethod>( results, rots );
	}

	template< class RotMeth >
	typename boost::disable_if< boost::is_same<void,RotMeth>, void >::type
	fill_rotamers(
		typename Objective::Results const & results,
		std::vector<std::pair<int32_t,int32_t> > & rotamers
	) const {
		rotamers = results.template get<_RotamerMethod>().rotamers_;
	}

	template< class RotMeth >
	typename boost::enable_if< boost::is_same<void,RotMeth>, void >::type
	fill_rotamers(
		typename Objective::Results const & results,
		std::vector<std::pair<int32_t,int32_t> > & rotamers
	) const {
		std::cerr << "rotamer output not supported if _RotamerMethod is void/unspecified" << std::endl;
		std::exit(-1);
		rotamers.clear();
	}


	virtual
	bool
	is_compatible( SceneBase const & s ) const {
		return nullptr != dynamic_cast<Scene const *>( &s );
	}

	virtual
	bool
	provides_rotamers() const {
		return boost::is_same<void,_RotamerMethod>::value;
	}
};


}
}
}

#endif
