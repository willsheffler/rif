
#include "riflib/kinematics/Director.hpp"

#include <list>

namespace scheme { namespace search {

template< class Xform, class Index >
struct BoundingFunction {
	typedef typename Xform::Scalar Float;
	virtual void set_resolution() = 0;
	virtual Float evaluate( kinematics::SceneBase<Xform,Index> const & scene ) const = 0;
};


template< class Scene, class Objective >
struct BoundingObjectiveFunction
 : public BoundingFunction< typename Scene::Position, typename Scene::Index >
{
	typedef std::list< std::pair< float, Objective > > ObjectiveList;
	ObjectiveList objectives_;
	Objective * active_objective_;

	void add_objective( float bounding_radius, Objective const & objective ){
		typename ObjectiveList::iterator i;
		for( i = objectives_.begin(); i != objectives_.end(); ++i ){
			if( bounding_radius > i->first ){
				if( i != objectives_.begin() ) --i;
				break;
			}
		}
		objectives_.insert( i, std::make_pair(bounding_radius,objective) );
	}
};


template< class BigIndex, class Float = float >
struct SpatialBandBResult {
	std::vector< std::pair< float, BigIndex > > results;
};

template<
	class _Xform,
	class _BigIndex = uint64_t,
	class _Index = uint64_t
>
struct SpatialBandB {
	typedef _Xform Xform;
	typedef _BigIndex BigIndex;
	typedef _Index Index;
	typedef kinematics::SceneBase<Xform,Index> Scene;
	typedef scheme::shared_ptr< Scene > SceneP;
	typedef scheme::shared_ptr< BoundingFunction<Xform,Index> > BoundP;
	typedef scheme::shared_ptr< kinematics::Director<Xform,Index,Index> > DirectorP;
	typedef SpatialBandBResult< BigIndex, float > Result;
	typedef scheme::shared_ptr< Result > ResultP;

	BoundP bounding_func_;
	SceneP scene_;
	DirectorP director_;


	SpatialBandB(){}

	ResultP search(){
		ResultP result = ResultP( new Result );

		return result;
	}

};

}}
