#ifndef INCLUDED_objective_voxel_FieldCache_HH
#define INCLUDED_objective_voxel_FieldCache_HH

#include "util/assert.hpp"
#include "objective/voxel/VoxelArray.hpp"
#include "io/cache.hpp"
// #include <boost/exception/all.hpp>
#include <exception>

namespace scheme { namespace objective { namespace voxel {

struct FieldException : public std::runtime_error {
	FieldException(std::string m) : std::runtime_error(m.c_str()) {}
};

template<class Float=float>
struct Field3D {
	virtual Float operator()(Float f, Float g, Float h) const = 0;
};

template<class Float=float>
struct FieldCache3D : public VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef typename BASE::Bounds Float3;
	typedef typename VoxelArray<3,Float>::Indices Indices;

	std::string cache_loc_;

	FieldCache3D() : cache_loc_("") {}

	template<class F1,class F2, class F3>
	FieldCache3D(
		Field3D<Float> const & field,
		F1 const & lb,
		F2 const & ub,
		F3 const & cs,
		std::string const & cache_loc="",
		bool no_init = false,
		int oversample=1
	) : BASE(lb,ub,cs), cache_loc_(cache_loc)
	{
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		#ifdef CEREAL
			if( io::read_cache(cache_loc_,*this) ){
				bool extents_eq = true;
				for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
				if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ){
					std::cout << "EXTENTS EQ, USING CACHE" << std::endl;
					check_against_field( field );
					return;
				}
				std::cout << "Warning, FieldCache3D bad cache " << cache_loc_ << " bounds mismatch, recomputing..." << std::endl;
				this->resize(extents);
			} // else {
		#endif
		// 	std::cout << "NO CACHE" << std::endl;
		// }
		if( !no_init ){
			for(int k = 0; k < this->shape()[2]; ++k){
			for(int j = 0; j < this->shape()[1]; ++j){
			for(int i = 0; i < this->shape()[0]; ++i){
				Float3 cen = this->indices_to_center( Indices(i,j,k) );
				this->operator[]( cen ) = this->sample_field( field, cen, oversample );
				// std::exit(-1);
			}}}
		}
		#ifdef CEREAL
			io::write_cache(cache_loc_,*this);
		#endif

	}

	Float
	sample_field(
		Field3D<Float> const & field,
		Float3 cen,
		int oversample
	) const {
		Float const om = 1.0/oversample;
		Float const oo = om/2.0 - 0.5;
		Float mn = std::numeric_limits<Float>::max();
		for(int o = 0; o < oversample; ++o){
		for(int p = 0; p < oversample; ++p){
		for(int q = 0; q < oversample; ++q){
			Float ox = cen[0] + ( o*om + oo ) * this->cs_[0];
			Float oy = cen[1] + ( p*om + oo ) * this->cs_[1];
			Float oz = cen[2] + ( q*om + oo ) * this->cs_[2];
			// std::cout << ox << " " << oy << " " << oz << " " << cen[0] << " " << cen[1] << " " << cen[2] << std::endl;
			Float fval = field( ox, oy, oz );
			mn = std::min( mn, fval );
		}}}
		return mn;
	}

	double check_against_field( Field3D<Float> const & field, int oversample=1, float tolerance=0.0001 ) const {
		int nerror(0), nsamp(0);
		{
			// sanity check
			// std::exit(-1);
			// std::cout << this->shape()[0] << std::endl;
			// std::cout << this->shape()[1] << std::endl;
			// std::cout << this->shape()[2] << std::endl;
			for(int i = 0; i < this->shape()[0]; i += 8){
			for(int j = 0; j < this->shape()[1]; j += 8){
			for(int k = 0; k < this->shape()[2]; k += 8){
				++nsamp;
				Float3 cen = this->indices_to_center( Indices(i,j,k) );
				Float test1 = this->operator[]( cen );
				// Float test2 = field( cen[0], cen[1], cen[2] );
				Float test2 = sample_field(field,cen,oversample);
				if( fabs(test1-test2) > tolerance ){
					std::cout << "FIELD MISMATCH stored: " << test1 << " recalculated: " << test2 << std::endl;
					++nerror;
					// if(!permissive) throw FieldException("field check fails");
				} else {

				}
			}}}
		}
		// if( nerror==0 ){
		// 	#ifdef USE_OPENMP
		// 	#pragma omp critical
		// 	#endif
		// 	std::cout << "    check_against_field pass " << nsamp-nerror << " of " << nsamp << std::endl;
		// }
		return (double)nerror/(double)nsamp;
	}

};

template<class Float> struct AggMin {
	static Float initval() { return std::numeric_limits<Float>::max(); }
	static Float aggregate( Float agg, Float newval ) { return std::min(agg,newval); }
};
template<class Float> struct AggMax {
	static Float initval() { return std::numeric_limits<Float>::min(); }
	static Float aggregate( Float agg, Float newval ) { return std::max(agg,newval); }
};

template<class Float=float,template<class> class AGG = AggMin>
struct BoundingFieldCache3D : public VoxelArray<3,Float> {
	typedef VoxelArray<3,Float> BASE;
	typedef AGG<Float> Aggregator;
	typedef typename BASE::Bounds Float3;
	// Float spread_;
	// std::string cache_loc_;
	template<class F>
	BoundingFieldCache3D(
		VoxelArray<3,Float> const & ref,
		Float spread,
		F const & cs,
		std::string cache_loc="",
		bool no_init = false
	) : BASE(ref.lb_-spread,ref.ub_+spread,cs) //,
	    // ref_(ref),
	    // spread_(spread),
	    // cache_loc_(cache_loc)
	{
		if( ref.cs_.norm() > spread )
			std::cout << "BoundingFieldCache3D warning: spread " << spread
		              << " less than ref cell_size.norm() " << ref.cs_.norm() << std::endl;
//		Float3 lb = ref.lb_-spread;
//		Float3 ub = ref.ub_+spread;
		typename BASE::Indices extents;
		for(size_t i = 0; i < BASE::DIM; ++i) extents[i] = this->shape()[i];
		#ifdef CEREAL
			if( io::read_cache(cache_loc,*this) ){
				bool extents_eq = true;
				for(size_t i = 0; i < BASE::DIM; ++i) extents_eq &= extents[i] == this->shape()[i];
				if( Float3(lb)==this->lb_ && Float3(ub)==this->ub_ && Float3(cs)==this->cs_ && extents_eq ) return;
				std::cout << "Warning, BoundingFieldCache3D bad cache " << cache_loc << " bounds mismatch, recomputing..." << std::endl;
				this->resize(extents);
			}
		#endif
		if( !no_init ){
			// size_t ncalls = 0;
			for(Float h = this->lb_[2]+this->cs_[2]/2.0; h < this->ub_[2]+this->cs_[2]/2.0; h += this->cs_[2]){
			for(Float g = this->lb_[1]+this->cs_[1]/2.0; g < this->ub_[1]+this->cs_[1]/2.0; g += this->cs_[1]){
			for(Float f = this->lb_[0]+this->cs_[0]/2.0; f < this->ub_[0]+this->cs_[0]/2.0; f += this->cs_[0]){
				// ++ncalls;
				this->operator[]( Float3(f,g,h) ) = calc_agg_val( ref, spread, Float3(f,g,h) );
			}}}
		}
		#ifdef CEREAL
			io::write_cache(cache_loc,*this);
		#endif
	}
	Float calc_agg_val(
		VoxelArray<3,Float> const & ref,
		float spread,
		Float3 const & f3
	) const {
		Float3 beg = ref.lb_.max(f3-spread) + ref.cs_/2.0;
		Float3 end = ref.ub_.min(f3+spread);
		Float3 idx;
		Float val = Aggregator::initval();
		int count = 0;
		for(idx[0] = beg[0]; idx[0] < end[0]; idx[0] += ref.cs_[0]){
		for(idx[1] = beg[1]; idx[1] < end[1]; idx[1] += ref.cs_[1]){
		for(idx[2] = beg[2]; idx[2] < end[2]; idx[2] += ref.cs_[2]){
			if( (f3-idx).squaredNorm() <= spread*spread ){
				val = Aggregator::aggregate( val, ref[idx] );
				++count;
			}
		}}}
		// std::cout << f3 << " | " << beg << " => " << end << std::endl;
		return count==0 ? 0.0 : val;
	}
};


}}}

#endif
