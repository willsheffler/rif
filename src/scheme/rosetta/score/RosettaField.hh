#ifndef INCLUDED_rosetta_objective_RosettaField_HH
#define INCLUDED_rosetta_objective_RosettaField_HH

#include "scheme/objective/voxel/FieldCache.hh"
#include "scheme/rosetta/score/AnalyticEvaluation.hh"
#include "scheme/rosetta/score/EtableParams.hh"
#include "scheme/numeric/util.hh"
#include "scheme/types.hh"
#include <vector>

namespace scheme { namespace rosetta { namespace score {



template< class Atom, class EtableInit >
struct RosettaField {
	typedef util::SimpleArray<3,int> I3;
	typedef util::SimpleArray<3,float> F3;

	EtableParams<float> params;
	std::vector<Atom> atoms_;
	boost::multi_array< std::vector<Atom> , 3 > atom_bins_;
	F3 atom_bins_lb_, atom_bins_ub_;
	I3 atom_bins_dim_;

	float const bin_witdh_ = 6.001f;

	RosettaField() { EtableInit::init_EtableParams(params); }

	RosettaField(
		std::vector<Atom> const & atm
	) : atoms_(atm) {
		EtableInit::init_EtableParams(params); init_atom_bins();
	}

	void init_atom_bins(){
		atom_bins_lb_.fill( std::numeric_limits<float>::max() );
		atom_bins_ub_.fill( std::numeric_limits<float>::min() );
		for( auto const & a : atoms_ ){
			atom_bins_lb_ = atom_bins_lb_.min(a.position());
			atom_bins_ub_ = atom_bins_ub_.max(a.position());
		}
		// std::cout << atom_bins_lb_ << std::endl;
		// std::cout << atom_bins_ub_ << std::endl;
		atom_bins_dim_[0] = std::max( 1, (int)std::ceil ( (atom_bins_ub_[0] - atom_bins_lb_[0]) / bin_witdh_ ) );
		atom_bins_dim_[1] = std::max( 1, (int)std::ceil ( (atom_bins_ub_[1] - atom_bins_lb_[1]) / bin_witdh_ ) );
		atom_bins_dim_[2] = std::max( 1, (int)std::ceil ( (atom_bins_ub_[2] - atom_bins_lb_[2]) / bin_witdh_ ) );
		atom_bins_.resize( atom_bins_dim_ );
		for( auto const & a : atoms_ ){
			I3 i = position_to_atombin(a.position());
			// std::cout << "add atom " << i << std::endl;
			atom_bins_(i).push_back(a);
		}

		int tot = 0;
		for( int i = 0; i < atom_bins_.num_elements(); ++i){
			tot += atom_bins_.data()[i].size();
		}
		BOOST_VERIFY( tot == atoms_.size() );
	}
	I3 position_to_atombin( F3 p ) const {
		I3 i = ( p - atom_bins_lb_ ) / bin_witdh_;
		// std::cout << p << std::endl;
		// std::cout << i << std::endl;
		return i;
	}


	float
	compute_rosetta_energy_one(
		Atom const & a,
		float x, float y, float z,
		int atype
	) const {
		float const dx = x-a.position()[0];
		float const dy = y-a.position()[1];
		float const dz = z-a.position()[2];
		float const dis2 = dx*dx+dy*dy+dz*dz;
		if( dis2 > 36.0 ) return 0; //  103s vs 53s
		float const dis = std::sqrt(dis2);
		float const inv_dis2 = 1.0f/dis2;
		float atr0=0,rep0=0,sol0=0;
		int at = a.type();
		bool very_repulsive = at==-12345;
		if( very_repulsive ) at = 5;
		bool neg_only = at < 0;
		at = abs(at);
		EtableParamsOnePair<float> const & p = params.params_for_pair( at, atype );
		lj_evaluation( p, dis, dis2, inv_dis2, atr0, rep0);
		lk_evaluation( p, dis, inv_dis2, sol0 );
		atr0 = neg_only ? 0.0 : atr0;
		if( very_repulsive ){
			// so hacky!!! in rif apo hsearch, this is the minimal setting that'll prevent
			// anything occulding the hydrophobic fake-biotin tail
			rep0 = std::max( rep0, ::scheme::numeric::sigmoidish( dis2, 3.0f, 4.6f ) * 1.0f );
		}
		return 0.8*atr0 + 0.44*rep0 + 0.75*sol0;
	}

	float compute_rosetta_energy_safe(float x, float y, float z, int atype) const
	{
		float E = 0;
		for( auto const & a : atoms_ ){
			E += compute_rosetta_energy_one( a, x, y, z, atype );
		}
		return E;
	}

	float compute_rosetta_energy(float x, float y, float z, int atype) const
	{
		I3 i = position_to_atombin( F3(x,y,z) );
		I3 lb = I3(  0, 0, 0 ).max( i-1 );
		I3 ub = atom_bins_dim_.min( i+2 );
		// std::cout << lb << " " << i << " " << ub << std::endl;
		BOOST_VERIFY( lb[0] <= ub[0] );
		BOOST_VERIFY( lb[1] <= ub[1] );
		BOOST_VERIFY( lb[2] <= ub[2] );
		float E = 0;
		I3 ii;
		for( ii[0] = lb[0]; ii[0] < ub[0]; ++ii[0] ){
		for( ii[1] = lb[1]; ii[1] < ub[1]; ++ii[1] ){
		for( ii[2] = lb[2]; ii[2] < ub[2]; ++ii[2] ){
			// std::cout << i << " " << ii << " " << std::endl;
			for( auto const & a : atom_bins_(ii) ){
				E += compute_rosetta_energy_one( a, x, y, z, atype );
			}
		}}}
		return E;
	}

	template<class F>
	float compute_rosetta_energy( F const & f, int atype ) const
	{
		return compute_rosetta_energy( f[0], f[1], f[2], atype );
	}
	template<class F>
	float compute_rosetta_energy_safe( F const & f, int atype ) const
	{
		return compute_rosetta_energy_safe( f[0], f[1], f[2], atype );
	}

};

template< class Atom, class EtableInit >
struct RosettaFieldAtype : objective::voxel::Field3D<float> {
	RosettaField<Atom,EtableInit> const & rf_;
	int atype_;
	RosettaFieldAtype(RosettaField<Atom,EtableInit> const & rf, int atype) : rf_(rf),atype_(atype) {}
	float operator()(float x, float y, float z) const {
		return rf_.compute_rosetta_energy(x,y,z,atype_);
	}
};


}}}

#endif
