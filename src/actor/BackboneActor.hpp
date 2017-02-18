#ifndef INCLUDED_actor_BackboneActor_HH
#define INCLUDED_actor_BackboneActor_HH

#include <Eigen/Dense>
#include "io/dump_pdb_atom.hpp"

namespace scheme {
namespace actor {




	template<class _Position>
	struct BackboneActor {

		///@brief Position type, leave out to make actor "Fixed"
		typedef _Position Position;
		typedef typename Position::Scalar Float;
		typedef BackboneActor<Position> THIS;
		typedef Eigen::Matrix<Float,3,1> V3;

		Position position_;
		char aa_,ss_;
		int index_;

		BackboneActor() : position_(Position::Identity()), aa_('-'), ss_('-'), index_(0) {}

		BackboneActor(Position const & p, char aa, char ss, int i=0) : position_(p), aa_(aa), ss_(ss),index_(i) {}

		template <class V>
		BackboneActor(
			V n, 
			V ca, 
			V c, 
			char aa='-', 
			char ss='-', 
			int i=0
		) : aa_(aa), ss_(ss), index_(i) {
			from_n_ca_c(n,ca,c);
		}

		BackboneActor(
			BackboneActor const & actor0,
			Position const & moveby
		){
			aa_ = actor0.aa_;
			ss_ = actor0.ss_;
			index_ = actor0.index_;			
			set_position(moveby*actor0.position());
		}

		template <class Vin>
		void from_n_ca_c( Vin _n, Vin _ca, Vin _c ){
			
			V3 n (  _n[0],  _n[1],  _n[2] );
			V3 ca( _ca[0], _ca[1], _ca[2] );
			V3 c (  _c[0],  _c[1],  _c[2] );						
			// V3 e1 = n - ca; 
			V3 e1 = (c+n)/2.0 - ca; // from old motif stuff to maintain compatibility			
			e1.normalize();
			V3 e3 = e1.cross( c - ca );
			e3.normalize();
			V3 e2 = e3.cross( e1 );
			Eigen::Matrix<Float,3,3> m;
			m(0,0) = e1[0];   m(0,1) = e2[0];   m(0,2) = e3[0];
			m(1,0) = e1[1];   m(1,1) = e2[1];   m(1,2) = e3[1];
			m(2,0) = e1[2];   m(2,1) = e2[2];   m(2,2) = e3[2];			

			V3 t = m * V3( -1.952799123558066, -0.2200069625712990, 1.524857 ) + ca;
			// V t = m * V(-0.865810,-1.764143,1.524857) + ca;// - (c+n)/2.0; // average 'CEN' icoor

			// std::cout << e1 << std::endl;
			// std::cout << e2 << std::endl;
			// std::cout << e3 << std::endl;						

			position_.setIdentity();
			position_.rotate( m );
			position_.translation()[0] = t[0];
			position_.translation()[1] = t[1];
			position_.translation()[2] = t[2];						
	
	// from_four_points(c,u,v,w)
	// 	Vector e1( u - v);
	// e1.normalize();
	// Vector e3( cross( e1, w - v ) );
	// e3.normalize();
	// Vector e2( cross( e3,e1) );
	// R.col_x( e1 ).col_y( e2 ).col_z( e3 );
	// t = c;	
		}

		template<class V>
		void get_n_ca_c( V & n, V & ca, V & c ) const {
			V3 tmpn  = position_ * V( 2.80144, -0.992889, -1.52486 );
		 	V3 tmpca = position_ * V( 1.95280,  0.220007, -1.52486 );
		 	V3 tmpc  = position_ * V( 2.87767,  1.4329  , -1.52486 );
		 	n [0] = tmpn [0]; n [1] = tmpn [1]; n [2] = tmpn [2];
		 	ca[0] = tmpca[0]; ca[1] = tmpca[1]; ca[2] = tmpca[2];
		 	c [0] = tmpc [0]; c [1] = tmpc [1]; c [2] = tmpc [2];
		}

		template<class V>
		void get_ca( V & ca ) const {
		 	V3 tmp = position_ * V3( 1.95280,  0.220007, -1.52486 );
		 	ca[0] = tmp[0]; ca[1] = tmp[1]; ca[2] = tmp[2];		 			 	
		}

		void 
		set_position(
			Position const & pos
		){ position_ = pos; }

		void 
		moveby(
			Position const & pos
		){ position_ = pos * position_; }

		Position const &
		position() const { return position_; }

		bool operator==(THIS const & o) const {
			return o.position_==position_ && 
			       o.aa_      ==aa_       &&
			       o.ss_      ==ss_       &&
			       o.index_   ==index_;
		}

		///@brief necessary for testing only
		// bool operator<(THIS const & o) const { return std::make_pair(position_,aa_) < std::make_pair(o.position_,o.aa_); }

		///@brief necessary for testing only
	  	template<class Archive> void serialize(Archive & ar, const unsigned int ){
	  		ar & position_;
	  		ar & aa_;
	  		ar & ss_;	  		
	  		ar & index_;
	  	}

	};


template< class P, class MetaData >
void write_pdb( std::ostream & out, BackboneActor<P> const & a, MetaData const & meta ){
	typedef Eigen::Matrix<typename P::Scalar,3,1> V3;
	V3 n,ca,c;
	a.get_n_ca_c(n,ca,c);
	int rnum = 1;
	int anum = 1;
	chemical::AtomData  ndata( " N  ", "GLY", 'A', rnum, anum+1, "N" );
	chemical::AtomData cadata( " CA ", "GLY", 'A', rnum, anum+2, "C" );
	chemical::AtomData  cdata( " C  ", "GLY", 'A', rnum, anum+3, "C" );
	io::dump_pdb_atom( out,  n,  ndata );
	io::dump_pdb_atom( out, ca, cadata );
	io::dump_pdb_atom( out,  c,  cdata );

	// AtomData(
	// 	std::string const & _atomname = AtomData::default_atomname(),        
	// 	std::string const & _resname  = AtomData::default_resname(),       
	// 	char                _chain    = AtomData::default_chain(),     
	// 	int                 _resnum   = AtomData::default_resnum(),      
	// 	int                 _atomnum  = AtomData::default_atomnum(),       
	// 	std::string const & _elem     = AtomData::default_elem(),    
	// 	bool                _ishet    = AtomData::default_ishet(),     
	// 	float               _occ      = AtomData::default_occ(),   
	// 	float               _bfac     = AtomData::default_bfac()
	// );

}



template<class X>
std::ostream & operator<<(std::ostream & out,BackboneActor<X> const& a){
	return out << "BackboneActor " << a.ss_ << " " << a.aa_ << " " << a.index_;
}

}
}

#endif
