#ifndef INCLUDED_objective_hash_XformHashNeighbors_HH
#define INCLUDED_objective_hash_XformHashNeighbors_HH

#include "riflib/util/SimpleArray.hh"
#include "riflib/numeric/rand_xform.hh"
#include "riflib/util/dilated_int.hh"

#include <fstream>
#include "riflib/io/dump_pdb_atom.hh"

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <random>
#include <set>
#include <map>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/foreach.hpp>


namespace scheme { namespace objective { namespace hash {


template< class XformHash, bool UNIQUE=false > struct XformHashNeighbors;


/// may contain duplicate keys!
template< class XformHash, bool UNIQUE >
struct XformHashNeighborCrappyIterator : boost::iterator_facade<
	XformHashNeighborCrappyIterator<XformHash,UNIQUE>,
	typename XformHash::Key,
	boost::forward_traversal_tag,
	typename XformHash::Key
>
{
	typedef XformHashNeighborCrappyIterator<XformHash,UNIQUE> THIS;
	typedef typename XformHash::Key Key;
	int i1, i2, i3, ix, iy, iz;
	bool end;
	XformHash const * xh;
	std::vector<uint64_t> const * ori_nbrs_;
	std::vector< util::SimpleArray<3,int16_t> > const * shifts_;
	// std::set<Key> seenit_;
	// google::dense_hash_set<Key> seenit_;

	XformHashNeighborCrappyIterator( XformHashNeighbors<XformHash,UNIQUE> & xhn, Key key, bool _end=false) : 
		xh( &xhn.hasher_ ),
		ori_nbrs_( &xhn.get_ori_neighbors( key ) ),
		shifts_( &xhn.get_cart_shifts() )
	{
		// if( UNIQUE ) seenit_.set_empty_key( std::numeric_limits<Key>::max() );
		i1 = i2 = i3 = 0;
		end = false;
		ix = (int)((util::undilate<7>( key>>1 ) & 63) | ((key>>57)&127)<<6);
		iy = (int)((util::undilate<7>( key>>2 ) & 63) | ((key>>50)&127)<<6);
		iz = (int)((util::undilate<7>( key>>3 ) & 63) | ((key>>43)&127)<<6);
		end = _end;
	}
private:
    friend class boost::iterator_core_access;
    void increment(){
    	++i3;
    	if( i3 == 2                 ){ i3 = 0; ++i2; }
    	if( i2 == shifts_->size()   ){ i2 = 0; ++i1; }
    	if( i1 == ori_nbrs_->size() ){ end = true; }
    }
    Key dereference() const {
    	if( end ) return std::numeric_limits<Key>::max();
		Key ori_key = (*ori_nbrs_)[i1];
		ori_key = xh->cart_shift_key( ori_key, ix, iy, iz );
		Key k = xh->cart_shift_key( ori_key, (*shifts_)[i2][0], (*shifts_)[i2][1], (*shifts_)[i2][2], i3 );
		// if( UNIQUE ){
		// 	if( seenit_.find(k) != seenit_.end() ){
		// 		const_cast<THIS*>(this)->increment();
		// 		return dereference();
		// 	} else {
		// 		const_cast<THIS*>(this)->seenit_.insert(k);
		// 	}		
		// }
		return k;
    }
    bool equal( THIS const & o) const {
    	if( end && o.end ) return true;
    	else return false;
    	std::cerr << "will needs to implement XformHashNeighborCrappyIterator comparison" << std::endl;
    }
};
template< class XformHash, bool UNIQUE >
std::ostream & operator<<( std::ostream & out, XformHashNeighborCrappyIterator<XformHash,UNIQUE> const & i ){
	return out << "XformHashNeighborCrappyIterator( " << i.i1 << " " << i.i2 << " " << i.i3 << " " << i.end << " )";
}


template< class H > 	
struct InitHash {
	static void init_hash( H & ){}
};

template< class K, class V >
struct InitHash< google::dense_hash_map<K,V> > {
	static void init_hash( google::dense_hash_map<K,V> & h ){
		h.set_empty_key( std::numeric_limits<K>::max() );
	}	
};

template< class XformHash, bool UNIQUE >
struct XformHashNeighbors {
	typedef typename XformHash::Key Key;
	typedef typename XformHash::Float Float;
	typedef typename XformHash::Xform Xform;
	typedef XformHashNeighborCrappyIterator<XformHash,UNIQUE> crappy_iterator;

	XformHash const hasher_;
	Float cart_bound_, ang_bound_, quat_bound_;
	int nsamp_;
	// google::sparse_hash_map< Key, std::vector<Key> > ori_cache_;
	typedef std::map< Key, std::vector<Key> > OriCache;
	// typedef google::dense_hash_map< Key, std::vector<Key> > OriCache;
	OriCache ori_cache_;	
	std::vector< util::SimpleArray<3,int16_t> > cart_shifts_;

	size_t n_queries_, n_cache_miss_;


	XformHashNeighbors(
		Float cart_bound,
		Float ang_bound, // degrees
		XformHash xh,
		double xcov_target=100.0
	) : hasher_(xh) {
		InitHash<OriCache>::init_hash(ori_cache_);
		cart_bound_ = cart_bound;
		ang_bound_ = ang_bound;
		quat_bound_ = numeric::deg2quat(ang_bound);
		Float ang_nside = 2.0 * quat_bound_ / hasher_.ang_width();
		nsamp_ = ang_nside*ang_nside*ang_nside*3.0*xcov_target; // ~100 fold coverage
		fill_cart_shifts();
		// std::cout << "XformHashNeighbors aw=" << hasher_.ang_width() << ", quat_bound=" << quat_bound_ 
		          // << " ang_nside=" << ang_nside << " nsamp=" << nsamp_ << " cart_shifts=" << cart_shifts_.size() << std:: endl;
	}

	void fill_cart_shifts(){
		Float thresh = cart_bound_+hasher_.cart_width();
		int16_t di = std::ceil( cart_bound_ / hasher_.cart_width() ) ;
		assert( di >= 0 );
		for(int16_t i = -di; i <= di; ++i){
		for(int16_t j = -di; j <= di; ++j){
		for(int16_t k = -di; k <= di; ++k){					
			if( 
				Eigen::Vector3d( i-0.5, j-0.5, k-0.5 ).norm()*hasher_.cart_width() < thresh || 
				Eigen::Vector3d( i    , j    , k     ).norm()*hasher_.cart_width() < thresh || 
				Eigen::Vector3d( i+0.5, j+0.5, k+0.5 ).norm()*hasher_.cart_width() < thresh ||
				false
			){
				cart_shifts_.push_back( util::SimpleArray<3,int16_t>(i,j,k) );
			}
		}}}
	}

	crappy_iterator neighbors_begin( Key key ) {
		return crappy_iterator(*this,key);
	}
	crappy_iterator neighbors_end( Key key ) {
		return crappy_iterator(*this,key,true);
	}
	std::pair<crappy_iterator,crappy_iterator> neighbors(Key key) {
		return std::make_pair( neighbors_begin(key), neighbors_end(key) );
	}

	std::vector< util::SimpleArray<3,int16_t> > const & get_cart_shifts() const {
		return cart_shifts_;
	}

	std::vector<Key> const & get_ori_neighbors( Key key ){
		++n_queries_;
		Key ori_key = key & XformHash::ORI_MASK;
		if( ori_cache_.find(ori_key) == ori_cache_.end() ){
			++n_cache_miss_;
			Key ksym, asym_key = hasher_.asym_key( ori_key, ksym );
			// if( asym_key==ori_key ){
			// TODO: figure out how to get quat key symmetries working
			if( true ){
				// std::cout << "get_asym nbrs the hard way, store in " << ori_key << std::endl;
				std::mt19937 rng((unsigned int)time(0) + 23058704);
				Xform c = hasher_.get_center(key);
				// c.translation()[0] = c.translation()[1] = c.translation()[2] = 0;
				std::set<Key> keys;
				for(int i = 0; i < nsamp_; ++i){
					Xform p; numeric::rand_xform_quat(rng,p,cart_bound_,quat_bound_);
					p.translation()[0] = p.translation()[1] = p.translation()[2] = 0;
					Key nbkey = hasher_.get_key( p * c ) & XformHash::ORI_MASK;
					// assert( ( nbkey & ~XformHash::ORI_MASK ) == 0 );
					keys.insert(nbkey);
				}
				assert( keys.size() > 0 );
				// #ifdef USE_OPENMP
				// #pragma omp critical
				// #endif
				{
					ori_cache_.insert( std::make_pair( ori_key, std::vector<Key>(keys.size()) ) );
					std::copy(keys.begin(),keys.end(),ori_cache_[ori_key].begin());
				}

				// std::ofstream out("nbrs_asym.pdb");
				// for(int i = 0; i < ori_cache_[ori_key].size(); ++i){
				// 	Xform x = hasher_.get_center( ori_cache_[ori_key][i] );
				// 	Eigen::Quaterniond q(x.rotation());
				// 	io::dump_pdb_atom( out, "C" ,i, 50*q.coeffs() );
				// }
				// out.close();

			} else {
				std::vector<Key> const & asym_nbrs = get_ori_neighbors( asym_key );
				// std::cout << "set sym_key nbrs from asym_nbrs " << asym_nbrs.size() << ", store in " << ori_key << " " << ori_cache_.size() << std::endl;
				ori_cache_.insert( std::make_pair( ori_key, std::vector<Key>( asym_nbrs.size() ) ) );
				for(int i = 0; i < asym_nbrs.size(); ++i){
					Key new_key = hasher_.sym_key( asym_nbrs[i], ksym );
					ori_cache_[ori_key][i] = new_key;
					assert( new_key == (new_key&hasher_.ORI_MASK) );
				}

				std::ofstream out("nbrs_sym.pdb");
				std::cout << "FLIP: " << ((ksym>>2)&1) << " " << ((ksym>>1)&1) << " " << ((ksym>>0)&1) << std::endl;
				for(int i = 0; i < asym_nbrs.size(); ++i){
					Xform x = hasher_.get_center( ori_cache_[ori_key][i] );
					Eigen::Quaternion<Float> q(x.rotation());
					Eigen::Matrix<Float,4,1> v( q.coeffs() );
					// v *= q.coeffs()[3];
					// std::cout << q.w() << std::endl;
					io::dump_pdb_atom( out, "C" ,i, 50.0*v );
					// std::cout << std::showpos;
					// std::cout << Eigen::Quaterniond(hasher_.get_center(ori_cache_[ori_key][i]).rotation()).coeffs().transpose() << std::endl
					//           << Eigen::Quaterniond(hasher_.get_center(          asym_nbrs[i]).rotation()).coeffs().transpose() << std::endl;
				}
				out.close();

			}
			// for(int i = 0; i < ori_cache_[ori_key].size(); ++i){
			// 	std::cout << ori_cache_[ori_key][i] << std::endl;
			// }
			// std::cout << "added ori_cahce_ " << ori_key << " " << (float)nsamp_/keys.size() << std::endl;
		}
		return ori_cache_[ori_key];
	}

	void merge( XformHashNeighbors<XformHash,UNIQUE> const & other ){
		BOOST_FOREACH( typename OriCache::value_type const & v, other.ori_cache_ ){
			if( ori_cache_.find( v.first ) == ori_cache_.end() ){
				ori_cache_.insert( v );
			} else {
				std::vector<Key> & mykeys( ori_cache_.find( v.first )->second );
				BOOST_FOREACH( Key k, v.second ){
					if( std::find( mykeys.begin(), mykeys.end(), k ) == mykeys.end() ){
						mykeys.push_back(k);
					}
				}
			}
		}
	}

	bool save( std::ostream & out ) {
		// no way to check if the stream was opened binary!
		// if( ! (out.flags() & std::ios::binary) ){
		// 	std::cerr << "XformMap::save must be binary ostream" << std::endl;
		// 	return false;
		// }
		std::ostringstream oss;
		oss << std::endl;
		oss << "=========== description ===========" << std::endl;
		oss << "Scheme XformHashNeighbors" << std::endl;
		oss << "Hasher: " << hasher_.name() << std::endl;
		oss << "Cart Bound: " << cart_bound_ << std::endl;
		oss << "Angular Bound: " << ang_bound_ << std::endl;
		oss << "=========== begin binary data ===========" << std::endl;
		size_t s = oss.str().size();
		out.write((char*)&s,sizeof(size_t));
		out.write(oss.str().c_str(),s);
		// begin binary data
		s = hasher_.name().size();
		// std::cout << "SIZE " << s << std::endl;
		out.write( (char*)&s, sizeof(size_t) );
		out.write( hasher_.name().c_str(), hasher_.name().size()*sizeof(char) );

		out.write( (char*)&hasher_, sizeof(XformHash) );
		out.write( (char*)&cart_bound_, sizeof(Float) );		
		out.write( (char*)&ang_bound_, sizeof(Float) );
		out.write( (char*)&quat_bound_, sizeof(Float) );
		out.write( (char*)&nsamp_, sizeof(int) );		

		// std::map< Key, std::vector<Key> > ori_cache_;
		// std::vector< util::SimpleArray<3,int16_t> > cart_shifts_;
		size_t nentries = ori_cache_.size();
		out.write( (char*)&nentries, sizeof(size_t) );
		for(typename OriCache::iterator j = ori_cache_.begin(); j != ori_cache_.end(); ++j){
			Key key = j->first;
			out.write( (char*)&key, sizeof(Key) );
			s = j->second.size();
			out.write( (char*)&s, sizeof(size_t) );
			for(size_t k = 0; k < s; ++k){
				Key key2 = j->second[k];
				out.write( (char*)&key2, sizeof(Key) );
			}
		}

		// shifts computed by c'tor, no need to save them

		return true;
	}
	bool load( std::istream & in ) {
		// no way to check if the stream was opened binary!
		// if( ! (in.flags() & std::ios::binary) ){
		// 	std::cerr << "XformMap::save must be binary ostream" << std::endl;
		// 	return false;
		// }
		size_t s;
		in.read((char*)&s,sizeof(size_t));
		char buf[9999];
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read(buf,s);
		std::cout << "XformHashNeighbors load, description: " << std::endl;
		std::cout << std::string(buf).substr(37,s-80) << std::endl;
		in.read( (char*)&s, sizeof(size_t) );
		// std::cout << "SIZE IN " << s << std::endl;
		for(int i = 0; i < 9999; ++i) buf[i] = 0;
		in.read( buf, s );
		// std::cout << "BUF: " <<  buf << " " << std::string(buf) <<  std::endl;
		if( hasher_.name() != std::string(buf) ){
			std::cerr << "XformHashNeighbors::load, hasher type mismatch, expected " << hasher_.name() << " got "  << buf << std::endl;
			return false;
		}
		XformHash test;
		in.read( (char*)&test, sizeof(XformHash) );
		if( test != hasher_ ){
			std::cerr << "XformHashNeighbors::load, hasher mismatch!" << std::endl;
			return false;			
		}
		Float cart_bound,ang_bound,quat_bound; int nsamp;
		in.read( (char*)&cart_bound, sizeof(Float) );		
		if( cart_bound != cart_bound_ ){
			std::cerr << "XformHashNeighbors::load, cart_bound mismatch, expected " << cart_bound_ << " got "  << cart_bound << std::endl;
			return false;
		}
		in.read( (char*)&ang_bound, sizeof(Float) );
		if( ang_bound != ang_bound_ ){
			std::cerr << "XformHashNeighbors::load, ang_bound mismatch, expected " << ang_bound_ << " got "  << ang_bound << std::endl;
			return false;
		}
		in.read( (char*)&quat_bound, sizeof(Float) );
		if( quat_bound != quat_bound_ ){
			std::cerr << "XformHashNeighbors::load, quat_bound mismatch, expected " << quat_bound_ << " got "  << quat_bound << std::endl;
			return false;
		}
		in.read( (char*)&nsamp, sizeof(int) );
		if( nsamp != nsamp_ ){
			std::cerr << "XformHashNeighbors::load, nsamp mismatch, expected " << nsamp_ << " got "  << nsamp << std::endl;
			return false;
		}

		size_t nentries;
		in.read( (char*)&nentries, sizeof(size_t) );
		for(size_t i = 0; i < nentries; ++i){
			Key key;
			in.read( (char*)&key, sizeof(Key) );
			std::vector<Key> & v = ori_cache_.insert( std::make_pair(key,std::vector<Key>()) ).first->second;
			in.read( (char*)&s, sizeof(size_t) );
			for(size_t k = 0; k < s; ++k){
				Key key2;
				in.read( (char*)&key2, sizeof(Key) );
				v.push_back(key2);
			}
		}

		// shifts computed by c'tor

		// std::cout << "SIZE IN " << map_.size() << std::endl;
		return true;

	}



};




}}}

#endif
