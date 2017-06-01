// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <numeric/geometry/hashing/xyzStripeHash.hh>

namespace numeric {
namespace geometry {
namespace hashing {

using std::cout;
using std::cerr;
using std::endl;

xyzStripeHash::xyzStripeHash(
	float grid_size,
    utility::vector1<Ball> const & balls
):
	grid_size_(grid_size),
	grid_size2_(grid_size*grid_size),
	grid_balls_(NULL),
	grid_stripe_(NULL)//,
	// neighbor_end_(*this)
{
	if(balls.size()) init(balls);
 }

void
xyzStripeHash::init(
	utility::vector1<Ball> const & balls
){
	if( balls.size() < 1 || balls.size() > 65535 ){
		// std::cout << std::endl << "nballs: " << balls.size() << std::endl;
		// utility_exit_with_message("xyzStripeHash requires 0 < N < 65535 balls!");
		nballs_ = 0;
		grid_size_ = 0;
		grid_size2_ = 0;
		xmx_ = 0;
		ymx_ = 0;
		zmx_ = 0;
		return;
	}
	if(grid_size_ == 0.0 ){
		BOOST_FOREACH(Ball const & b, balls) grid_size_ = max(grid_size_,2.0f*b.lj_radius());
		grid_size2_ = grid_size_*grid_size_;
	}
	nballs_ = balls.size();
	// neighbor_end_.end();

	if( grid_size_ <= 0.0 ) utility_exit_with_message("grid_size_ <= 0");

	float xmn= 9e9,ymn= 9e9,zmn= 9e9;
	float xmx=-9e9,ymx=-9e9,zmx=-9e9;
	for(int i = 1; i <= nballs_; ++i) {
		xmn = numeric::min(xmn,balls[i].x());
		ymn = numeric::min(ymn,balls[i].y());
		zmn = numeric::min(zmn,balls[i].z());
		xmx = numeric::max(xmx,balls[i].x());
		ymx = numeric::max(ymx,balls[i].y());
		zmx = numeric::max(zmx,balls[i].z());
	}

	xdim_ = static_cast<int>((xmx-xmn+0.0001)/grid_size_+0.999999);
	ydim_ = static_cast<int>((ymx-ymn+0.0001)/grid_size_+0.999999);
	zdim_ = static_cast<int>((zmx-zmn+0.0001)/grid_size_+0.999999);
	assert(xdim_ < 9999); assert(ydim_ < 9999); assert(zdim_ < 9999);
	int const gsize = xdim_*ydim_*zdim_;
	ushort2 *gindex  = new ushort2[gsize];
	ushort2 *gstripe = new ushort2[gsize];
	for(int i = 0; i < gsize; ++i) { gindex[i].y = 0; gindex[i].x = 0; }
	//TR<<"atom "<<nballs_<<" grid1 "<<xdim_*ydim_*zdim_<<" "<<xdim_<<" "<<ydim_<<" "<<zdim_<<std::endl;

	for(int i = 1; i <= nballs_; ++i) {
		int ix = static_cast<int>((balls[i].x()-xmn/*+FUDGE*/)/grid_size_);
		int iy = static_cast<int>((balls[i].y()-ymn/*+FUDGE*/)/grid_size_);
		int iz = static_cast<int>((balls[i].z()-zmn/*+FUDGE*/)/grid_size_);
		assert(ix >= 0); assert(iy >= 0); assert(iz >= 0); assert(ix < xdim_); assert(iy < ydim_); assert(iz < zdim_);
		int ig = ix+xdim_*iy+xdim_*ydim_*iz;
		assert(ig>=0);assert(ig<9999999);
		++(gindex[ig].y);
	}
	for(int i = 1; i < gsize; ++i) gindex[i].x = gindex[i-1].x + gindex[i-1].y;
	for(int i = 1; i < gsize; ++i) gindex[i].y = gindex[i  ].x + gindex[i  ].y;
	for( int iz = 0; iz < zdim_; ++iz) for( int iy = 0; iy < ydim_; ++iy) for( int ix = 0; ix < xdim_; ++ix) {
				int const ixl = (int)numeric::max(      0 ,(int)ix-1 );
				int const ixu =       numeric::min(xdim_-1u,     ix+1u);
				int const ig0 = xdim_*iy+xdim_*ydim_*iz;
				gstripe[ix+ig0].x = gindex[ixl+ig0].x;
				gstripe[ix+ig0].y = gindex[ixu+ig0].y;
			}
	grid_stripe_ = gstripe;
	// for(int iz = 0; iz < zdim_; ++iz) for(int iy = 0; iy < ydim_; ++iy) for(int ix = 0; ix < xdim_; ++ix) {
	//       int i = ix+xdim_*iy+xdim_*ydim_*iz;
	//       TR<<ix<<" "<<iy<<" "<<iz<<" "<<I(3,gindex[i].x)<<" "<<I(3,gindex[i].y) <<" "<<I(3,grid_stripe_[i].x)<<" "<<I(3,grid_stripe_[i].y)<<std::endl;
	//     }
	Ball *gatom = new Ball[nballs_+4]; // space for 4 overflow balls
	for(int i=0;i<4;++i) {gatom[nballs_+i].x()=9e9;gatom[nballs_+i].y()=9e9;gatom[nballs_+i].z()=9e9;}
	ushort *gridc = new ushort[gsize];
	for(int i = 0; i < gsize; ++i) gridc[i] = 0;
	for(int i = 1; i <= nballs_; ++i) {
		int const ix = static_cast<int>((balls[i].x()-xmn/*+FUDGE*/)/grid_size_);
		int const iy = static_cast<int>((balls[i].y()-ymn/*+FUDGE*/)/grid_size_);
		int const iz = static_cast<int>((balls[i].z()-zmn/*+FUDGE*/)/grid_size_);
		int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
		int const idx = gindex[ig].x + gridc[ig];
		gatom[ idx ].x() = balls[i].x()-xmn/*+FUDGE*/;
		gatom[ idx ].y() = balls[i].y()-ymn/*+FUDGE*/;
		gatom[ idx ].z() = balls[i].z()-zmn/*+FUDGE*/;
		gatom[ idx ].atype( balls[i].atype() );
		gatom[ idx ].resid_ = balls[i].resid_;
		gatom[ idx ].atomno_ = balls[i].atomno_;
		++(gridc[ig]);
	}
	grid_balls_ = gatom;
	translation_.x() =/* FUDGE*/ - xmn;
	translation_.y() =/* FUDGE*/ - ymn;
	translation_.z() =/* FUDGE*/ - zmn;
	xmx_ = xmx-xmn/*+FUDGE*/+grid_size_;
	ymx_ = ymx-ymn/*+FUDGE*/+grid_size_;
	zmx_ = zmx-zmn/*+FUDGE*/+grid_size_;
	// for(int iz = 0; iz < zdim(); ++iz) for(int iy = 0; iy < ydim(); ++iy) for(int ix = 0; ix < xdim(); ++ix) {
	//       int i = ix+xdim_*iy+xdim_*ydim_*iz;
	//       TR<<"GRID CELL "<<ix<<" "<<iy<<" "<<iz<<std::endl;
	//       for(int ig = gindex[i].x; ig < gindex[i].y; ++ig) {
	//       TR<<F(7,3,gatom[ig].x)<<" "<<F(7,3,gatom[ig].y)<<" "<<F(7,3,gatom[ig].z)<<std::endl;
	//     }
	//   }
	delete gridc;
	delete gindex;
 }

bool xyzStripeHash::sanity_check() const {
	using namespace ObjexxFCL::format;
	for(int ix = 0; ix < xdim_; ++ix) {
		for(int iy = 0; iy < ydim_; ++iy) {
			for(int iz = 0; iz < zdim_; ++iz) {
				//std::cout << ix << " " << iy << " " << iz << endl;
				ushort const ig  = ix+xdim_*iy+ydim_*xdim_*iz;
				ushort const igl = grid_stripe_[ig].x;
				ushort const igu = grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					// float const & x(grid_balls_[i].x);
					float const & y(grid_balls_[i].y());
					float const & z(grid_balls_[i].z());
				 // if(i==igl) std::cout << endl;
					// bool xc = grid_size_*(float)ix <= x && x <= grid_size_*(float)(ix+1);
					bool yc = grid_size_*(float)iy <= y && y <= grid_size_*(float)(iy+1);
					bool zc = grid_size_*(float)iz <= z && z <= grid_size_*(float)(iz+1);
					if(/*!xc||*/!yc||!zc) utility_exit_with_message("INSANE!");
					//std::cout<<I(2,ix)<<" "<<I(2,iy)<<" "<<I(2,iz)<<" "<<F(8,3,x)<<" "<<F(8,3,y)<<" "<<F(8,3,z)<<" "<<xc<<" "<<yc<<" "<<zc<<std::endl;
				}
			}
			return true;
		}
	}
	return true;
 }

int
xyzStripeHash::nbcount( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
	if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
	int count = 0;
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= grid_size2_ ) {
					++count;
				}
			}
		}
	}
	return count;
 }

int
xyzStripeHash::nbcount_raw( Vec const & v ) const {
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return 0; // worth it iff
	if( x > xmx_ || y > ymx_ || z > zmx_ ) return 0;                      // worth it iff
	int count = 0;
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_), iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_), static_cast<int>(iz0+2));
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= grid_size2_ ) {
					++count;
				}
			}
		}
	}
	return count;
 }

bool
xyzStripeHash::clash( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					return true;
				}
			}
		}
	}
	return false;
 }

float
xyzStripeHash::clash_amount( Vec const & v_in ) const {
	Vec const v = v_in+translation_;
	float clash_amount = grid_size2_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					cout << grid_size_-sqrt(d2) << " " << grid_size2_ << " " << a2.resi() << " " << a2.atomno() << " " <<  endl;
					clash_amount = numeric::min(d2,clash_amount);
			}
		}
	}
 }
	clash_amount = grid_size_ - sqrt(clash_amount);
	return clash_amount;
 }

bool
xyzStripeHash::clash_not_resid( Vec const & v_in, int const & resid, int const & resid2 ) const {
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ && (int)a2.resid_ != resid && (int)a2.resid_ != resid2 ) {
					return true;
				}
			}
		}
	}
	return false;
 }

bool
xyzStripeHash::clash_raw( Vec const & v ) const {
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_), static_cast<int>(iz0+2));
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 < grid_size2_ ) {
					return true;
				}
			}
		}
	}
	return false;
 }

int
xyzStripeHash::clash_check_ball( Ball const & b_in ) const {
	// Input ball coordinates are in global coordinate space,
	// translate ball into hash frame.
	Vec b_in_hash_space = b_in.xyz() + translation_;
	float x = b_in_hash_space.x(); float y = b_in_hash_space.y(); float z = b_in_hash_space.z();

	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return false; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return false; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1, static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_), static_cast<int>(iz0+2));
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				
				if( d2 < (a2.lj_radius() + b_in.lj_radius()) ) {
					// Return type is 1-based ball number.
					return i + 1;
				}
			}
		}
	}

	return 0;
 }

bool xyzStripeHash::clash_check_residue_pairs(
		utility::vector1<Ball> const & test_balls,
		std::map<Size, Size> & residue_pairs
) const
{
	bool found_clash = false;

	for (Size i = 1; i <= test_balls.size(); i++)
	{
		if (residue_pairs.count(test_balls[i].resi()))
		{
			continue;
		}

		int clashing_ball = clash_check_ball(test_balls[i]);

		if (clashing_ball)
		{
			residue_pairs[test_balls[i].resi()] = ball(clashing_ball).resi();
			found_clash = true;
		}
	}

	return found_clash;
}

void
xyzStripeHash::fill_pairs(
	Vec const & v_in,
	int const & ir,
	utility::vector1<std::pair<int,int> > & pairs,
	float maxd2
) const {
	if(0.0==maxd2) maxd2=grid_size2_;
	Vec const v = v_in+translation_;
	float x = v.x(); float y = v.y(); float z = v.z();
	if( x < -grid_size_ || y < -grid_size_ || z < -grid_size_ ) return; // worth it iff
	if( x > xmx_        || y > ymx_        || z > zmx_        ) return; // worth it iff
	int const ix   = (x<0) ? 0 : numeric::min(xdim_-1,static_cast<int>(x/grid_size_));
	int const iy0  = (y<0) ? 0 : static_cast<int>(y/grid_size_);
	int const iz0  = (z<0) ? 0 : static_cast<int>(z/grid_size_);
	int const iyl = numeric::max(0,iy0-1);
	int const izl = numeric::max(0,iz0-1);
	int const iyu = numeric::min(static_cast<int>(ydim_),iy0+2);
	int const izu = numeric::min(static_cast<int>(zdim_),static_cast<int>(iz0)+2);
	for(int iy = iyl; iy < iyu; ++iy) {
		for(int iz = izl; iz < izu; ++iz) {
			int const ig = ix+xdim_*iy+xdim_*ydim_*iz;
			assert(ig < xdim_*ydim_*zdim_);
			assert(ix < xdim_);
			assert(iy < ydim_);
			assert(iz < zdim_);
			int const & igl = grid_stripe_[ig].x;
			int const & igu = grid_stripe_[ig].y;
			for(int i = igl; i < igu; ++i) {
				Ball const a2 = grid_balls_[i];
				float const d2 = (x-a2.x())*(x-a2.x()) + (y-a2.y())*(y-a2.y()) + (z-a2.z())*(z-a2.z());
				if( d2 <= maxd2 ) {
					pairs.push_back(std::make_pair(ir,a2.resi()));
				}
			}
		}
	}
 }

std::string xyzStripeHash::debug_pdb(Xform const & x) const {
	using namespace ObjexxFCL::format;
	std::ostringstream out;
	int atomno = 0;
	for(int i = 0; i < nballs_; ++i){
		Ball const b = grid_balls_[i];
		Vec v = x*(b.xyz() - translation_);
		++atomno;
		out<<"ATOM  "<<I(5,atomno)<<' '<<" XSH"<<' '<<"XSH"<<' '<<'~'<<I(4,atomno/100)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
	}
	return out.str();

}

void Ball::atype(Size at){
	if( at < 52 ) atype_ = at;
}


// NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME
// CNH2     C    2.0000    0.1200    0.0000    3.5000   14.7000
// COO      C    2.0000    0.1200   -1.4000    3.5000    8.3000 
// CH1      C    2.0000    0.0486   -0.2500    3.5000   23.7000
// CH2      C    2.0000    0.1142    0.5200    3.5000   22.4000
// CH3      C    2.0000    0.1811    1.5000    3.5000   30.0000
// aroC     C    2.0000    0.1200    0.0800    3.5000   18.4000 
// Ntrp     N    1.7500    0.2384   -8.9000    3.5000    4.4000 
// Nhis     N    1.7500    0.2384   -4.0000    3.5000    4.4000 
// NH2O     N    1.7500    0.2384   -7.8000    3.5000   11.2000 
// Nlys     N    1.7500    0.2384  -20.0000    6.0000   11.2000 
// Narg     N    1.7500    0.2384  -10.0000    6.0000   11.2000 
// Npro     N    1.7500    0.2384   -1.5500    3.5000    0.0000 
// OH       O    1.5500    0.1591   -6.7000    3.5000   10.8000 
// ONH2     O    1.5500    0.1591   -5.8500    3.5000   10.8000 
// OOC      O    1.5500    0.2100  -10.0000    6.0000   10.8000 
// Oaro     O    1.5500    0.1591   -4.0000    3.5000   10.8000 
// S        S    1.9000    0.1600   -4.1000    3.5000   14.7000 
// Nbb      N    1.7500    0.2384   -5.0000    3.5000    4.4000 
// CAbb     C    2.0000    0.0486    1.0000    3.5000   23.7000
// CObb     C    2.0000    0.1400    1.0000    3.5000   14.7000
// OCbb     O    1.5500    0.1591   -5.0000    3.5000   10.8000 
// Phos     P    2.1500    0.5850  -24.0000    3.5000   34.8000
// Pbb      P    2.1500    0.5850  -24.0000    3.5000   34.8000 
// Hpol     H    1.0000    0.0500    0.0000    3.5000    0.0000 
// Hapo     H    1.2000    0.0500    0.0000    3.5000    0.0000
// Haro     H    1.2000    0.0500    0.0000    3.5000    0.0000
// HNbb     H    1.0000    0.0500    0.0000    3.5000    0.0000 
// HOH      O    1.4000    0.0500    0.0000    3.5000   10.8000 
// Bsp2     B    1.9800    0.0340   -4.0000    3.5000    4.4000 
// F        F    1.7100    0.0750   -6.0000    3.5000   11.5000 
// Cl      CL    2.0700    0.2400   -5.0000    3.5000   24.4000 
// Br      BR    2.2200    0.3200   -4.0000    3.5000   35.5000 
// I        I    2.3600    0.4240   -3.0000    3.5000   44.6000 
// Zn2p    ZN    1.0900    0.2500   -5.0000    3.5000    5.4000 
// Co2p    CO    1.5680    0.2500   -5.0000    3.5000   16.1483 
// Cu2p    CU    1.0978    0.2500   -5.0000    3.5000    5.5419 
// Fe2p    FE    0.7800    0.0000   -5.0000    3.5000    2.0000 
// Fe3p    FE    0.6500    0.0000   -5.0000    3.5000    1.2000 
// Mg2p    MG    1.1850    0.0150   -5.0000    3.5000    7.0000 
// Ca2p    CA    1.3670    0.1200    0.0000    2.0000   10.7000 
// Pha      P    2.1500    0.5850   -4.0000    6.0000   37.9000 
// OPha     O    1.7000    0.1200   -6.0000    1.7000   10.8000 
// OHha     O    1.7700    0.1521   -6.7700    3.5000   22.1000 
// Hha      H    0.2245    0.0460    0.0000    3.5000    2.0000 
// CO3      C    1.7000    0.0880   -2.1000    8.0000    8.3    
// OC3      O    1.6000    0.1370  -10.0000    8.0000   10.8    
// Si       Si   1.8532    0.4020    3.9000    2.2000   34.7    
// OSi      O    1.7766    0.0600   -3.9000    1.6000   20.4    
// Oice     O    1.6000    0.1591   -6.7000    3.5000   14.3    
// Hice     H    0.8000    0.0498    0.0000    3.500     2.0    
// Na1p    NA    1.3638    0.0469   -5.0000    3.5000   10.6000 
// K1p      K    1.7638    0.0870   -5.0000    3.5000   23.0000 

static std::string const ATOM_NAME[53] = {
"NONE",
"CNH2",
"COO",
"CH1",
"CH2",
"CH3",
"aroC",
"Ntrp",
"Nhis",
"NH2O",
"Nlys",
"Narg",
"Npro",
"OH",
"ONH2",
"OOC",
"Oaro",
"S",
"Nbb",
"CAbb",
"CObb",
"OCbb",
"Phos",
"Pbb",
"Hpol",
"Hapo",
"Haro",
"HNbb",
"HOH",
"Bsp2",
"F",
"Cl",
"Br",
"I",
"Zn2p",
"Co2p",
"Cu2p",
"Fe2p",
"Fe3p",
"Mg2p",
"Ca2p",
"Pha",
"OPha",
"OHha",
"Hha",
"CO3",
"OC3",
"Si",
"OSi",
"Oice",
"Hice",
"Na1p",
"K1p"
};

std::string 
Ball::atom_type_name() const {
	return ATOM_NAME[atype_];
}

static float const ATOM_LJ_RADIUS[53] = {
1.0,
2.0000,
2.0000,
2.0000,
2.0000,
2.0000,
2.0000,
1.7500,
1.7500,
1.7500,
1.7500,
1.7500,
1.7500,
1.5500,
1.5500,
1.5500,
1.5500,
1.9000,
1.7500,
2.0000,
2.0000,
1.5500,
2.1500,
2.1500,
1.0000,
1.2000,
1.2000,
1.0000,
1.4000,
1.9800,
1.7100,
2.0700,
2.2200,
2.3600,
1.0900,
1.5680,
1.0978,
0.7800,
0.6500,
1.1850,
1.3670,
2.1500,
1.7000,
1.7700,
0.2245,
1.7000,
1.6000,
1.8532,
1.7766,
1.6000,
0.8000,
1.3638,
1.7638
};

float 
Ball::lj_radius() const {
	return ATOM_LJ_RADIUS[atype_];
}


static float const ATOM_LJ_WDEPTH[53] = {
	0.0,
   0.1200,
   0.1200,
   0.0486,
   0.1142,
   0.1811,
   0.1200,
   0.2384,
   0.2384,
   0.2384,
   0.2384,
   0.2384,
   0.2384,
   0.1591,
   0.1591,
   0.2100,
   0.1591,
   0.1600,
   0.2384,
   0.0486,
   0.1400,
   0.1591,
   0.5850,
   0.5850,
   0.0500,
   0.0500,
   0.0500,
   0.0500,
   0.0500,
   0.0340,
   0.0750,
   0.2400,
   0.3200,
   0.4240,
   0.2500,
   0.2500,
   0.2500,
   0.0000,
   0.0000,
   0.0150,
   0.1200,
   0.5850,
   0.1200,
   0.1521,
   0.0460,
  0.0880 ,
  0.1370 ,
  0.4020 ,
  0.0600 ,
  0.1591 ,
  0.0498 ,
   0.0469,
   0.0870
};

float 
Ball::lj_depth() const {
	return ATOM_LJ_WDEPTH[atype_];
}

static float const ATOM_LK_DGFREE[53] = {
	0.0,
   0.0000,
  -1.4000,
  -0.2500,
   0.5200,
   1.5000,
   0.0800,
  -8.9000,
  -4.0000,
  -7.8000,
 -20.0000,
 -10.0000,
  -1.5500,
  -6.7000,
  -5.8500,
 -10.0000,
  -4.0000,
  -4.1000,
  -5.0000,
   1.0000,
   1.0000,
  -5.0000,
 -24.0000,
 -24.0000,
   0.0000,
   0.0000,
   0.0000,
   0.0000,
   0.0000,
  -4.0000,
  -6.0000,
  -5.0000,
  -4.0000,
  -3.0000,
  -5.0000,
  -5.0000,
  -5.0000,
  -5.0000,
  -5.0000,
  -5.0000,
   0.0000,
  -4.0000,
  -6.0000,
  -6.7700,
   0.0000,
 -2.1000 ,
 -10.0000 ,
  3.9000 ,
 -3.9000 ,
 -6.7000 ,
  0.0000 ,
  -5.0000,
  -5.0000
};


float 
Ball::lk_dgfree() const {
	return ATOM_LK_DGFREE[atype_];
}


static float const ATOM_LK_LAMBDA [53] = {
	1.0,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   6.0000 ,
   6.0000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   6.0000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   3.5000 ,
   2.0000 ,
   6.0000 ,
   1.7000 ,
   3.5000 ,
   3.5000 ,
   8.0000 ,
   8.0000 ,
   2.2000 ,
   1.6000 ,
   3.5000 ,
   3.500  ,
   3.5000 ,
   3.5000 
};

float 
Ball::lk_lambda() const {
	return ATOM_LK_LAMBDA[atype_];
}

static float const ATOM_LK_VOLUME[53] = {
	0.0,
  14.7000,
   8.3000 ,
  23.7000,
  22.4000,
  30.0000,
  18.4000 ,
   4.4000 ,
   4.4000 ,
  11.2000 ,
  11.2000 ,
  11.2000 ,
   0.0000 ,
  10.8000 ,
  10.8000 ,
  10.8000 ,
  10.8000 ,
  14.7000 ,
   4.4000 ,
  23.7000,
  14.7000,
  10.8000 ,
  34.8000,
  34.8000 ,
   0.0000 ,
   0.0000,
   0.0000,
   0.0000 ,
  10.8000 ,
   4.4000 ,
  11.5000 ,
  24.4000 ,
  35.5000 ,
  44.6000 ,
   5.4000 ,
  16.1483 ,
   5.5419 ,
   2.0000 ,
   1.2000 ,
   7.0000 ,
  10.7000 ,
  37.9000 ,
  10.8000 ,
  22.1000 ,
   2.0000 ,
   8.3    ,
  10.8    ,
  34.7    ,
  20.4    ,
  14.3    ,
   2.0    ,
  10.6000 ,
  23.0000 
};


float 
Ball::lk_volume() const {
	return ATOM_LK_VOLUME[atype_];
}



 } // namespace hashing
 } // namespace geometry
 } // namespace numeric

