#ifndef INCLUDED_scheme_numeric_geom_4d_HH
#define INCLUDED_scheme_numeric_geom_4d_HH

#include "riflib/numeric/util.hh"
#include "riflib/util/SimpleArray.hh"
#include <cmath>

namespace scheme { namespace numeric {

static bool is_not_0(double a){ return fabs(a) > 0.00000001; }
static bool is_not_0(float  a){ return fabs(a) > 0.0001; } // need something higher than numeric_limits::epsilon
template<class Q> Q to_half_cell(Q const & q){
	return is_not_0(q.w()) ? ( q.w()>0 ? q : Q(-q.w(),-q.x(),-q.y(),-q.z()) ) : (
				// q
	       		is_not_0(q.x()) ? ( q.x()>0 ? q : Q(-q.w(),-q.x(),-q.y(),-q.z()) ) : (
		       		is_not_0(q.y()) ? ( q.y()>0 ? q : Q(-q.w(),-q.x(),-q.y(),-q.z()) ) : (
			       		( q.z()>0 ? q : Q(-q.w(),-q.x(),-q.y(),-q.z()) )
			       	)
	       		)
           );

}



template<class Float>
static Float const *
get_raw_48cell(){
	static Float const r = sqrt(2)/2;
	static Float const h = 0.5;
	static Float const raw48[48*4] = {
	        1, 0, 0, 0, // 0
	        0, 1, 0, 0, // 1
	        0, 0, 1, 0, // 2
	        0, 0, 0, 1, // 3
	       -1, 0, 0, 0, // 4
	        0,-1, 0, 0, // 5
	        0, 0,-1, 0, // 6
	        0, 0, 0,-1, // 7
	        h, h, h, h, // 8
	        h, h, h,-h, // 9
	        h, h,-h, h, // 10
	        h, h,-h,-h, // 11
	        h,-h, h, h, // 12
	        h,-h, h,-h, // 13
	        h,-h,-h, h, // 14
	        h,-h,-h,-h, // 15
	       -h, h, h, h, // 16
	       -h, h, h,-h, // 17
	       -h, h,-h, h, // 18
	       -h, h,-h,-h, // 19
	       -h,-h, h, h, // 20
	       -h,-h, h,-h, // 21
	       -h,-h,-h, h, // 22
	       -h,-h,-h,-h, // 23
	        r, r, 0, 0, // 24
	        r, 0, r, 0, // 25
	        r, 0, 0, r, // 26
	        0, r, r, 0, // 27
	        0, r, 0, r, // 28
	        0, 0, r, r, // 29
	       -r, r, 0, 0, // 30
	       -r, 0, r, 0, // 31
	       -r, 0, 0, r, // 32
	        0,-r, r, 0, // 33
	        0,-r, 0, r, // 34
	        0, 0,-r, r, // 35
	        r,-r, 0, 0, // 36
	        r, 0,-r, 0, // 37
	        r, 0, 0,-r, // 38
	        0, r,-r, 0, // 39
	        0, r, 0,-r, // 40
	        0, 0, r,-r, // 41
	       -r,-r, 0, 0, // 42
	       -r, 0,-r, 0, // 43
	       -r, 0, 0,-r, // 44
	        0,-r,-r, 0, // 45
	        0,-r, 0,-r, // 46
	        0, 0,-r,-r  // 47
	    };
	return raw48;
}


template<class V4,class Index> void
get_cell_48cell( V4 const & quat, Index & cell ){
	typedef typename V4::Scalar Float;
	V4 const quat_pos = quat.cwiseAbs();
	V4 tmpv = quat_pos;

	Float hyperface_dist; // dist to closest face
	Index hyperface_axis; // closest hyperface-pair
	Float edge_dist;      // dist to closest edge
	Index edge_axis_1;    // first axis of closest edge
	Index edge_axis_2;    // second axis of closest edge
	Float corner_dist;    // dist to closest corner

	numeric::max2( quat_pos, hyperface_dist, edge_dist, hyperface_axis, edge_axis_2 );
	edge_dist = sqrt(2)/2 * ( hyperface_dist + edge_dist );
	corner_dist = quat_pos.sum()/2;
	edge_axis_1 = hyperface_axis<edge_axis_2 ? hyperface_axis : edge_axis_2;
	edge_axis_2 = hyperface_axis<edge_axis_2 ? edge_axis_2 : hyperface_axis;
	assert( edge_axis_1 < edge_axis_2 );

	// cell if closest if of form 1000 (add 4 if negative)
	Index facecell = hyperface_axis | (quat[hyperface_axis]<0 ? 4 : 0);

	// cell if closest is of form 1111, bitwise by ( < 0)
	Index bit0 = quat[3]<0;
	Index bit1 = quat[2]<0;
	Index bit2 = quat[1]<0;
	Index bit3 = quat[0]<0;						
	Index cornercell = bit0 | bit1<<1 | bit2<<2 | bit3<<3;

	// cell if closest is of form 1100
	Index perm_shift[3][4] = { {9,0,1,2}, {0,9,3,4}, {1,3,9,5} };
	Index sign_shift = (quat[edge_axis_1]<0)*1*6 + (quat[edge_axis_2]<0)*2*6;
	Index edgecell = sign_shift + perm_shift[edge_axis_1][edge_axis_2];

	// pick case 1000 1111 1100 without if statements
	Index swtch; util::SimpleArray<3,Float>(hyperface_dist,corner_dist,edge_dist).maxCoeff(&swtch);
	cell = swtch==0 ? facecell : (swtch==1 ? cornercell+8 : edgecell+24);
		// this is slower !?!
		// Float mx = std::max(std::max(hyperface_dist,corner_dist),edge_dist);
		// cell2[i] = hyperface_dist==mx ? facecell : (corner_dist==mx ? cornercell+8 : edgecell+24);
}


template<class Float>
static Float const *
get_raw_48cell_half(){
	static Float const r = sqrt(2)/2;
	static Float const h = 0.5;
	static Float const raw48[24*4] = {
	        1, 0, 0, 0, // 0
	        0, 1, 0, 0, // 1
	        0, 0, 1, 0, // 2
	        0, 0, 0, 1, // 3
	        h, h, h, h, // 8
	       -h, h, h, h, // 10
	        h,-h, h, h, // 12
	       -h,-h, h, h, // 14
	        h, h,-h, h, // 16
	       -h, h,-h, h, // 18
	        h,-h,-h, h, // 20
	       -h,-h,-h, h, // 22
	        r, r, 0, 0, // 24
	        r, 0, r, 0, // 25
	        r, 0, 0, r, // 26
	        0, r, r, 0, // 27
	        0, r, 0, r, // 28
	        0, 0, r, r, // 29
	       -r, r, 0, 0, // 30
	       -r, 0, r, 0, // 31
	       -r, 0, 0, r, // 32
	        0,-r, r, 0, // 33
	        0,-r, 0, r, // 34
	        0, 0,-r, r  // 35
	    };
	return raw48;
}


template<class V4,class Index> void
get_cell_48cell_half( V4 const & quat, Index & cell ){
	typedef typename V4::Scalar Float;
	V4 const quat_pos = quat.cwiseAbs();
	V4 tmpv = quat_pos;

	Float hyperface_dist; // dist to closest face
	Index hyperface_axis; // closest hyperface-pair
	Float edge_dist;      // dist to closest edge
	Index edge_axis_1;    // first axis of closest edge
	Index edge_axis_2;    // second axis of closest edge
	Float corner_dist;    // dist to closest corner

	// std::cout << quat_pos.transpose() << std::endl;
	numeric::max2( quat_pos, hyperface_dist, edge_dist, hyperface_axis, edge_axis_2 );
	edge_dist = sqrt(2)/2 * ( hyperface_dist + edge_dist );
	corner_dist = quat_pos.sum()/2;
	// std::cout << hyperface_axis << " " << edge_axis_2 << std::endl;
	edge_axis_1 = hyperface_axis<edge_axis_2 ? hyperface_axis : edge_axis_2;
	edge_axis_2 = hyperface_axis<edge_axis_2 ? edge_axis_2 : hyperface_axis;
	assert( edge_axis_1 < edge_axis_2 );

	// cell if closest if of form 1000 (add 4 if negative)
	Index facecell = hyperface_axis;// | (quat[hyperface_axis]<0 ? 4 : 0);

	// cell if closest is of form 1111, bitwise by ( < 0)
	Index bit0 = quat[0]<0;
	Index bit1 = quat[1]<0;
	Index bit2 = quat[2]<0;
	Index cornercell = quat[3] > 0 ? bit0 | bit1<<1 | bit2<<2 : (!bit0) | (!bit1)<<1 | (!bit2)<<2;

	// cell if closest is of form 1100
	Index perm_shift[3][4] = { {9,0,1,2}, {0,9,3,4}, {1,3,9,5} };
	Index sign_shift = ( quat[edge_axis_1]<0 != quat[edge_axis_2]<0 )*1*6;
	Index edgecell = sign_shift + perm_shift[edge_axis_1][edge_axis_2];

	// pick case 1000 1111 1100 without if statements
	Index swtch; util::SimpleArray<3,Float>(hyperface_dist,corner_dist,edge_dist).maxCoeff(&swtch);
	cell = swtch==0 ? facecell : (swtch==1 ? cornercell+4 : edgecell+12);
		// this is slower !?!
		// Float mx = std::max(std::max(hyperface_dist,corner_dist),edge_dist);
		// cell2[i] = hyperface_dist==mx ? facecell : (corner_dist==mx ? cornercell+8 : edgecell+24);

}


}}

#endif
