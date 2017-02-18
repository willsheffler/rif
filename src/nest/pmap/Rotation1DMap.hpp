#ifndef INCLUDED_scheme_nest_maps_Rotation1DMap_HH
#define INCLUDED_scheme_nest_maps_Rotation1DMap_HH

#include "util/SimpleArray.hpp"

#include <Eigen/Dense>

#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme { namespace nest { namespace pmap {

	template<
		int DIM=1,
		class Value=Eigen::Matrix3d,
		class Index=uint64_t,
		class Float=double
	>
	struct Rotation1DMap {
		BOOST_STATIC_ASSERT_MSG(DIM==1,"Rotation1D DIM must be == 1");

		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef util::SimpleArray<DIM,Index> Indices;
		typedef util::SimpleArray<DIM,Float> Params;		
		typedef Eigen::Matrix<Float,3,1> V;
		typedef Eigen::Matrix<Float,3,3> M;

		Index nside;
		Float lb, ub;
		V axis, flipaxis;
		M flip;


		Rotation1DMap() : nside(1), lb(-M_PI), ub(M_PI), axis(0,0,1), flipaxis(0,0,0) {}
		Rotation1DMap( Params p, Params q, Indices i ) : lb(p[0]), ub(q[0]), nside(i[0]), axis(0,0,1), flipaxis(0,0,0) {}

		void set_axis(V a){ axis = a.normalized(); }

		void set_flip_axis(V a){
			BOOST_ASSERT( fabs(axis.dot(a)) < 0.00001 );
			a.normalize();
			flipaxis = a;
			flip = Eigen::AngleAxis<Float>(M_PI,a).toRotationMatrix();
		}

		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Index resl,
			Value & value
		) const {
			BOOST_ASSERT( fabs(axis.norm()-1.0) < 0.00001 );
			BOOST_ASSERT( lb < ub );
			bool doflip = cell_index >= nside;
			cell_index -= doflip ? nside : 0;
			BOOST_ASSERT( cell_index < nside );
			Float ang = (Float)cell_index + params[0];
			ang *= (ub-lb) / (Float)nside;
			ang += lb;
			value = Eigen::AngleAxis<Float>( ang, axis ).matrix();
			value = value * ( doflip? flip : M::Identity() );
			BOOST_ASSERT( lb <= ang && ang <= ub );
			// std::cout << "SET " << doflip <<  std::endl;//axis.transpose() << " " << ang*180/M_PI << std::endl;
			// std::cout << value << std::endl;
			return true;
		}

		void get_axis_angle_and_flip( Value value, Eigen::Vector3d & ax, Float & ang, bool & doflip ) const {
			Eigen::AngleAxis<Float> aa(value);
			ang = aa.angle();
			ax = aa.axis();
			if( ax.dot(axis) < 0 ){ ax = -ax; ang = -ang; }
			doflip = false;
			if( fabs(ang) > 0.00001 && ax.dot(axis) < 0.999 ){
				doflip = true;
				value = value * flip.transpose();
				Eigen::AngleAxis<Float> aa(value);
				ax = aa.axis();
				ang = aa.angle();
				if( ax.dot(axis) < 0 ){ ax = -ax; ang = -ang; }
				// std::cout << ax.transpose() << " / " << axis.transpose() << " / " << ang*180/M_PI << std::endl;
				// std::cout << value << std::endl;
			}

		}

		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value0,
			Index /*resl*/,
			Params & params,
			Index & cell_index
		) const {
			M value = value0;
			BOOST_ASSERT( fabs(axis.norm()-1.0) < 0.00001 );
			BOOST_ASSERT( lb < ub );

			Float ang;
			bool doflip;
			Eigen::Vector3d vax;
			get_axis_angle_and_flip( value0, vax, ang, doflip );

			BOOST_ASSERT( fabs(ang) < 0.00001 || vax.dot(axis) > 0.99999 );
			if( ang >  M_PI ) ang -= 2.0*M_PI;
			if( ang < -M_PI ) ang += 2.0*M_PI;			
			// std::cout << lb << " " << ang << " " << ub << " /  " << nside << " " << aa.axis().transpose() << " " << aa.angle() << std::endl;
			BOOST_ASSERT( lb <= ang && ang <= ub );
			Float p = (ang-lb)/(ub-lb) * (Float)nside;
			cell_index = (Index)p;
			params[0] = p - cell_index;
			BOOST_ASSERT( cell_index < nside );
			cell_index += doflip ? nside : 0;
			BOOST_ASSERT( 0 <= params[0] && params[0] <= 1 );
			return true;
		}

		///@brief get parameter space repr of Value for particular cell
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Params & params
		) const {
			std::cerr << "Not Implemented" << std::endl;
			std::exit(-1);			
		}

		///@brief return the cell_index of neighboring cells within radius of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(Value const & value, Float radius, OutIter out) const {
			std::cerr << "Not Implemented" << std::endl;
			std::exit(-1);
		}

		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const {
			return (ub-lb) / nside / (Float)(1<<resl);
		}

		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const {
			return (ub-lb) / nside / (Float)(1<<resl);
		}

		///@brief cell size
		Index num_cells() const { return nside * ( flipaxis==V(0,0,0)? 1 : 2 ); }

	};



}}}

#endif
