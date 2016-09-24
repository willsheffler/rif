#include <gtest/gtest.h>

#include "scheme/nest/NEST.hh"
#include "scheme/nest/pmap/EulerAnglesMap.hh"
#include "scheme/io/dump_pdb_atom.hh"

#include <Eigen/Geometry>

#include <random>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <fstream>

namespace scheme { namespace nest { namespace pmap { namespace test {

using std::cout;
using std::endl;

TEST( EulerAnglesMap, covering ){
	using namespace Eigen;
	std::mt19937 rng((unsigned int)time(0));
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;

	cout << "EulerAnglesMap Covrad" << endl;
	int NRES = 7;
	// int const NITER = 1000000;
	int NITER = 100000;	
	NEST<3,Matrix3d,EulerAnglesMap> nest;
	for(int r = 1; r <= NRES; ++r){
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < NITER; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			Matrix3d m = nest.set_and_get( nest.get_index(q.matrix(),r) , r );
			Quaterniond qcen(m);
			// if( q.angularDistance(qcen) > maxdiff ){
			// 	RowVector3d euler; numeric::euler_angles(q.matrix(),euler);
			// 	euler[0] /= M_PI*2.0;
			// 	euler[1] /= M_PI*2.0;
			// 	euler[2] /= M_PI;
			// 	cout << r << " " << maxdiff << " " << euler << endl;
			// }
			avgdiff += q.angularDistance(qcen);
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
		}
		avgdiff /= NITER;
		// size/2 because half samples are ignored
		double volfrac = (double)nest.size(r)/2*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)nest.size(r)/2*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		printf("%2i %16llu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
			r, nest.size(r)/2, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );
		// cout << boost::format("%2i %20i %.7d %.7d") % r % nest.size(r) % (maxdiff*180.0/M_PI) % volfrac << endl;
	}

}


TEST( EulerAnglesMap, shapes ){
	std::mt19937 rng((unsigned int)time(0));
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;

	/// inspect maxdiff/avgdiff for some simple shapes
	int NITER = 10000;
	{
		util::SimpleArray<3,double> b(1,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < NITER; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			if(samp.norm() > 0.5){ --i; continue; }
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= NITER;
		cout << "sphere: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(1,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < NITER; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= NITER;
		cout << "square: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(2,1,1);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < NITER; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			avgdiff += samp.norm();
			maxdiff = std::max(samp.norm(),maxdiff);
		}
		avgdiff /= NITER;
		cout << "rect211: " << maxdiff / avgdiff << endl;
	}
	{
		util::SimpleArray<3,double> b(1,1,1), cen(0.135022,0.135022,0.135022);
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i < NITER; ++i){
			util::SimpleArray<3,double> samp(uniform(rng),uniform(rng),uniform(rng));
			samp = (samp-0.5) * b;
			if( samp.sum() < 0 ){ --i; continue; }
			avgdiff += (samp-cen).norm();
			maxdiff = std::max((samp-cen).norm(),maxdiff);
		}
		avgdiff /= NITER;
		cout << "triang: " << maxdiff / avgdiff << endl;
	}

}

// TEST( EulerAnglesMap, visualize ){
// 	using namespace Eigen;

// 	std::mt19937 rng((unsigned int)time(0));
// 	std::normal_distribution<> gauss;
// 	std::uniform_real_distribution<> uniform;

// 	Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng) );
// 	// Quaterniond qrand( 1,0,0,0 );
// 	qrand.normalize();
// 	Vector3d X = qrand*Vector3d(1,0  ,0  );
// 	Vector3d Y = qrand*Vector3d(0,1.2,0  );
// 	Vector3d Z = qrand*Vector3d(0,0  ,1.4);

// 	NEST<3,Matrix3d,EulerAnglesMap> nest;
// 	// size_t beg = 0;
// 	// while(!nest.set_state(beg,10)) beg = std::max<size_t>(uniform(rng)*(nest.size(10)-1000),0);

// 	for(size_t r = 0; r <= 8; ++r){
// 		int N = 8*8*8;
// 		// int beg = std::max( 0, (int)nest.size(r)/12 - N/2 );
// 		int beg = 0;
// 		std::ofstream out(("euler_"+boost::lexical_cast<std::string>(r)+".pdb").c_str());
// 		io::dump_pdb_atom(out,  "Z" ,0,Vector3d(0,0,0));
// 		int count1 = 0, count2 = 0;
// 		// cout << r << " " << nest.size(r) << " " << (beg>>(4*(10-r))) << endl;
// 		// continue;
// 		// for(size_t i = beg>>(4*(10-r)); i < nest.size(r); ++i){
// 		for(size_t i = beg; i < nest.size(r); ++i){		
// 			++count1;
// 			if( nest.set_state(i,r) ){
// 				++count2;
// 				if( count1 > N) break;
// 				Matrix3d m = nest.value();
// 				// cout << r << " " << i << " " << q.coeffs().transpose() << endl;
// 				Vector3d ximg = m * X;
// 				Vector3d yimg = m * Y;
// 				Vector3d zimg = m * Z;;
// 				// cout << r << " " << nest.cell_index(i,r) << endl;
// 				// out << "MODEL" << std::endl;
// 				// io::dump_pdb_atom(out,  "Z" ,count2,Vector3d(0,0,0));
// 				io::dump_pdb_atom(out, "O" ,count2,50*ximg);
// 				io::dump_pdb_atom(out, "NI",count2,50*yimg);
// 				io::dump_pdb_atom(out, "N" ,count2,50*zimg);
// 				// out << "ENDMDL" << std::endl;				 
// 			}
// 		}
// 		out.close();
// 		cout << r << " " << count2 / (double)count1 << endl;
// 	}

// }


}}}}

