#include <gtest/gtest.h>

#include "riflib/rosetta/score/AnalyticEvaluation.hh"

#include "riflib/rosetta/score/EtableParams_init.hh"


namespace scheme { namespace rosetta { namespace score { namespace test {

using std::cout;
using std::endl;

TEST(RosettaEtable,init_params){
	EtableParams<float> params;
	EtableParamsInit::init_EtableParams(params);
	ASSERT_EQ( params.size(), 325ul );

	double sol,atr,rep;
	double dis,dis2,inv_dis2;
	int pair;

	pair = 13; dis = 3.7; dis2 = dis*dis; inv_dis2 = 1.0/dis2;
	lj_evaluation( params[pair], dis, dis2, inv_dis2, atr, rep);
	lk_evaluation( params[pair], dis, inv_dis2, sol ); 
	ASSERT_NEAR( atr, -0.13149315, 0.0001 );
	ASSERT_NEAR( rep,  0.0, 0.0001 );
	ASSERT_NEAR( sol,  0.161156, 0.0001 );
}


}}}}

