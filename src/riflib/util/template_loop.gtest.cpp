#include <gtest/gtest.h>
#include "riflib/util/template_loop.hpp"
#include <iostream>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "riflib/util/SimpleArray.hpp"

namespace scheme {
namespace util {

using std::cout;
using std::endl;

void testfunc(SimpleArray<3,int> const & /*a*/){}

TEST(NESTED_FOR,interface_test){
    SimpleArray<3,int> lb(0,0,3), ub(4,2,3);
    NESTED_FOR<3>(lb,ub,testfunc);
}

template<class Indices>
struct TestFun {
    size_t sum,ncalls;
    TestFun() : sum(0),ncalls(0) {}
    void func_to_call(int dummyi, int dummyj, Indices const & indices){
        // std::cerr << "func_to_call " << dummyi << " " << dummyj << " " << indices.transpose() << std::endl;
        sum += indices.sum() + dummyi + dummyj;
        ++ncalls;
    }
};

TEST(NESTED_FOR,boost_bind_functor){
    static const size_t DIM=2;
    typedef SimpleArray<DIM,int> IDX;
    IDX lb(0,0),ub(4,2);
    {
        TestFun<IDX> mytest;
        boost::function<void(IDX)> functor = boost::bind( &TestFun<IDX>::func_to_call, &mytest, 0, 0, _1 );
        NESTED_FOR<DIM>( lb, ub, functor );
        ASSERT_EQ( mytest.ncalls, (ub-lb+1).prod() ); 
        ASSERT_EQ( mytest.sum   , 45 );
    }
    {
        TestFun<IDX> mytest;
        boost::function<void(IDX)> functor = boost::bind( &TestFun<IDX>::func_to_call, &mytest, 18, 329, _1 );
        NESTED_FOR<DIM>( lb, ub, functor );
        ASSERT_EQ( mytest.ncalls, (ub-lb+1).prod() ); 
        ASSERT_EQ( mytest.sum   , 5250 );
    }
}


template<class Indices>
struct RecordCalls {
    // std::vector<Indices,Eigen::aligned_allocator<Indices> > calls;
    std::vector<Indices> calls;
    void operator()(Indices const & i){
        calls.push_back(i);
    }
};


TEST(NESTED_FOR,check_all_output){
    {
        static const size_t DIM=1;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0),ub(4);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }
    }


    {
        static const size_t DIM=2;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10),ub(4,-4 );
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}
    }



    {
        static const size_t DIM=3;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3),ub(4,-4 ,3);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}
    }


    {
        static const size_t DIM=4;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6),ub(4,-4 ,7,8);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}
    }

    {
        static const size_t DIM=5;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2),ub(4,-4 ,7,8,4);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}
    }

    {
        static const size_t DIM=6;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2,7),ub(4,-4 ,7,8,4,8);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[5] = lb[5]; idx[5] <= ub[5]; ++idx[5]){        
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}}
    }

    {
        static const size_t DIM=7;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2,7,0),ub(4,-4 ,7,8,4,8,1);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[6] = lb[6]; idx[6] <= ub[6]; ++idx[6]){        
        for(idx[5] = lb[5]; idx[5] <= ub[5]; ++idx[5]){        
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}}}
    }

    {
        static const size_t DIM=8;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2,7,0,34),ub(4,-4 ,7,8,4,8,1,35);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[7] = lb[7]; idx[7] <= ub[7]; ++idx[7]){        
        for(idx[6] = lb[6]; idx[6] <= ub[6]; ++idx[6]){        
        for(idx[5] = lb[5]; idx[5] <= ub[5]; ++idx[5]){        
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}}}}
    }

    {
        static const size_t DIM=9;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2,7,0,34,1),ub(2,-7 ,7,8,4,8,1,35,2);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[8] = lb[8]; idx[8] <= ub[8]; ++idx[8]){        
        for(idx[7] = lb[7]; idx[7] <= ub[7]; ++idx[7]){        
        for(idx[6] = lb[6]; idx[6] <= ub[6]; ++idx[6]){        
        for(idx[5] = lb[5]; idx[5] <= ub[5]; ++idx[5]){        
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}}}}}
    }

    {
        static const size_t DIM=10;
        typedef SimpleArray<DIM,int> IDX;
        IDX lb(0,-10,3,6,2,7,0,34,1,12),ub(2,-8 ,5,8,4,8,1,35,2,13);
        RecordCalls<IDX> recorder;
        NESTED_FOR<DIM>( lb, ub, recorder );
        ASSERT_EQ( recorder.calls.size(), (ub-lb+1).prod() );
        size_t count = 0;
        IDX idx;
        for(idx[9] = lb[9]; idx[9] <= ub[9]; ++idx[9]){        
        for(idx[8] = lb[8]; idx[8] <= ub[8]; ++idx[8]){        
        for(idx[7] = lb[7]; idx[7] <= ub[7]; ++idx[7]){        
        for(idx[6] = lb[6]; idx[6] <= ub[6]; ++idx[6]){        
        for(idx[5] = lb[5]; idx[5] <= ub[5]; ++idx[5]){        
        for(idx[4] = lb[4]; idx[4] <= ub[4]; ++idx[4]){        
        for(idx[3] = lb[3]; idx[3] <= ub[3]; ++idx[3]){        
        for(idx[2] = lb[2]; idx[2] <= ub[2]; ++idx[2]){        
        for(idx[1] = lb[1]; idx[1] <= ub[1]; ++idx[1]){        
        for(idx[0] = lb[0]; idx[0] <= ub[0]; ++idx[0]){
            for(size_t j=0; j < DIM; ++j){
                ASSERT_EQ( recorder.calls[count][j], idx[j] );
            }
            ++count;
        }}}}}}}}}}
    }

}


}
}
