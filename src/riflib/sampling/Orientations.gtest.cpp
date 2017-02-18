#include <gtest/gtest.h>
#include <iostream>
#include <sampling/Orientations.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <fstream>

TEST(Orientation, read_karney_datasets){
    // todo: unzip data files
    // fill in data structure instead of stream
    std::cout << "TEST read_karney_datasets" << std::endl;

    std::ifstream file( "/home/sheffler/riflib/data/orientations/karney/c48n309.grid.gz" , std::ios_base::in|std::ios_base::binary );
    boost::iostreams::filtering_streambuf<boost::iostreams::input> fin;
    fin.push(boost::iostreams::gzip_decompressor());
    fin.push(file);
    std::istream in(&fin);

    // std::string s;
    // while(in >> s) std::cout << s << endl;
    dump_karney_orientation_file(in);

    // ASSERT_TRUE(false);


}