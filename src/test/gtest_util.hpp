#ifndef INCLUDED_main_gtest_util_HH
#define INCLUDED_main_gtest_util_HH

#include <gtest/gtest.h>
#include <iostream>
//#include <gmock/gmock.h>

void init_gtest_tests(std::vector<std::string> const &args) {
  int argc = (int)args.size();
  char **argv = new char *[argc];
  for (int i = 0; i < argc; ++i) {
    argv[i] = new char[args[i].length() + 1];
    strcpy(argv[i], args[i].c_str());
  }
  // std::cout << "init_gtest_tests(" << std::endl;
  // for (int i = 0; i < argc; ++i) std::cout << "    " << args[i] << std::endl;
  // std::cout << ") FROM " << __FILE__ << std::endl;
  // testing::InitGoogleMock(&argc, argv);
  testing::InitGoogleTest(&argc, argv);
}

int run_gtest_tests() {
  // std::cout << "RUN_ALL_TESTS FROM " << __FILE__ << std::endl;
  return RUN_ALL_TESTS();
}

#endif
