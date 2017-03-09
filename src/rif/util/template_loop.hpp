/// @file   util/template_math.hpp
/// @brief  a place for some template loops to live
/// @author will sheffler

#ifndef INCLUDED_util_template_loop_HH
#define INCLUDED_util_template_loop_HH

#include <boost/static_assert.hpp>
#include <cstddef>

namespace rif {
namespace util {

template <int DIM>
struct NESTED_FOR {
  // this is to make sure a compile errer happens if DIM not implemented
  BOOST_STATIC_ASSERT_MSG(DIM == 123456789,
                          "NESTED_FOR not implemented for requested DIM");
};

template <>
struct NESTED_FOR<1> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
      f(i);
    }
  }
};

template <>
struct NESTED_FOR<2> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
      for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
        f(i);
      }
    }
  }
};

template <>
struct NESTED_FOR<3> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
      for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
        for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
          f(i);
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<4> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
      for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
        for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
          for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
            f(i);
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<5> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
      for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
        for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
          for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
            for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
              f(i);
            }
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<6> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[5] = lb[5]; i[5] <= ub[5]; ++i[5]) {
      for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
        for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
          for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
            for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
              for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
                f(i);
              }
            }
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<7> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[6] = lb[6]; i[6] <= ub[6]; ++i[6]) {
      for (i[5] = lb[5]; i[5] <= ub[5]; ++i[5]) {
        for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
          for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
            for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
              for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
                for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
                  f(i);
                }
              }
            }
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<8> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[7] = lb[7]; i[7] <= ub[7]; ++i[7]) {
      for (i[6] = lb[6]; i[6] <= ub[6]; ++i[6]) {
        for (i[5] = lb[5]; i[5] <= ub[5]; ++i[5]) {
          for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
            for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
              for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
                for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
                  for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
                    f(i);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<9> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[8] = lb[8]; i[8] <= ub[8]; ++i[8]) {
      for (i[7] = lb[7]; i[7] <= ub[7]; ++i[7]) {
        for (i[6] = lb[6]; i[6] <= ub[6]; ++i[6]) {
          for (i[5] = lb[5]; i[5] <= ub[5]; ++i[5]) {
            for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
              for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
                for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
                  for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
                    for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
                      f(i);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

template <>
struct NESTED_FOR<10> {
  template <class Functor, class Indices>
  NESTED_FOR(Indices const &lb, Indices const &ub, Functor &f) {
    Indices i;
    for (i[9] = lb[9]; i[9] <= ub[9]; ++i[9]) {
      for (i[8] = lb[8]; i[8] <= ub[8]; ++i[8]) {
        for (i[7] = lb[7]; i[7] <= ub[7]; ++i[7]) {
          for (i[6] = lb[6]; i[6] <= ub[6]; ++i[6]) {
            for (i[5] = lb[5]; i[5] <= ub[5]; ++i[5]) {
              for (i[4] = lb[4]; i[4] <= ub[4]; ++i[4]) {
                for (i[3] = lb[3]; i[3] <= ub[3]; ++i[3]) {
                  for (i[2] = lb[2]; i[2] <= ub[2]; ++i[2]) {
                    for (i[1] = lb[1]; i[1] <= ub[1]; ++i[1]) {
                      for (i[0] = lb[0]; i[0] <= ub[0]; ++i[0]) {
                        f(i);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};
}
}

#endif  // INCLUDED_util_template_math_HH
