#ifndef INCLUDED_scheme_util_timer_HH
#define INCLUDED_scheme_util_timer_HH

#include <chrono>

namespace scheme {
namespace util {

template< class Clock = std::chrono::system_clock >
class Timer {
	std::chrono::time_point<Clock> start;
public:
	Timer(){
		start = Clock::now();
	}
	double elapsed() const {
		std::chrono::duration<double> elapsed_seconds = Clock::now()-start;
		return elapsed_seconds.count();
	}
	double elapsed_nano() const {
		std::chrono::duration<double,std::nano> elapsed_seconds = Clock::now()-start;
		return elapsed_seconds.count();
	}
};

}
}

#endif

