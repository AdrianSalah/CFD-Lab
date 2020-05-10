#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>


class Timer
{
private:
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1>>;

	std::chrono::time_point<clock_t> m_begin;

public:
	Timer() : m_begin{ clock_t::now() } {}
	void reset() {
		m_begin = clock_t::now();
	}
	double elapsed() const {
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_begin).count();
	}

	void printTimer() {
		std::cout << "----------------------------------------------------" << std::endl;
		std::cout << "Elapsed time, sec = " << elapsed() << "\n" << std::endl;
	}
};


void solutionProgress(double const& currentTime, double t_end) {
	static int num = 0;
	if (num == 0)
	{
		std::cout << "Solution progress:" << std::endl;
		std::cout << "\t" << round(currentTime / t_end * 100) << "%" << std::endl;
	}
	else
		std::cout << "\t" << round(currentTime / t_end * 100) << "%" << std::endl;

	++num;
}

#endif // !TIMER_H
