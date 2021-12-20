// copied from https://gist.github.com/mcleary/b0bf4fa88830ff7c882d

#ifndef TIMER_HPP
#define TIMER_HPP
#include <iostream>
#include <chrono>
#include <ctime>
#include <cmath>

class Timer
{
public:
	Timer() : m_StartTime(std::chrono::steady_clock::now()) {

	}

    void start()
    {
        m_StartTime = std::chrono::steady_clock::now();
        m_bRunning = true;
    }

    void stop()
    {
        m_EndTime = std::chrono::steady_clock::now();
        m_bRunning = false;
    }

	template<typename U>
	double elapsed()
	{
		std::chrono::time_point<std::chrono::steady_clock> endTime;

		if(m_bRunning)
		{
			endTime = std::chrono::steady_clock::now();
		}
		else
		{
			endTime = m_EndTime;
		}

		return std::chrono::duration_cast<U>(endTime - m_StartTime).count();
	}

private:
    std::chrono::time_point<std::chrono::steady_clock> m_StartTime;
    std::chrono::time_point<std::chrono::steady_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

#endif