#include "pch.h"
#include "LogPoisson.h"
#include <iostream>

std::map<size_t, double> LogPoisson::s_cache = { std::pair<size_t, double> {0, 0.0}, std::pair<size_t, double> {1, 0.0} };

LogPoisson::LogPoisson(double mu)
	:m_mu(mu), m_logMu(log(mu))
{ }

double 
LogPoisson::EvaluateLog(size_t n)
{
	double logFactorial{ 0 };
	auto itr = s_cache.lower_bound(n);
	if (itr == s_cache.end())  // all entries are less than n
	{
		auto itr = s_cache.rbegin();  // get the largest key entry
		logFactorial = itr->second;
		for (size_t j = itr->first + 1; j <= n; j++)
		{
			logFactorial += log((double)j);
			s_cache[j] = logFactorial;
		}
	}
	else  // either we have entry n, or we have the first one after it
	{
		logFactorial = itr->second;
		for (size_t j = itr->first; j > n; j--)  // if we have the right entry, does nothing
		{
			logFactorial -= log((double)j);
			if (j > 0) s_cache[j - 1] = logFactorial;
		}
	}
	return -m_mu + n * m_logMu - logFactorial;
}

double 
LogPoisson::EvaluateLogSlow(size_t n)
{
	double logFactorial{ 0 };
	for (size_t j = 2; j <= n; j++)
	{
		logFactorial += log(double(j));
	}
	return -m_mu + n * m_logMu - logFactorial;
}