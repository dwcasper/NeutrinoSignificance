#include "pch.h"
#include <iostream>
#include <mutex>
#include "FactorialCache.h"

std::mutex cacheMtx;

FactorialCache::FactorialCache()
	: m_cache { std::pair<size_t, double> {0, 0.0}, std::pair<size_t, double> {1, 0.0} } 
{ 
	std::cout << "Initialized cache with initial cache size = " << m_cache.size() << std::endl;
}

double
FactorialCache::Lookup(size_t n)
{
	// simple case - value already present
	auto itr = m_cache.lower_bound(n);
	if (itr != m_cache.end() && itr->first == n) return itr->second;

	// otherwise, add it
	try
	{
		std::lock_guard<std::mutex> cache_guard(cacheMtx);
		itr = m_cache.lower_bound(n);
		double logFactorial{ 0 };
		if (itr == m_cache.end())  // all entries are less than n
		{
			auto itr = m_cache.rbegin();  // get the largest key entry
			logFactorial = itr->second;
			for (size_t j = itr->first + 1; j <= n; j++)
			{
				logFactorial += log((double)j);
				m_cache[j] = logFactorial;
			}
		}
		else  // either we have entry n, or we have the first one after it
		{
			logFactorial = itr->second;
			for (size_t j = itr->first; j > n; j--)  // if we have the right entry, does nothing
			{
				logFactorial -= log((double)j);
				if (j > 0) m_cache[j - 1] = logFactorial;
			}
		}
		return logFactorial;
	}
	catch (std::exception e) 
	{
		std::cout << "Exception detected and rethrown." << std::endl;
		throw;
	}
}