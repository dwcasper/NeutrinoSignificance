#include "pch.h"
#include "LogPoisson.h"
#include <iostream>

thread_local std::mt19937 LogPoisson::s_twister{ std::random_device{}() };

LogPoisson::LogPoisson(double mu)
	:m_mu{mu}, m_logMu{log(mu)}, m_distribution{m_mu}, m_cache{FactorialCache::GetInstance()}
{ }

//double 
//LogPoisson::EvaluateLogSlow(size_t n) const
//{
//	double logFactorial{ 0 };
//	for (size_t j = 2; j <= n; j++)
//	{
//		logFactorial += log(double(j));
//	}
//	return -m_mu + n * m_logMu - logFactorial;
//}