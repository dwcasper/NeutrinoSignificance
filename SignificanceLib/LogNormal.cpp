#include "pch.h"
#include "LogNormal.h"
#include "LogPoisson.h"
#include <cmath>

thread_local std::mt19937 LogNormal::s_twister{ std::random_device{}()};


LogNormal::LogNormal(double mu)
	: m_mu{ mu }, m_sigma{ std::sqrt(mu) }, m_poisson{ new LogPoisson(m_mu) }, m_distribution{m_mu, m_sigma}
{ }

double 
LogNormal::Evaluate(size_t n) const
{
	double low = n - 0.5;
	double high = n + 0.5;
	return std::erf((high - m_mu) / (std::sqrt(2) * m_sigma))/2 - std::erf((low - m_mu) / (std::sqrt(2) * m_sigma))/2;
}

double 
LogNormal::LogWeight(size_t n) const
{
	return m_poisson->EvaluateLog(n) - EvaluateLog(n);
}