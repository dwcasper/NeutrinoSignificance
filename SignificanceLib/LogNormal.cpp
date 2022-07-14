#include "pch.h"
#include "LogNormal.h"
#include "LogPoisson.h"
#include <cmath>
#include <iostream>

thread_local std::mt19937 LogNormal::s_twister{ std::random_device{}()};


LogNormal::LogNormal(double mu)
	: m_mu{ mu }, m_sigma{ std::sqrt(mu) }, m_poisson{ new LogPoisson(m_mu) }, m_distribution{m_mu, m_sigma}
{
	m_normalize = 1 - std::erfc((-0.5 - m_mu) / (std::sqrt(2) * m_sigma))/2;
	std::cout.precision(std::numeric_limits<double>::max_digits10);
	std::cout << std::fixed << "For mu = " << m_mu << ", normalization = " << m_normalize << std::endl;
}

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
	return (n >= 0 ? m_poisson->EvaluateLog(n) - EvaluateLog(n) : std::numeric_limits<double>::lowest());
}

size_t
LogNormal::Generate() const
{
	size_t result {0};
	do
	{
		result = (size_t)std::round(m_distribution(s_twister));
	} while (result < 0);

	return result;
}