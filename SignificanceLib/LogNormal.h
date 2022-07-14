#pragma once
#include <random>
#include "IProposal.h"

class LogPoisson;

class LogNormal : public IProposal
{
public:
	LogNormal(double mu);
	~LogNormal() = default;
	double EvaluateLog(size_t n) const { return std::log(Evaluate(n)); };
	double Evaluate(size_t n) const;
	double LogWeight(size_t n) const;
	size_t Generate() const;
private:
	LogNormal(const LogNormal& rhs) = delete;
	LogNormal& operator=(const LogNormal& rhs) = delete;

	double m_mu;
	double m_sigma;
	double m_normalize;
	LogPoisson* m_poisson;
	thread_local static std::mt19937 s_twister;
	mutable std::normal_distribution<> m_distribution;
};

