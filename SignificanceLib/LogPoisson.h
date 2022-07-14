#pragma once
#include <map>
#include <random>
#include "IProposal.h"
#include "FactorialCache.h"
#include "pcg_random.hpp"

class LogPoisson : public IProposal
{
public:
	LogPoisson(double mu);
	~LogPoisson() = default;
	double EvaluateLog(size_t n) const { return -m_mu + n * m_logMu - m_cache.Lookup(n); }
	//double EvaluateLogSlow(size_t n) const;
	double Evaluate(size_t n) const { return exp(EvaluateLog(n)); }
	double LogWeight(size_t n) const { return 0.0; }
	size_t Generate() const { return (size_t) m_distribution(s_pcg); }
private:
	LogPoisson(const LogPoisson& right) = delete;
	LogPoisson& operator=(const LogPoisson& right) = delete;
	double m_mu;
	double m_logMu;
	//thread_local static std::mt19937 s_twister;
	static pcg_extras::seed_seq_from<std::random_device> s_seedSeq;
	thread_local static pcg64_unique s_pcg;
	std::poisson_distribution<> m_distribution;
	FactorialCache& m_cache;
};

