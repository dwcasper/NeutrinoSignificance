#pragma once
#include <map>
class LogPoisson
{
public:
	LogPoisson(double mu);
	~LogPoisson() = default;
	double EvaluateLog(size_t n);
	double EvaluateLogSlow(size_t n);
	double Evaluate(size_t n) { return exp(EvaluateLog(n)); }
private:
	LogPoisson(const LogPoisson& right) = delete;
	LogPoisson& operator=(const LogPoisson& right) = delete;
	double m_mu;
	double m_logMu;
	static std::map<size_t, double> s_cache;
};

