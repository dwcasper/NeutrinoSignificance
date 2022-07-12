#pragma once
class IProposal
{
public:
	virtual double Evaluate(size_t n) const = 0;
	virtual double EvaluateLog(size_t n) const = 0;
	virtual double LogWeight(size_t n) const = 0;
	virtual size_t Generate() const = 0;
};

