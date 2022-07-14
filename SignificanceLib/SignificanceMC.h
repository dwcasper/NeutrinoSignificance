#pragma once
#include <array>
#include <vector>
#include "IProposal.h"
using GridPoint = std::array<size_t, 3>;
class LogNormal;
class SignificanceMC
{
public:
	SignificanceMC(const GridPoint& observation);
	~SignificanceMC();
	void SetMinBatchSize(size_t size) { m_minBatchSize = size; }
	size_t GetMinBatchSize() const { return m_minBatchSize; }

	struct PResult
	{
	public:
		double p;
		double sigmaP;
		double pLow;
		double pHigh;
		double z;
		double zLow;
		double zHigh;
	};

	PResult GetPValue(size_t nWorkers = 1, double relativeErrorTarget = 0.01);

private:
	SignificanceMC(const SignificanceMC& rhs) = delete;
	SignificanceMC(const SignificanceMC&& rhs) = delete;
	SignificanceMC& operator=(const SignificanceMC& rhs) = delete;
	SignificanceMC& operator=(const SignificanceMC&& rhs) = delete;

	struct BatchResult
	{
	public:
		size_t lower;
		size_t upper;
	};

	BatchResult HandleBatch(size_t nWorkers, size_t lower, size_t upper) const;
	BatchResult RunBatch(size_t lower, size_t upper) const;

	GridPoint GenerateSample() const { return GridPoint{ m_bkgDist[0]->Generate(), m_bkgDist[1]->Generate(), m_bkgDist[2]->Generate() }; }
	double Calc_q0(const GridPoint grid) const;

	const GridPoint m_observation;
	double m_phatBkg;
	double m_bhatBkg;
	std::array<double, 3> m_means;

	std::vector<IProposal*> m_bkgDist;

	double m_q0Observed;
	size_t m_minBatchSize{ 1000000 }; // Minimum number of samples in batch, regardless of fraction
	const double kAlphaOneSigma;
};

