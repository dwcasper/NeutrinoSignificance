#pragma once
#include <array>
#include <vector>
using GridPoint = std::array<size_t, 3>;
using Sample = std::pair<GridPoint, double>;
class LogNormal;
class IProposal;
class SignificanceMC
{
public:
	SignificanceMC(const GridPoint& observation);
	~SignificanceMC();
	void SetBatchFraction(double fraction) { m_batchFraction = fraction; }
	double GetBatchFraction() const { return m_batchFraction; }
	void SetMinBatchSize(size_t size) { m_minBatchSize = size; }
	size_t GetMinBatchSize() const { return m_minBatchSize; }
	double GetPValue(size_t N);

private:
	SignificanceMC(const SignificanceMC& rhs) = delete;
	SignificanceMC(const SignificanceMC&& rhs) = delete;
	SignificanceMC& operator=(const SignificanceMC& rhs) = delete;
	SignificanceMC& operator=(const SignificanceMC&& rhs) = delete;

	struct BatchResult
	{
	public:
		size_t count;
		double lower;
		double upper;
	};

	BatchResult HandleBatch(double lower, double upper, bool relaunch) const;
	BatchResult RunBatch(double lower, double upper) const;

	Sample&& GenerateSample() const;
	double Calc_q0(const GridPoint grid) const;

	const GridPoint m_observation;
	//double m_shat;
	//double m_phat;
	//double m_bhat;
	double m_phatBkg;
	double m_bhatBkg;
	std::array<double, 3> m_means;

	//std::vector<IProposal*> m_globalDist;
	std::vector<IProposal*> m_bkgDist;

	double m_q0Observed;
	double m_batchFraction{ 0.01 }; // Fraction of total integral to accumulate before adding to get new total
	size_t m_minBatchSize{ 100000 }; // Minimum number of samples in batch, regardless of fraction
	const size_t k_maxWorkers = 20;
};

