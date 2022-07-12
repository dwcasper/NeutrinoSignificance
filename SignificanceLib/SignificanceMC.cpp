#include "pch.h"
#include <iostream>
#include <numeric>
#include "SignificanceMC.h"
#include "LogNormal.h"
#include "LogPoisson.h"

SignificanceMC::SignificanceMC(const GridPoint& point)
	: m_observation {point}, m_q0Observed{ Calc_q0(point) }
{
	std::cout << "Observed point: " << point[0] << "/" << point[1] << "/" << point[2] << std::endl;
	//m_shat = point[0] - pow(((double)point[1]) / 2, 2) / point[2];
	//m_bhat = pow(((double)point[1]) / 2 + point[2], 2) / point[2];
	//m_phat = ((double)point[1]) / (point[1] + 2 * point[2]);
	//m_globalDist.push_back(new LogPoisson(m_shat + m_phat * m_phat * m_bhat));
	//m_globalDist.push_back(new LogNormal(2 * m_phat * (1 - m_phat) * m_bhat));
	//m_globalDist.push_back(new LogNormal((1 - m_phat) * (1 - m_phat) * m_bhat));

	// Get background-only best parameters and create distributions
	m_bhatBkg = (double)std::accumulate(point.cbegin(), point.cend(), (size_t)0);
	m_phatBkg = ((double)point[1]) / (2 * (point[1] + point[2]));
	m_means = { m_phatBkg * m_phatBkg * m_bhatBkg, 2 * m_phatBkg * (1 - m_phatBkg) * m_bhatBkg, (1 - m_phatBkg) * (1 - m_phatBkg) * m_bhatBkg };
	m_bkgDist.push_back(new LogPoisson(m_phatBkg * m_phatBkg * m_bhatBkg));
	m_bkgDist.push_back(new LogNormal(2 * m_phatBkg * (1 - m_phatBkg) * m_bhatBkg));
	m_bkgDist.push_back(new LogNormal((1 - m_phatBkg) * (1 - m_phatBkg) * m_bhatBkg));
	// Get approximate best data point for background-only
	//size_t mu0 = (size_t)round(m_phatBkg * m_phatBkg * m_bhatBkg);
	//size_t mu1 = (size_t)round(2 * m_phatBkg * (1 - m_phatBkg) * m_bhatBkg);
	//size_t mu2 = (size_t)round((1 - m_phatBkg) * (1 - m_phatBkg) * m_bhatBkg);
}

SignificanceMC::~SignificanceMC()
{
	for (auto d : m_bkgDist)
	{
		if (d != nullptr) delete d;
	}
	m_bkgDist.clear();
}

double
SignificanceMC::GetPValue(size_t N)
{
	double totalUnder{ 0.0 };
	double totalOver{ 0.0 };
	double p{ 0.0 };
	size_t n{ 0 };
	std::cout.precision(std::numeric_limits<double>::max_digits10);
	while (n < N)
	{
		//size_t nBatch{ 0 };
		//double incUnder{ 0.0 };
		//double incOver{ 0.0 };
		//while (nBatch++ < m_minBatchSize || incUnder < m_batchFraction * totalUnder || incOver < m_batchFraction * totalOver)
		//{
		//	Sample sample = GenerateSample();
		//	if (Calc_q0(sample.first) < m_q0Observed)
		//	{
		//		incUnder += sample.second;
		//	}
		//	else
		//	{
		//		incOver += sample.second;
		//	}
		//}
		BatchResult result = RunBatch(totalUnder, totalOver);
		//std::cout << result.count << "/" << result.lower << "/" << result.upper << std::endl;
		totalOver += result.upper;
		totalUnder += result.lower;
		n += result.count;
		p = totalOver / (totalOver + totalUnder);
    	std::cout << "After " << n << " points, TotalProbability = " << std::fixed << (totalOver + totalUnder) / n << " Under: " << totalUnder / n << " Over: " << std::scientific << totalOver / n << " pValue: " << p << " +/- " << sqrt(n * p * (1 - p)) / n << std::endl;
	}
	return p;
}


SignificanceMC::BatchResult
SignificanceMC::RunBatch(double under, double over) const
{
	BatchResult result { 0, 0.0, 0.0 };
	while (result.count < m_minBatchSize || result.lower < m_batchFraction * under || result.upper < m_batchFraction * over)
	{
		Sample sample = GenerateSample();
		if (Calc_q0(sample.first) < m_q0Observed)
		{
			result.lower += sample.second;
		}
		else
		{
			result.upper += sample.second;
		}
		result.count++;
	}
	//std::cout << result.count << ":" << result.lower << ":" << result.upper << std::endl;

	return result;
}

Sample&&
SignificanceMC::GenerateSample() const
{
	std::vector<size_t> points;
	double logWeight{ 0 };
	for (IProposal* p : m_bkgDist)
	{
		size_t value = p->Generate();
		points.push_back(value);
		logWeight += p->LogWeight(value);
	}
	return std::move(std::pair<GridPoint, double> {GridPoint {points[0], points[1], points[2]}, exp(logWeight)});
}

double
SignificanceMC::Calc_q0(GridPoint point) const
{
	if (point[1] * point[1] >= 4 * point[0] * point[2]) return 0; // Negative signal -> no discovery significance
	return 2 * (point[1] * log(2.0 * point[1]) +
		point[0] * log(4.0 * point[0]) -
		(point[1] + 2 * point[0]) * log(point[1] + 2.0 * point[0]) +
		point[2] * log(4.0 * point[2]) +
		(point[0] + point[1] + point[2]) * log((double)(point[0] + point[1] + point[2])) -
		(point[1] + 2 * point[2]) * log(point[1] + 2.0 * point[2]));
}