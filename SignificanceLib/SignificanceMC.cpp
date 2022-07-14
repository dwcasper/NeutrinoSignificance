#include "pch.h"
#include <iostream>
#include <numeric>
#include <future>
#include <chrono>
#include <mutex>
#include "SignificanceMC.h"
#include "LogNormal.h"
#include "LogPoisson.h"
#include "boost/math/special_functions/beta.hpp"

SignificanceMC::SignificanceMC(const GridPoint& point)
	: m_observation {point}, m_q0Observed{ Calc_q0(point) }, kAlphaOneSigma { 1.0 - std::erf(1.0/std::sqrt(2.0)) }
{
	std::cout << "Observed point: " << point[0] << "/" << point[1] << "/" << point[2] << std::endl;

	// Get background-only best parameters and create distributions
	m_bhatBkg = (double)std::accumulate(point.cbegin(), point.cend(), (size_t)0);
	m_phatBkg = ((double)point[1]) / (2 * (point[1] + point[2]));
	m_means = { m_phatBkg * m_phatBkg * m_bhatBkg, 2 * m_phatBkg * (1 - m_phatBkg) * m_bhatBkg, (1 - m_phatBkg) * (1 - m_phatBkg) * m_bhatBkg };
	m_bkgDist.push_back(new LogPoisson(m_phatBkg * m_phatBkg * m_bhatBkg));
	m_bkgDist.push_back(new LogPoisson(2 * m_phatBkg * (1 - m_phatBkg) * m_bhatBkg));
	m_bkgDist.push_back(new LogPoisson((1 - m_phatBkg) * (1 - m_phatBkg) * m_bhatBkg));
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
SignificanceMC::GetPValue()
{
	double totalUnder{ 0.0 };
	double totalOver{ 0.0 };
	double p{ 0.0 };
	size_t n{ 0 };
	double confLow{ 0 };
	double confHigh{ 0 };
	std::cout.precision(std::numeric_limits<double>::max_digits10);
	while (p == 0 || (confHigh-confLow)/(2*p) > kRelativeError)
	{
		auto result = HandleBatch(totalUnder, totalOver);
		totalOver += result.upper;
		totalUnder += result.lower;
		n += result.count;
		p = totalOver / (totalOver + totalUnder);
		confLow = (totalUnder > 0 && totalOver > 0 ? boost::math::ibeta_inv(totalOver, totalUnder + 1, kAlphaOneSigma / 2) : p);
		confHigh = (totalUnder > 0 && totalOver > 0 ? boost::math::ibeta_inv(totalOver + 1, totalUnder, 1 - kAlphaOneSigma / 2) : p);
		std::cout << "After " << n << " points, TotalProbability = " << std::fixed << (totalOver + totalUnder) / n << " Under: " << totalUnder / n << " Over: " << std::scientific << totalOver / n << " pValue: " << p << " +/- " << sqrt(n * p * (1 - p)) / n << " {" << confLow << ", " << confHigh << "}" << std::endl;
	}
	return p;
}

SignificanceMC::BatchResult
SignificanceMC::HandleBatch(double under, double over) const
{
	static std::vector<std::future<BatchResult> > workers;

	// Make sure the requested number of workers have been started
	while (workers.size() < k_maxWorkers) workers.emplace_back(std::async(std::launch::async, &SignificanceMC::RunBatch, this, under, over));

	// Wait until a thread finishes, and return its results
	while (true)
	{
		for (auto itr = workers.begin(); itr != workers.end(); ++itr)
		{
			if (itr->wait_for(std::chrono::milliseconds{ 500 }) == std::future_status::ready)
			{
				auto result = itr->get();
				workers.erase(std::remove_if(itr, workers.end(), [](std::future<BatchResult> const& f) {return !f.valid(); }), workers.end());
				return result;
			}
		}
	}


}

SignificanceMC::BatchResult
SignificanceMC::RunBatch(double under, double over) const
{
	BatchResult result { 0, 0.0, 0.0 };
	while (result.count < m_minBatchSize || result.lower < std::sqrt(under) || result.upper < std::sqrt(over))
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

	return result;
}

Sample
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
	return std::pair<GridPoint, double> {GridPoint {points[0], points[1], points[2]}, exp(logWeight)};
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