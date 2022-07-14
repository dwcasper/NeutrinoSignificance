#include <iostream>
#include <numeric>
#include <future>
#include <chrono>
#include <mutex>
#include "SignificanceMC.h"
#include "LogNormal.h"
#include "LogPoisson.h"
#include "boost/math/special_functions/beta.hpp"
#include "boost/math/special_functions/erf.hpp"

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

SignificanceMC::PResult
SignificanceMC::GetPValue(size_t nWorkers, double relativeErrorTarget)
{
	size_t totalUnder{ 0 };
	size_t totalOver{ 0 };
	double p{ 0.0 };
	double sigmaP{ 0.0 };
	size_t n{ 0 };
	double confLow{ 0 };
	double confHigh{ 0 };
	std::cout.precision(std::numeric_limits<double>::max_digits10);
	while (p == 0 || sigmaP == 0 || sigmaP/p > relativeErrorTarget)
	{
		auto result = HandleBatch(nWorkers, totalUnder, totalOver);
		totalOver += result.upper;
		totalUnder += result.lower;
		n += (result.upper + result.lower);
		p = ((double)totalOver) / n;
		sigmaP = sqrt(n * p * (1 - p)) / n;
		confLow = (totalUnder > 0 && totalOver > 0 ? boost::math::ibeta_inv((double)totalOver, (double)totalUnder + 1, kAlphaOneSigma / 2) : p);
		confHigh = (totalUnder > 0 && totalOver > 0 ? boost::math::ibeta_inv((double)totalOver + 1, (double)totalUnder, 1 - kAlphaOneSigma / 2) : p);
		std::cout << "After " << n << " points, Under: " << totalUnder  << " Over: " << totalOver << " pValue: " << p << " +/- " << sigmaP << " {" << confLow << ", " << confHigh << "} significance: " << (p > 0 ? std::sqrt(2.0) * boost::math::erfc_inv(p / 2) : 0) << std::endl;
	}
	return PResult{ p, sigmaP, confLow, confHigh, std::sqrt(2.0)*boost::math::erfc_inv(p/2), std::sqrt(2.0) * boost::math::erfc_inv(confLow / 2), std::sqrt(2.0) * boost::math::erfc_inv(confHigh / 2) };
}

SignificanceMC::BatchResult
SignificanceMC::HandleBatch(size_t nWorkers, size_t under, size_t over) const
{
	static std::vector<std::future<BatchResult> > workers;

	// Make sure the requested number of workers have been started
	while (workers.size() < nWorkers) workers.emplace_back(std::async(std::launch::async, &SignificanceMC::RunBatch, this, under, over));

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

// This is the function run by worker threads
SignificanceMC::BatchResult
SignificanceMC::RunBatch(size_t under, size_t over) const
{
	BatchResult result { 0, 0 };
	while (result.lower + result.upper < m_minBatchSize || result.lower < std::sqrt((double)under) || result.upper < std::sqrt((double)over))
	{
		GridPoint sample = GenerateSample();
		if (Calc_q0(sample) < m_q0Observed)
		{
			result.lower++;
		}
		else
		{
			result.upper++;
		}
	}
	return result;
}

double
SignificanceMC::Calc_q0(GridPoint point) const
{
	if (point[1] * point[1] >= 4 * point[0] * point[2]) return 0; // Negative signal -> no discovery significance
	// From Mathematica solution for best background-only fit
	return 2 * (point[1] * log(2.0 * point[1]) +
		point[0] * log(4.0 * point[0]) -
		(point[1] + 2 * point[0]) * log(point[1] + 2.0 * point[0]) +
		point[2] * log(4.0 * point[2]) +
		(point[0] + point[1] + point[2]) * log((double)(point[0] + point[1] + point[2])) -
		(point[1] + 2 * point[2]) * log(point[1] + 2.0 * point[2]));
}