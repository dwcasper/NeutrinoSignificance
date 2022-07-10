#include "pch.h"
#include "Significance.h"
#include "LogPoisson.h"
#include <stdexcept>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ios>
#include <exception>

Significance::Significance(double tail, const GridPoint& observed)
	: m_observed {observed}, m_q0Observed(observed)
{ 
	// Get global best parameters and create distributions
	std::array<size_t, 3> point{ (observed & grid0) >> shift0, (observed & grid1) >> shift1, (observed & grid2) >> shift2 };
	std::cout << point[0] << "/" << point[1] << "/" << point[2] << std::endl;
	m_shat = point[0] - pow(((double) point[1]) / 2, 2) / point[2];
	m_bhat = pow(((double) point[1]) / 2 + point[2], 2) / point[2];
	m_phat = ((double) point[1]) / (point[1] + 2 * point[2]);
	m_globalDist.push_back(new LogPoisson(m_shat + m_phat * m_phat * m_bhat));
	m_globalDist.push_back(new LogPoisson(2 * m_phat * (1 - m_phat) * m_bhat));
	m_globalDist.push_back(new LogPoisson((1 - m_phat)* (1 - m_phat) * m_bhat));

	// Get background-only best parameters and create distributions
	m_bhat_bkg = (double) std::accumulate(point.cbegin(), point.cend(), (size_t)0);
	m_phat_bkg = ((double) point[1]) / (2 * (point[1] + point[2]));
	m_means = { m_phat_bkg * m_phat_bkg * m_bhat_bkg, 2 * m_phat_bkg * (1 - m_phat_bkg) * m_bhat_bkg, (1 - m_phat_bkg) * (1 - m_phat_bkg) * m_bhat_bkg };
	m_bkgDist.push_back(new LogPoisson(m_phat_bkg * m_phat_bkg * m_bhat_bkg));
	m_bkgDist.push_back(new LogPoisson(2 * m_phat_bkg * (1 - m_phat_bkg) * m_bhat_bkg));
	m_bkgDist.push_back(new LogPoisson((1 - m_phat_bkg) * (1 - m_phat_bkg) * m_bhat_bkg));
	// Get approximate best data point for background-only
	size_t mu0 = (size_t) round(m_phat_bkg * m_phat_bkg * m_bhat_bkg);
	size_t mu1 = (size_t) round(2 * m_phat_bkg * (1 - m_phat_bkg) * m_bhat_bkg);
	size_t mu2 = (size_t) round((1 - m_phat_bkg) * (1 - m_phat_bkg) * m_bhat_bkg);
	if (mu0 > (grid0 >> shift0) || mu1 > (grid1 >> shift1) || mu2 > (grid2 >> shift2))
		throw std::out_of_range("starting point out of range");
	GridPoint start { (mu0 << shift0) |
					  (mu1 << shift1) |
					  (mu2 << shift2) };
	std::cout << round(m_phat_bkg * m_phat_bkg * m_bhat_bkg) << "/" << round(2 * m_phat_bkg * (1 - m_phat_bkg) * m_bhat_bkg) << "/" << round((1 - m_phat_bkg) * (1 - m_phat_bkg) * m_bhat_bkg) << std::endl;
	std::cout << std::hex << start << std::dec << std::endl;
	m_tail = Integrate(tail, start);
}

Significance::~Significance()
{}

double 
Significance::Evaluate()
{
	if (m_tail >= 1.0)
	{
		throw std::logic_error("Attempt to evaluate uninitialized Significance.");
	}
	return 0.0;
}

double 
Significance::Integrate(double tail, GridPoint start)
{
	double safety = 0.95;
	size_t index = AddToGrid(start);
	double probabilitySum = Include(index);
	std::cout << std::fixed;
	std::cout.precision(14);
	while (1 - probabilitySum > safety * tail)
	{
		if (m_in.size() % 100000 == 1)
			std::cout << "Points in sum: " << m_in.size() << " Probability integral: " << probabilitySum /* << " or " << std::accumulate(m_in.rbegin(), m_in.rend(), 0.0, [this](double val, size_t index) { return val + m_probability[index]; }) */ <<
			" Total points: " << m_grid.size() << " Surface points: " << m_out.size() << std::endl;
		auto itr = m_out.begin();
		probabilitySum += Include(*itr);
	}
	m_out.clear();
	return 1 - std::accumulate(m_in.rbegin(), m_in.rend(), 0.0, [this](double val, size_t index) { return val + m_probability[index]; });
}

size_t
Significance::AddToGrid(GridPoint point)
{
	size_t index = m_grid.size();
	double probability = CalcProbability(point);
	m_probability.push_back(probability);
	//m_q0.push_back(Calc_q0(point));
	m_lookup[point] = index;
	m_grid.push_back(point);
	m_out.insert(index);
	return index;
}


double
Significance::Include(size_t index)
{
	//std::cout << "Including: " << m_grid[index][0] << "/" << m_grid[index][1] << "/" << m_grid[index][2] << std::endl;
	m_in.insert(index);
	m_out.erase(index);
	FindNeighbors(index);

	return m_probability[index];
}

std::vector<size_t>&&
Significance::FindNeighbors(size_t index)
{
	std::vector<size_t> neighbors;
	size_t grid = m_grid[index];
	std::array<size_t, 3> center{ (grid & grid0) >> shift0, (grid & grid1) >> shift1, (grid & grid2) >> shift2};
	if (center[0] + 1 > (grid0 >> shift0) || center[1] + 1 > (grid1 >> shift1) || center[2] + 1 > (grid2 >> shift2))
		throw std::out_of_range("Neighbor point out of range");
	std::map<GridPoint, size_t>::iterator itr;
	GridPoint trial1{ ((center[0] - 1) << shift0) | (center[1] << shift1) | (center[2] << shift2) };
	GridPoint trial2{ ((center[0] + 1) << shift0) | (center[1] << shift1) | (center[2] << shift2) };
	GridPoint trial3{ (center[0] << shift0) | ((center[1] - 1) << shift1) | (center[2] << shift2) };
	GridPoint trial4{ (center[0] << shift0) | ((center[1] + 1) << shift1) | (center[2] << shift2) };
	GridPoint trial5{ (center[0] << shift0) | (center[1] << shift1) | ((center[2] - 1) << shift2) };
	GridPoint trial6{ (center[0] << shift0) | (center[1] << shift1) | ((center[2] + 1) << shift2) };

	if (center[0] > 0 && (itr = m_lookup.find(trial1)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial1)));
	}
	if ((itr = m_lookup.find(trial2)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial2)));
	}
	if (center[1] > 0 && (itr = m_lookup.find(trial3)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial3)));
	}
	if ((itr = m_lookup.find(trial4)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial4)));
	}
	if (center[2] > 0 && (itr = m_lookup.find(trial5)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial5)));
	}
	if ((itr = m_lookup.find(trial6)) == m_lookup.end())
	{
		neighbors.push_back(AddToGrid(std::move(trial6)));
	}

	return std::move(neighbors);
}


double 
Significance::CalcProbability(GridPoint grid) const
{
	return exp(m_bkgDist[0]->EvaluateLog((grid & grid0) >> shift0) + m_bkgDist[1]->EvaluateLog((grid & grid1) >> shift1) + m_bkgDist[2]->EvaluateLog((grid & grid2) >> shift2));
}

double
Significance::Calc_q0(const GridPoint grid) const
{
	std::array<size_t, 3> point { (grid & grid0) >> shift0, (grid & grid1) >> shift1, (grid & grid2) >> shift2 };
	if (point[1] * point[1] >= 4 * point[0] * point[2]) return 0; // Negative signal -> no discovery significance
	return 2 * (point[1] * log(2.0 * point[1]) +
		point[0] * log(4.0 * point[0]) -
		(point[1] + 2 * point[0]) * log(point[1] + 2.0 * point[0]) +
		point[2] * log(4.0 * point[2]) +
		(point[0] + point[1] + point[2]) * log((double)(point[0] + point[1] + point[2])) -
		(point[1] + 2 * point[2]) * log(point[1] + 2.0 * point[2]));
}

