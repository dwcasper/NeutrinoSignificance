#pragma once
#include <vector>
#include <map>
#include <set>
#include <array>
class LogPoisson;

//using GridPoint = std::array<size_t, 3>;
using GridPoint = size_t;
class Significance
{
public:
	Significance(double tail, const GridPoint& observed);
	~Significance();

	double Evaluate();
	GridPoint GetObserved() const { return m_observed; }
	std::array<double, 3> GetMeans()    const { return m_means; }
	double GetTail() const { return m_tail; }
	double Get_q0Observed() const { return m_q0Observed; }
	static const size_t grid2 = 0x00000000FFFFFFFF;
	static const size_t grid1 = 0x000FFFFF00000000;
	static const size_t grid0 = 0xFFF0000000000000;
	static const size_t shift2 = 0;
	static const size_t shift1 = 32;
	static const size_t shift0 = 52;
private:
	struct ByProbability
	{
		ByProbability(const Significance* parent)
		{
			m_parent = parent;
		}
		bool operator() (const size_t& lhs, const size_t& rhs) const
		{
			return m_parent->GetProbability(lhs) > m_parent->GetProbability(rhs);
		}
		const Significance* m_parent{ nullptr };
	};

	double Integrate(double tail, GridPoint start);
	double CalcProbability(const GridPoint grid) const;
	double Calc_q0(const GridPoint grid) const;
	double Include(size_t index);
	size_t AddToGrid(GridPoint point);
	std::vector<size_t>&& FindNeighbors(size_t index);
	double GetProbability(size_t index) const { return m_probability[index]; }

	double m_shat, m_bhat, m_phat;
	double m_bhat_bkg, m_phat_bkg;

	double m_tail {1.0};
	GridPoint m_observed;
	double m_q0Observed;
	std::array<double, 3> m_means;

	std::vector<LogPoisson*> m_globalDist;
	std::vector<LogPoisson*> m_bkgDist;

	std::vector<GridPoint> m_grid;
	std::map<GridPoint, size_t> m_lookup;
	std::vector<double> m_probability;
	std::vector<double> m_q0;
	std::set<size_t, ByProbability> m_done {{}, ByProbability(this)};
	std::set<size_t, ByProbability> m_in   {{}, ByProbability(this)};
	std::set<size_t, ByProbability> m_out  {{}, ByProbability(this)};
};

