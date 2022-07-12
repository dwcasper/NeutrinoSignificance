#pragma once
#include <map>
class FactorialCache
{
public:
	static FactorialCache& GetInstance()
	{
		static FactorialCache instance;
		return instance;
	}

	double Lookup(size_t n);

private:
	FactorialCache();
	~FactorialCache() = default;
	FactorialCache(const FactorialCache& rhs) = delete;
	FactorialCache(const FactorialCache&& rhs) = delete;
	FactorialCache& operator=(const FactorialCache& rhs) = delete;
	FactorialCache& operator=(const FactorialCache&& rhs) = delete;
	static FactorialCache* s_instance;
	std::map<size_t, double> m_cache;
};

