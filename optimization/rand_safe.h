#ifndef _RAND_SAFE_H_
#define _RAND_SAFE_H_

#include<random>

//在[min,max)范围内生成随机实数
inline double realRand(const double& min, const double& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}
//在[min,max]范围内生成随机整数
inline int intRand(const int& min, const int& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}

//double realRand(const double& min, const double& max) {
//	static thread_local std::mt19937_64 generator(std::random_device{}());
//	static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
//	return distribution(generator);
//}

#endif
