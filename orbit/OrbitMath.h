#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_

#include<math.h>
#include<iostream>
#include<assert.h>
//#include"constant.h"


//符号函数
template<class T> inline int Sign(const T& InputValue)
{
	if (InputValue > 0)
		return 1;
	else if (InputValue < 0)
		return -1;
	else
		return 0;
}

//求最大值
template <class T>
inline T Max(T x, T y)
{
	return (x > y) ? x : y;
}


//求最小值
template <class T>
inline T Min(T x, T y)
{
	return (x < y) ? x : y;
}


#endif
