#ifndef _VECMAT_H_
#define _VECMAT_H_
#include <iostream>
#include <math.h>
#include <assert.h>
//using namespace std;
//释放内存
template<class T> inline void V_Dele(T* B)
{
	delete[] B;
	B=NULL;
}

//将向量A的值赋给B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

//将三个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

//将六个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//向量B与A每个元素都相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++)	if(B[I_]!=A[I_])return false;
	return true;
}

//B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}

//向量C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

//向量C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}

//向量C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

//向量C[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B;
}

//向量C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

//向量C[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B*A[I_];
}

//向量C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

//向量C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//求内积
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

//求外积C=AXB,不能用V_Cross(B,B,A)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

//求A[i]=|B[i]|
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			A[I_]=B[I_];
		else
			A[I_]=-B[I_];
	}
}

//求1-范数
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			result+=B[I_];
		else
			result-=B[I_];
	}
	return result;
}

//求2-范数
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

//求无穷-范数
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}

//求最大元素
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]>result) result=B[I_];
	return result;
}

//求最大元素
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

//求最小元素
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

//求最小元素
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

//find max abs place 
template<class T> T* MaxAbs(T *A, int N, int incx = 1) {
	if (N < 1)
		return NULL;
	T tmp = abs(A[0]);
	int index = 0;
	for (int i = 0; i < N; i++) {
		if (abs(A[i*incx]) > tmp) {
			tmp = abs(A[i*incx]);
			index = i;
		}
	}
	return A + index*incx;
}

//Swap of two vectors
template<class T> void Swap(T *X, T *Y, int N, int incx = 1, int incy = 1) {
	if (N < 1)
		return;
	//if equal to 1 ,use faster way??
	if (incx == 1 && incy == 1) {
		for (int i = 0; i < N; i++) {
			T tmp = X[i];
			X[i] = Y[i];
			Y[i] = tmp;
		}
		return;
	}
	//otherwise use ordinary way
	int ix = 0, iy = 0;
	for (int i = 0; i < N; i++) {
		T tmp = X[ix];
		X[ix] = Y[iy];
		Y[iy] = tmp;
		ix += incx;
		iy += incy;
	}
	return;
}

//Scale by a num
template<class T> void Scale(T *A, T b, int N, int incx = 1) {
	if (0 == b) {
		int indx = 0;
		for (int i = 0; i < N; i++) {
			A[indx] = 0;
			indx += incx;
		}
	}
	if (incx == 1) {
		for (int i = 0; i < N; i++)
			A[i] *= b;
		return;
	}
	//else use stupid way
	int ix = 0;
	for (int i = 0; i < N; i++) {
		A[ix] *= b;
		ix += incx;
	}
}

//y = ax + y
//basic operation of matrix col or row
template<class T> void aXpY(T a, T *X, T *Y, int N, int incx = 1, int incy = 1) {
	if (N < 1)
		return;
	if (incx == 1 && incy == 1) {
		for (int i = 0; i < N; i++)
			Y[i] += a*X[i];
		return;
	}
	int ix = 0, iy = 0;
	for (int i = 0; i < N; i++) {
		Y[iy] += a*X[ix];
		ix += incx;
		iy += incy;
	}
	return;
}

//Copy from X to Y
template<class T> void Copy(T *Y, T *X, int N, int incy = 1, int incx = 1) {
	if (incx == 1 && incy == 1) {
		for (int i = 0; i < N; i++)
			Y[i] = X[i];
		return;
	}
	//other inc
	int ix = 0, iy = 0;
	for (int i = 0; i < N; i++) {
		Y[iy] = X[ix];
		ix += incx;
		iy += incy;
	}
}

//输入流函数
template<class T> inline void V_Input(std::istream  &input, T* Vec, int N)
{
    for(int I_=0;I_<N;I_++)	input>>Vec[I_];
}

//输出流函数
template<class T> inline void V_Output(std::ostream &output, const T* Vec, int N)
{
    for(int I_=0;I_<N;I_++) output<<Vec[I_]<<std::endl;   
}

/**********************************************************************************/
//矩阵，以一维数组表示

//释放内存，与矢量相同

//将九个值依次赋给B
template<class T> inline void M_Copy(T* B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	T temp[9]={a11, a12, a13, a21, a22, a23, a31, a32, a33};
	for(int I_=0;I_<9;I_++)
		B[I_]=temp[I_];
}

template<class T> inline void M_Copy(T* B, T angle, int axis)
{
	assert(axis==1||axis==2||axis==3);
	for(int I_=0;I_<9;I_++) B[I_]=0;
	if(axis==1)
	{
		B[0]=1.0;
		B[4]=B[8]=cos(angle);
		B[5]=sin(angle);
		B[7]=-B[5];
	}
	if(axis==2)
	{
		B[4]=1.0;
		B[0]=B[8]=cos(angle);
		B[6]=sin(angle);
		B[2]=-B[6];
	}
	if(axis==3)
	{
		B[8]=1.0;
		B[0]=B[4]=cos(angle);
		B[1]=sin(angle);
		B[3]=-B[1];
	}
}

//C[i][j]=A[i][k]*B[k][j],C:NXM,A:NXK,B:KXM,不能M_Multi(C,C,B...
template<class T> inline void M_Multi(T* C, const T* A, const T* B, int N, int M, int K)
{
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		C[s]=0;
		for(int i=0;i<K;i++) C[s]+=A[I_*K+i]*B[i*M+J_];
	}
}

//B[i][j]=A[j][i],B:NXM,A:MXN
template<class T> inline void M_Tranpose(T* B, const T* A, int N, int M)
{
	int index=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[index++]=A[J_*N+I_];
}

//求最大元素
template<class T> inline T M_Max(int& row, int& col, const T* B, int N, int M)
{
	T maximal=B[0];
	row=0;
	col=0;
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		if(B[s]>maximal) {row=I_; col=J_;maximal=B[s];}
	}
	return maximal;
}

//求最小元素
template<class T> inline T M_Min(int& row, int& col, const T* B, int N, int M)
{
	T minimal=B[0];
	row=0;
	col=0;
	int s=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		s=I_*M+J_;
		if(B[s]<minimal) {row=I_; col=J_;minimal=B[s];}
	}
	return minimal;
}

//输出流函数
template<class T> inline void M_Output(const T* Vec, int N, int M)//ostream &output, 
{
    for(int I_=0;I_<N;I_++)
	{
		//for(int J_=0;J_<M;J_++)	output<<setprecision(8)<<Vec[I_*M+J_]<<",";
		//output<<endl;
		for(int J_=0;J_<M;J_++)	printf("%.6f%s",Vec[I_*M+J_],",");
		printf("\n");
	}
}

//求逆
void M_Inverse(double* B, const double* A, int N, double* WA);

#endif