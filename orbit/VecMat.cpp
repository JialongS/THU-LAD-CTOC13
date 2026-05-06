#include "VecMat.h"

//NXM的矩阵A中的行row0与row1元素互换
template<class T> inline void swaprows(T* A, int row0, int row1, int N, int M)
{
    T temp=0;
	int s0=row0*M,s1=row1*M;
	for(int i=0;i<M;i++)
	{
		temp=A[s0];
		A[s0++]=A[s1];
		A[s1++]=temp;
	}
}
//	gjelim 
//A,B:NXN,WA:NX2N
void M_Inverse(double* B, const double* A, int N, double* WA)
{
	int row,col,dindex, N2=2*N;
	int swa=0,sa=0;
    for (row=0; row<N; ++row)
	{
        for (col=0; col<N; ++col)
            WA[swa++]=A[sa++];
        for (col=N; col<N2; ++col)
		{
			if(row==col-N) WA[swa++]=1.0;
			else WA[swa++]=0.0;
		}
    }
    //	perform forward elimination to get WA in row-echelon form
    for (dindex=0; dindex<N; ++dindex)
	{
		int ss=dindex*N2+dindex;
        //	run along diagonal, swapping rows to move zeros in working position 
        //	(along the diagonal) downwards
        if ( (dindex==(N-1)) && (WA[ss]==0)) return; //  no solution
	//	else if (WA[ss]==0) swaprows(WA, dindex, dindex+1,N,N2);
		int maxrow=dindex;
		double maxvalue=fabs(WA[ss]);
		for(int kk=dindex+1;kk<N;kk++)
		{
			if (fabs(WA[kk*N2+dindex])>maxvalue) {maxvalue=fabs(WA[kk*N2+dindex]);maxrow=kk;}
		}
		if (maxrow!=dindex) swaprows(WA, dindex, maxrow,N,N2);
        //	divide working row by value of working position to get a 1 on the
        //	diagonal
        if (WA[ss] == 0.0) return;
		else
		{
            double tempval=WA[ss];
			int index=dindex*N2;
            for (col=0; col<N2; ++col) WA[index++]/=tempval;
        }

        //	eliminate value below working position by subtracting a multiple of 
        //	the current row
        for (row=dindex+1; row<N; ++row)
		{
            double wval=WA[row*N2+dindex];
            for (col=0; col<N2; ++col) WA[row*N2+col]-=wval*WA[dindex*N2+col];
        }
    }

    //	backward substitution steps
    for (dindex=N-1; dindex>=0; --dindex)
	{
        //	eliminate value above working position by subtracting a multiple of 
        //	the current row
        for (row=dindex-1; row>=0; --row)
		{
            double wval=WA[row*N2+dindex];
            for (col=0; col<N2; ++col)
                WA[row*N2+col]-=wval*WA[dindex*N2+col];
        }
    }

    //	assign result to replace B
    for (row=0; row<N; ++row)
	{
        for (col=0; col<N; ++col)
			B[row*N+col]=WA[row*N2+col+N];
    }
}