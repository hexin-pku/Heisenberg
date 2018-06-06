/*
	this is a Matrix structure of high dimensions!
	this can bind with Eigen's default data type.
	but Eigen can use <Tensor> for high dimensions, so it boud to be abandoned.
*/

#ifndef Hsbg_Others_H
#define Hsbg_Others_H
#include <Eigen/Dense>

class HMatrix
{
	double *M;
	int 	dim;
	int		sz;
	int		ds;
	int*	w;
	
	HMatrix(int my_dim, int my_sz)
	{
		this->dim = my_dim;
		this->sz = my_sz;
		this->ds = 1;
		this->w = new int(this->dim);
		for(int i=0; i<my_dim; i++)
		{
			if(i==0)
			{
				this->w[this->dim-1-i] = 1;
			}
			else
			{
				this->w[this->dim-1-i] = this->w[this->dim - i]*this->sz;
			}
			this->ds *= this->size;
		}
		if(this->ds > 10000000000)
		{
			cerr << "the Matrix is too big" << endl; exit(-1);
		}
		this->M = new double[this->ds];
	}
	
	int set_DS(int my_dim, int my_sz)
	{
		this->dim = my_dim;
		this->sz = my_sz;
		this->ds = 1;
		this->w = new int(this->dim);
		for(int i=0; i<my_dim; i++)
		{
			if(i==0)
			{
				this->w[this->dim-1-i] = 1;
			}
			else
			{
				this->w[this->dim-1-i] = this->w[this->dim - i]*this->sz;
			}
			this->ds *= this->size;
		}
		if(this->ds > 10000000000)
		{
			cerr << "the Matrix is too big" << endl; exit(-1);
		}
		this->M = new double[this->ds];
		return 0;
	}
	
	friend double& operator[](int idxs[])
	{
		int idim=getArrayLen(idxs);
		if(idim!=this->dim) {cerr << "matrix not fix with dim" << endl; exit(-1);}
		int cnt=0;
		for(int i=0; i< this->sz; i++)
		{
			if(idxs[i]>=this->dim) {cerr << "matrix not fix with size" << endl; exit(-1);}
			cnt+=this->w[i]*idxs[i];
		}
		return this->M[cnt];
	}
	
	Matrix2d& operator=(const HMatrix& HM)
	{
		if(HM.dim!=2) {cerr << "matrix interface with eigen should in dim-2" << endl; exit(-1);}
		Matrix2d eigen_M(HM.sz,HM.sz);
		for(int i=0; i<HM.sz; i++)
		{
			for(int j=0; j<HM.sz; j++)
			{
				eigen_M(i,j) = HM[i,j];
			}
		}
	}
}

#endif
