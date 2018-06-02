#ifndef GTOBASIS_H
#define GTOBASIS_H

class GTOBasis{
private:
	double* m_alphas;
	double* m_coeffs;
	
	int m_L, m_M, m_N;
	int m_count;
	double m_x, m_y, m_z;
	
	void freeMemory();
	void allocateMemory(int n);
	void copyData(const GTOBasis& gb);
public:
	GTOBasis();
	GTOBasis(int n);
	GTOBasis(const GTOBasis& gb);
	GTOBasis& operator=(const GTOBasis& gb);
	~GTOBasis();
public:	
	void setOrbital(double alpha, double coeff, int index);
	void getOrbital(double& alpha, double& coeff, int index) const;
	void setExpoent(int L, int M, int N);
	void getExpoent(int& L, int& M, int& N) const;
	void setPosition(double x, double y, double z);
	void getPosition(double& x, double& y, double& z) const;
	int getCount() const;
	void printBasis() const;
	
};
#endif


