#include "GTOBasis.h"
#include <string>
#include <iostream>
#include <exception>
#include <iomanip>
void GTOBasis::freeMemory(){
	
	//std::cout << "freeMemory" << std::endl;

	if(m_alphas) delete [] m_alphas;
	if(m_coeffs) delete [] m_coeffs;
	m_alphas = nullptr;
	m_coeffs = nullptr;
	
}

void GTOBasis::allocateMemory(int n){
	//std::cout << "allocateMemory" << std::endl;
	
	m_alphas = new double[n];
	m_coeffs = new double[n];
}

void GTOBasis::copyData(const GTOBasis& gb){
	//std::cout << "copyData" << std::endl;
	m_count = gb.m_count;
	m_L = gb.m_L;
	m_M = gb.m_M;
	m_N = gb.m_N;
	m_x = gb.m_x;
	m_y = gb.m_y;
	m_z = gb.m_z;
	for(int i = 0; i < m_count; ++i){
		m_alphas[i] = gb.m_alphas[i];
		m_coeffs[i] = gb.m_coeffs[i];
	}
}

GTOBasis::GTOBasis(){
	//std::cout << "mo ren" << std::endl;
	m_alphas = nullptr;
	m_coeffs = nullptr;
	m_count = m_L = m_M = m_N = 0;
	m_x = m_y = m_z = 0;
}

GTOBasis::GTOBasis(int n){
	//std::cout << "n" << std::endl;
	m_count = n;
	m_L = m_M = m_N = 0;
	m_x = m_y = m_z = 0;
	allocateMemory(n); 
}

GTOBasis::GTOBasis(const GTOBasis& gb){
	//std::cout << "copy construct" << std::endl;
	allocateMemory(gb.m_count);
	copyData(gb);
}


GTOBasis& GTOBasis::operator=(const GTOBasis& gb){
	//std::cout << "operator=" << std::endl;
	if(&gb == this) return *this;
	if(m_count != gb.m_count){
		freeMemory();
		allocateMemory(gb.m_count);
	}	
	copyData(gb);
}

void GTOBasis::setOrbital(double alpha, double coeff, int index){
	if(index >= 0 && index < m_count){		
		m_alphas[index] = alpha;
		m_coeffs[index] = coeff;
	}
	else{
		std::string errorMsg = "setOrbital: out of range, m_count = " + std::to_string(m_count) + " ,index = " + std::to_string(index); 
	
		throw std::out_of_range(errorMsg);
	}
}
void GTOBasis::getOrbital(double& alpha, double& coeff, int index) const{
	if(index >= 0 && index < m_count){
		alpha = m_alphas[index];
		coeff = m_coeffs[index];
	}
	else{
		std::string errorMsg = "setOrbital: out of range, m_count = " + std::to_string(m_count) + " ,index = " + std::to_string(index); 
	
		throw std::out_of_range(errorMsg);
	}
}
void GTOBasis::setExpoent(int L, int M, int N){
	m_L = L;
	m_M = M;
	m_N = N;
}
void GTOBasis::getExpoent(int& L, int& M, int& N) const{
	L = m_L;
	M = m_M;
	N = m_N;
}
int GTOBasis::getCount() const { return m_count; }

GTOBasis::~GTOBasis(){
	//std::cout << "destroy" <<std::endl;
	freeMemory();
}

void GTOBasis::printBasis() const{
	std::cout<<m_L<<"\t"<<m_M<<"\t"<<m_N<<std::endl;
	std::cout<<m_x<<"\t"<<m_y<<"\t"<<m_z<<std::endl;	
	for(int i = 0; i < m_count; ++i){
		std::cout<<m_alphas[i]<<"   "<<m_coeffs[i]<<"\n";
	}
	std::cout<<std::endl;
}
	
void GTOBasis::setPosition(double x, double y, double z){
	m_x = x;
	m_y = y;
	m_z = z;
}
void GTOBasis::getPosition(double& x, double& y, double& z) const{
	x = m_x;
	y = m_y;
	z = m_z;
}
	
	
	
	
	
	