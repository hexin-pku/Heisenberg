#ifndef Hsbg_Spfunc_H
#define Hsbg_Spfunc_H

#include <cstdlib>
#include <iostream>
#include <cmath>
/*
    cmath provide special function since c++17
*/


unsigned int fact(unsigned int n)
{
    if(n==0)
    {
        return 1;
    }else
    {
        return n*fact(n-1);
    }
}

unsigned int dfact(unsigned int n)
{
    if(n==0 || n==1)
    {
        return 1;
    }else
    {
        return n*dfact(n-2);
    }
}

unsigned int comb(unsigned int n, unsigned int k)
{
    if(k==0 || k>= n)
    {
        return 1;
    }else
    {
        return comb(n-1,k) + comb(n-1,k-1);
    }
}

unsigned int comb2(unsigned int n, unsigned int k)
{
    // used k?
    if( k>n/2) k=n-k;
    if( k<=0 )
    {
        return 1;
    }else
    {
        return ( comb(n,k-1) * (n-k+1) )/(k);
    }
}

double legendre(unsigned int l, double x)
{
    if(l==0)
    {
        return 1;
    }else if(l==1)
    {
        return x;
    }else
    {
        return ( (2*l-1)*x*legendre(l-1,x) - (l-1)*legendre(l-2,x) )/double(l);
    }
}


// NOTE: only recieve positive m
double assoc_legendre(unsigned int l, unsigned int m, double x)
{
    if(m==0)
    {
        return legendre(l,x);
    }else if(l==m)
    {
        return - (2*l-1) * sqrt(1-x*x) * assoc_legendre(l-1,l-1,x);
    }else
    {
        return ( (2*l-1)*x*assoc_legendre(l-1,m,x) - (l+m-1)*assoc_legendre(l-2,m,x) )/ double(l-m);
    }
}

double hermite(unsigned int n, double x)
{
    if(n==0)
    {
        return 1;
    }else if(n==1)
    {
        return x;
    }else
    {
        return 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x);
    }
}

#endif






