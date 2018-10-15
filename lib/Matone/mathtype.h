#ifndef MATHTYPE_H
#define MATHTYPE_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

#define PI 3.14159265358979323846
#define EE 2.71828182845904523536

////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DEFINITION OF REAL TYPE

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////

class Real
{
    public:
    double val;
    
    // Constructors and destructor
    Real();
    Real(const double& d);
    Real(const Real& r);
    ~Real();
    
    // Basic operator
    Real operator+(const int& i);
    Real operator-(const int& i);
    Real operator*(const int& i);
    Real operator/(const int& i);
    
    Real operator+(const double& d);
    Real operator-(const double& d);
    Real operator*(const double& d);
    Real operator/(const double& d);
    
    Real operator+(const Real& r);
    Real operator-(const Real& r);
    Real operator-();
    Real operator*(const Real& r);
    Real operator/(const Real& r);
        
    // Type conversion
    operator int() {return (int)this->val;}
    operator double() {return (double)this->val;}
    
    // Assignment
    Real& operator=(const int& i)       {this->val = (double)i; return *this;}
    Real& operator=(const double& d)    {this->val = d; return *this;}
    Real& operator=(const Real& r)      {this->val = r.val; return *this;}
    
    // IOstream
    friend std::ostream& operator<<(std::ostream& output, const Real& r)
    {
        output << std::setiosflags(std::ios::left) << std::setprecision(8);
        output << r.val;
        return output;
    }
    friend std::istream& operator>>(std::istream& input, Real& r)
    {
        input >> r.val;
        return input;
    }
    
    // friend functions
    friend inline Real operator+(const int& i, const Real& r)
    {
        return Real((double)i+r.val);
    }
    friend inline Real operator-(const int& i, const Real& r)
    {
        return Real((double)i-r.val);
    }
    friend inline Real operator*(const int& i, const Real& r)
    {
        return Real((double)i * r.val);
    }
    friend inline Real operator/(const int& i, const Real& r)
    {
        return Real((double)i / r.val);
    }
    friend inline Real operator+(const double& d, const Real& r)
    {
        return Real(d+r.val);
    }
    friend inline Real operator-(const double& d, const Real& r)
    {
        return Real(d-r.val);
    }
    friend inline Real operator*(const double& d, const Real& r)
    {
        return Real(d*r.val);
    }
    friend inline Real operator/(const double& d, const Real& r)
    {
        return Real(d/r.val);
    }
    friend inline Real power(const Real& r, const int& i )
    {
        return Real(pow(r.val, (double)i));
    }
    friend inline Real power(const Real& r, const double& d)
    {
        return Real(pow(r.val, d));
    }
    friend inline Real power(const int& i, const Real& r)
    {
        return pow((double)i, r.val);
    }
    friend inline Real power(const double& d, const Real& r)
    {
        return pow(d, r.val);
    }
    friend inline Real power(const Real& r1, const Real& r2)
    {
        return Real(pow(r1.val, r2.val));
    }
    // Not friend power function
    Real power(const double& d) { return Real(pow(this->val, d )); }
    Real power(const Real& r)   { return Real(pow(this->val, r.val)) ;}
};


Real::Real()
{
    this->val = 1.0;
}

Real::Real(const double& d)
{
    this->val = d;
}

Real::Real(const Real& r)
{
    this->val = r.val;
}

Real::~Real(){};

// basic operator with ints
Real Real::operator+(const int& i)
{
    return Real(this->val + (double)i);
}

Real Real::operator-(const int& i)
{
    return Real(this->val - (double)i);
}

Real Real::operator*(const int& i)
{
    return Real(this->val * (double)i);
}

Real Real::operator/(const int& i)
{
    return Real(this->val / (double)i);
}

// basic operator with doubles
Real Real::operator+(const double& d)
{
    return Real(this->val + d);
}

Real Real::operator-(const double& d)
{
    return Real(this->val - d);
}

Real Real::operator*(const double& d)
{
    return Real(this->val * d);
}

Real Real::operator/(const double& d)
{
    return Real(this->val / d);
}

// basic operator with Real
Real Real::operator+(const Real& r)
{
    return Real(this->val + r.val);
}

Real Real::operator-(const Real& r)
{
    return Real(this->val - r.val);
}

Real Real::operator-()
{
    return Real(-this->val);
}

Real Real::operator*(const Real& r)
{
    return Real(this->val * r.val);
}

Real Real::operator/(const Real& r)
{
    return Real(this->val / r.val);
}



////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DEFINITION OF COMPLEX TYPE

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////

class Complex
{
    public:
    double  real, imag;
    Complex();
    Complex(const double& r, const double& im);
    Complex(const double& r, const double& im, const bool& ispolar);
    Complex(const Real& r, const Real& im);
    Complex(const Real& r, const Real& im, const bool& ispolar);
    Complex(const Complex& c);
    ~Complex();
    
    Complex Conj();
    double  Angle();
    double  const_Angle() const;
    double  Norm();
    double  const_Norm() const;
    double  Norm2();
    
    Complex  operator+(const int& i);
    Complex  operator-(const int& i);
    Complex  operator*(const int& i);
    Complex  operator/(const int& i);
    
    Complex  operator+(const double& d);
    Complex  operator-(const double& d);
    Complex  operator*(const double& d);
    Complex  operator/(const double& d);
    
    Complex  operator+(const Real& r);
    Complex  operator-(const Real& r);
    Complex  operator*(const Real& r);
    Complex  operator/(const Real& r);
    
    Complex  operator+(const Complex& c);
    Complex  operator-(const Complex& c);
    Complex  operator-();
    Complex  operator*(const Complex& c);
    Complex  operator/(const Complex& c);
    
    Complex& operator=(const int& i)
    {
        this->real = (double)i;
        this->imag = 0.0;
        return *this;
    }
    Complex& operator=(const double& d)
    {
        this->real = d;
        this->imag = 0.0;
        return *this;
    }
    Complex& operator=(const Real& r)
    {
        this->real = r.val;
        this->imag = 0.0;
        return *this;
    }
    Complex& operator=(const Complex& c)
    {
        this->real = c.real;
        this->imag = c.imag;
        return *this;
    }
    
    operator double() {return this->real;}
    operator Real()   {return Real(this->real);}

    // IOstream
    friend std::ostream& operator<<(std::ostream& output, const Complex& c)
    {
        output << std::setiosflags(std::ios::left) << std::setprecision(8);
        if(c.imag > 0.0)
        {
            output << c.real << "+" << c.imag << "(im)" ;
        }else{
            output << c.real << "-" << -c.imag << "(im)" ;
        }
        return output;
    }
    friend std::istream& operator<<(std::istream& input, Complex& c)
    {
        input >> c.real >> c.imag;
        return input;
    }
    // friend functions
    friend inline Complex operator+(const int& i, const Complex& c)
    {
        return Complex((double)i+c.real, c.imag);
    }
    friend inline Complex operator-(const int& i, const Complex& c)
    {
        return Complex((double)i-c.real, -c.imag);
    }
    friend inline Complex operator*(const int& i, const Complex& c)
    {
        return Complex(c.real * (double)i, c.imag * (double)i );
    }
    friend inline Complex operator/(const int& i, const Complex& c) 
    {
        double m = (double)i/(c.real*c.real+c.imag*c.imag);
        return Complex(c.real*m, -c.imag*m);
    }
    friend inline Complex operator+(const double& d, const Complex& c)
    {
        return Complex(d+c.real, c.imag);
    }
    friend inline Complex operator-(const double& d, const Complex& c)
    {
        return Complex(d-c.real, -c.imag);
    }
    friend inline Complex operator*(const double& d, const Complex& c)
    {
        return Complex(c.real*d, c.imag*d);
    }
    friend inline Complex operator/(const double& d, const Complex& c)
    {
        double m = d/(c.real*c.real+c.imag*c.imag);
        return Complex(c.real*m, -c.imag*m);
    }
    friend inline Complex operator+(const Real& r, const Complex& c)
    {
        return Complex(r.val+c.real, c.imag);
    }
    friend inline Complex operator-(const Real& r, const Complex& c)
    {
        return Complex(r.val-c.real, -c.imag);
    }
    friend inline Complex operator*(const Real& r, const Complex& c)
    {
        return Complex(r.val*c.real, r.val*c.imag);
    }
    friend inline Complex operator/(const Real& r, const Complex& c) 
    {
        double m = r.val/(c.real*c.real+c.imag*c.imag) ;
        return Complex(c.real*m, -c.imag*m) ;
    }
    friend inline Complex power(const int& i, const Complex& c)
    {
        return (i==0)? Complex(0,0) : Complex(pow((double)i,c.real), log((double)i)*c.imag, true);
    }
    friend inline Complex power(const double& d, const Complex& c)
    {
        return (d==0)? Complex(0,0) : Complex(pow(d,c.real), log(d)*c.imag, true);
    }
    friend inline Complex power(const Real& r, const Complex& c)
    {
        return (r.val==0)? Complex(0,0) : Complex(pow(r.val,c.real), log(r.val)*c.imag, true);
    }
    friend inline Complex power(const Complex& c, const int& i)
    {
        double n,a;
        n = c.const_Norm(); a = c.const_Angle();
        return (c.real==0 && c.imag==0)? Complex(0,0) : Complex(pow(n,(double)i), (double)i*a, true);
    }
    friend inline Complex power(const Complex& c, const double& d)
    {
        double n,a;
        n = c.const_Norm(); a = c.const_Angle();
        return (c.real==0 && c.imag==0)? Complex(0,0) : Complex(pow(n,d), d*a, true);
    }
    friend inline Complex power(const Complex& c, const Real& r)
    {
        double n,a;
        n = c.const_Norm(); a = c.const_Angle();
        return (c.real==0 && c.imag==0)? Complex(0,0) : Complex(pow(n,r.val), r.val*a, true);
    }
    friend inline Complex power(const Complex& c1, const Complex& c2)
    {
        double n,a;
        n = c1.const_Norm(); a = c1.const_Angle();
        return (c1.real==0 && c1.imag==0) ? Complex(0,0) : ( 
        Complex( pow(n,c2.real)*exp(-a*c2.imag),  a*c2.real+log(n)*c2.imag,  true) );
    }
    Complex  power(const int& i)
    {
        return (this->real==0 && this->imag==0)? Complex(0,0) : Complex(pow(this->Norm(),(double)i), (double)i*this->Angle(), true);
    }
    Complex  power(const double& d)
    {
        return (this->real==0 && this->imag==0)? Complex(0,0) : Complex(pow(this->Norm(),d), d*this->Angle(), true);
    }
    Complex  power(const Real& r)
    {
        return (this->real==0 && this->imag==0)? Complex(0,0) : Complex(pow(this->Norm(),r.val), r.val*this->Angle(), true);
    }
    Complex  power(const Complex& c)
    {
        double n,a;
        return (this->real==0 && this->imag==0)? Complex(0,0) : ( n=this->Norm(), a=this->Angle(),
        Complex(pow(n,c.real)*exp(-a*c.imag), a*c.real+log(n)*c.imag, true) );
    }
    private:
        double norm;
        double angle;
};

// constructor anf destructors
Complex::Complex()
{
    this->real = 0.0;
    this->imag = 1.0;
}

Complex::Complex(const double& r, const double& im)
{
    this->real = r;
    this->imag = im;
}

Complex::Complex(const double& r, const double& im, const bool& ispolar)
{
    if(ispolar)
    {
        this->real = r*cos(im);
        this->imag = r*sin(im);
    }else
    {
        this->real = r;
        this->imag = im;
    }
}

Complex::Complex(const Real& r, const Real& im)
{
    this->real = r.val;
    this->imag = im.val;
}

Complex::Complex(const Real& r, const Real& im, const bool& ispolar)
{
    if(ispolar)
    {
        this->real = r.val*cos(im.val);
        this->imag = r.val*sin(im.val);
    }else
    {
        this->real = r.val;
        this->imag = im.val;
    }
}

Complex::Complex(const Complex& c)
{
    this->real = c.real;
    this->imag = c.imag;
}

Complex::~Complex(){};

// basic function
Complex Complex::Conj()
{
    return Complex( this->real, -this->imag);
}

double  Complex::Angle()
{
    double ratio;
    ratio = this->imag / this->real;
    if(this->real > 0){
        return atan(ratio);
    }else if(this->real < 0)
    {
        return atan(ratio) + PI;
    }else
    {
        if(this->imag > 0)
        {
            return PI/2;
        }else if(this->imag < 0)
        {
            return -PI/2;
        }else{
            return NAN;
        }
    }
}

double  Complex::const_Angle() const
{
    double ratio;
    ratio = this->imag / this->real;
    if(this->real > 0){
        return atan(ratio);
    }else if(this->real < 0)
    {
        return atan(ratio) + PI;
    }else
    {
        if(this->imag > 0)
        {
            return PI/2;
        }else if(this->imag < 0)
        {
            return -PI/2;
        }else{
            return NAN;
        }
    }
}

double  Complex::Norm()
{
    return sqrt( this->imag * this->imag + this->real*this->real );
}

double  Complex::const_Norm() const
{
    return sqrt( this->imag * this->imag + this->real*this->real );
}

double  Complex::Norm2()
{
    return this->imag * this->imag + this->real*this->real;
}

// overload operators
// ints

Complex Complex::operator+(const int& i)
{
    return Complex(this->real+(double)i, this->imag);
}

Complex Complex::operator-(const int& i)
{
    return Complex(this->real-(double)i, this->imag);
}


Complex Complex::operator*(const int& i)
{
    return Complex(this->real*(double)i, this->imag*(double)i );
}


Complex Complex::operator/(const int& i)
{
    return Complex(this->real/(double)i, this->imag/(double)i);
}

// doubles
Complex Complex::operator+(const double& d)
{
    return Complex(this->real+d, this->imag);
}


Complex Complex::operator-(const double& d)
{
    return Complex(this->real-d, this->imag);
}


Complex Complex::operator*(const double& d)
{
    return Complex(this->real*d, this->imag*d);
}


Complex Complex::operator/(const double& d)
{
    return Complex(this->real/d, this->imag/d);
}


// Real

Complex Complex::operator+(const Real& r)
{
    return Complex(this->real+r.val, this->imag);
}

Complex Complex::operator-(const Real& r)
{
    return Complex(this->real-r.val, this->imag);
}

Complex Complex::operator*(const Real& r)
{
    return Complex(this->real*r.val, this->imag*r.val);
}

Complex Complex::operator/(const Real& r)
{
    return Complex(this->real/r.val, this->imag/r.val);
}

// Complex
Complex Complex::operator+(const Complex& c)
{
    return Complex(this->real + c.real, this->imag + c.imag);
}

Complex Complex::operator-(const Complex& c)
{
    return Complex(this->real - c.real, this->imag - c.imag);
}

Complex Complex::operator-()
{
    return Complex(-this->real, -this->imag);
}

Complex Complex::operator*(const Complex& c)
{
    return Complex(this->real * c.real - this->imag * c.imag, this->real * c.imag + this->imag * c.real);
}

Complex Complex::operator/(const Complex& c)
{
    return Complex(this->real * c.real + this->imag * c.imag, -this->real * c.imag + 
    this->imag * c.real)/( c.real * c.real + c.imag* c.imag);
}

#endif


////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##### ####  #### #####
      #   #___  #___   #
      #   #        #   #
      #   ####  ####   #
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////

/*
int main(int argc, char* argv[])
{
    int i=1;
    double d=EE;
    std::cout << "DOING TEST--------------------------------------" << std::endl;
    std::cout << "DEFINES i     = " << i << std::endl;
    std::cout << "DEFINES d     = " << d << std::endl;

    // TEST INITIALIZATION
    Real r1 = Real(i);
    std::cout << "REAL r1       = " << r1 << std::endl;
    Real r2 = Real(d);
    std::cout << "REAL r2       = " << r2 << std::endl;
    Complex c1 = Complex();
    std::cout << "COMPLEX c1    = " << c1 << std::endl;
    Complex c2 = Complex(i,d);
    std::cout << "COMPLEX c2    = " << c2 << std::endl;
    Complex c3 = Complex(3,4);
    std::cout << "COMPLEX c3    = " << c3 << std::endl;
    
    // TEST BAISC OPERATOR OF REAL
    std::cout << "REAL d +r1    = " << d +r1 << std::endl;
    std::cout << "REAL i -r2    = " << i -r2 << std::endl;
    std::cout << "REAL d *r1    = " << d *r1 << std::endl;
    std::cout << "REAL i /r2    = " << i /r2 << std::endl;
    std::cout << "REAL r1+r2    = " << r1+r2 << std::endl;
    std::cout << "REAL r1-r2    = " << r1-r2 << std::endl;
    std::cout << "REAL r1*r2    = " << r1*r2 << std::endl;
    
    // TEST BASIC FUCTION OF COMPLEX
    std::cout << "Norm  Ang c1  = " << c1.Norm() << " " << c1.Angle() << std::endl;
    std::cout << "Norm2 Ang c3  = " << c3.Norm2() << " " << c3.Angle() << std::endl;
    std::cout << "c2, c2*       = " << c2 << ", " << c2.Conj() << std::endl;
    std::cout << "COMPLEX d +c3 = " << d +c3 << std::endl;
    std::cout << "COMPLEX i -c3 = " << i -c3 << std::endl;
    std::cout << "COMPLEX d *c3 = " << d *c3 << std::endl;
    std::cout << "COMPLEX i /c3 = " << i /c3 << std::endl;
    std::cout << "COMPLEX c2+c3 = " << c2+c3 << std::endl;
    std::cout << "COMPLEX c2-c3 = " << c2-c3 << std::endl;
    std::cout << "COMPLEX c2*c3 = " << c2*c3 << std::endl;
    std::cout << "COMPLEX c2/c3 = " << c2/c3 << std::endl;
    std::cout << "COMPLEX (-c3+c1)*c2   = " << (-c3+c1)*c2 << std::endl;
    
    // TEST POWER OPERATOR
    std::cout << "REAL r2.power(EE)      = " << r2.power(EE) << std::endl;
    std::cout << "REAL r2.power(r2)     = " << r2.power(r2) << std::endl;
    std::cout << "REAL power(r2, r2)    = " << power(r2, r2) << std::endl;
    std::cout << "REAL power(d, r2)     = " << power(d, r2) << std::endl;
    std::cout << "COMPLEX c2.power(d)   = " << c2.power(d) << std::endl;
    std::cout << "COMPLEX c2.power(r2)  = " << c2.power(r2) << std::endl;
    std::cout << "COMPLEX power(c2, r2) = " << power(c2, r2) << std::endl;
    std::cout << "COMPLEX c2.power(c2)  = " << c2.power(c2) << std::endl;
    std::cout << "COMPLEX power(c2, PI) = " << power(c2, PI) << std::endl;
    std::cout << "COMPLEX power(d, c2)  = " << power(d, c2) << std::endl;
    std::cout << "COMPLEX power(EE, c2) = " << power(EE, c2) << std::endl;
    
    
    std::cout << "DONE TEST!--------------------------------------" << std::endl;
    return 0;
} */










