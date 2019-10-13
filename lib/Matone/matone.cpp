#ifndef MATONE_H
#define MATONE_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <exception>
//#include <random>
#include "mathtype.h"

////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DEFINITION OF VECTOR TYPE

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////

template <typename T>
class Vec
{
    public:
        Vec();
        Vec(unsigned int d);
        Vec(int d);
        Vec(const Vec& v_copy);
        //Vec(Vec&& v_move);
        ~Vec();
        //
        unsigned int const_len() const;
        T dot(const Vec& v);
        T mean();
        T sum();
        T max();
        T min();
        static Vec& range(const int& d, const double& shift=0, const double& interval=1)
        {
            Vec* pv = new Vec(d);
            for(int i=0; i<d; i++) *(pv->vals+i) = (double)i;
            (*pv) = shift + (*pv)*interval;
            return *pv;
        }
        static Vec& random(int d);

        Vec& sort(const int& i_from, const int& i_end, const bool& rev = false);
        
        Vec& operator+(const T& a);
        Vec& operator-(const T& a);
        Vec& operator*(const T& a);
        Vec& operator/(const T& a);
        Vec& operator=(const T& a);
        Vec& operator+(const Vec& v);
        Vec& operator-(const Vec& v);
        Vec& operator-();
        Vec& operator*(const Vec& v);
        Vec& operator/(const Vec& v);
        Vec& operator=(const Vec& v);
        T&   operator()(const int& d);
        Vec& operator()(const int& d1, const int& d2, const bool& iscopy=false);
        // Vec& operator=(Vec&& v); // move assignment
        friend std::ostream& operator<<(std::ostream& output, const Vec& v)
        {
            output << std::setiosflags(std::ios::left) << std::setprecision(8)<<"[ ";
            if(v.const_len() <= 20)
            {
                for(int i=0; i< v.const_len(); i++) output << *(v.vals+i) << "  ";
            }else
            {
                for(int i=0; i<4; i++) output << *(v.vals+i) << "  ";
                output << "... ";
                for(int i=v.const_len()-4; i<v.const_len(); i++) output << *(v.vals+i) << "  ";
            }
            output << "]";
            return output;
        }
        friend std::istream& operator>>(std::istream& input, Vec& v)
        {
            for(int i=0; i< v.const_len(); i++) input >> *(v.vals+i) ;
            return input;
        }
        
        inline friend Vec& operator+(const T& a, const Vec& v)
        {
            if(v.const_len() > 0)
            {
                Vec* pv1 = new Vec(v);
                for(int i=0; i< v.const_len(); i++) *(pv1->vals+i) = a + *(v.vals+i) ;
                return *pv1;
            }else
            {
                std::cerr << "Vec type without initization" << std::endl; exit(-1);
            }    
        }
        inline friend Vec& operator-(const T& a, const Vec& v)
        {
            if(v.const_len()>0)
            {
                Vec* pv1 = new Vec(v);
                for(int i=0; i< v.const_len(); i++) *(pv1->vals+i) = a - *(v.vals+i) ;
                return *pv1;
            }else
            {
                std::cerr << "Vec type without initization" << std::endl; exit(-1);
            }
        }
        inline friend Vec& operator*(const T& a, const Vec& v)
        {
            if(v.const_len()>0)
            {
                Vec* pv1 = new Vec(v);
                for(int i=0; i< v.const_len(); i++) *(pv1->vals+i) = *(v.vals+i) * a ;
                return *pv1;
            }else
            {
                std::cerr << "Vec type without initization" << std::endl; exit(-1);
            }
        }
    protected:
        unsigned int len;
        T* vals;
};

template <typename T>
Vec<T>::Vec()
{
    this->len = 0;
    this->vals = NULL;
}

template <typename T>
Vec<T>::Vec(unsigned int d)
{
    if(d > 0)
    {
        this->len = d;
        this->vals = new T[this->len];
        for(int i=0; i<this->len; i++) *(this->vals+i) = (T)0.0;
    }else
    {
        std::cout << "Vec initial length error" << std::endl;
    }
}

template <typename T>
Vec<T>::Vec(int d)
{
    if( d > 0 ){
        this->len = d;
        this->vals = new T[this->len];
        for(int i=0; i<this->len; i++) *(this->vals+i) = (T)0.0;
    }else
    {
         std::cout << "Vec initial length error" << std::endl;
    }
}

template <typename T>
Vec<T>::Vec(const Vec& v_copy)
{
    unsigned int common_len = v_copy.const_len();
    if(common_len > 0)
    {
        this->len = common_len ;
        this->vals = new T[this->len];
        for(int i=0; i<common_len ; i++) *(this->vals+i) = *(v_copy.vals+i);
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

/*
template <typename T>
Vec<T>::Vec(Vec&& v_move)
{
    unsigned int common_len = v_move.const_len();
    if(common_len > 0 && this->len == common_len) // how about inequal?
    {
        this->vals = v_move.vals;
        v_move.vals = NULL;
        return *this;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}*/

template <typename T>
Vec<T>::~Vec()
{
    delete[] this->vals;
}

// Functions
template <typename T>
unsigned int Vec<T>::const_len() const
{
    return this->len;
}

template <typename T>
T Vec<T>::dot(const Vec& v)
{
    T s = (T)0.0;
    if(this->len > 0 && this->len == v.const_len())
    {
        for(int i=0; i<this->len; i++) s += (*(this->vals+i)) * (*(v.vals+i));
        return s;
    }else
    {
        std::cerr << "Vec type without initization or Vec size dismarch" << std::endl; exit(-1);
    }
}

template <typename T>
T Vec<T>::mean()
{
    T s = T(0.0);
    if(this->len > 0)
    {
        for(int i=0; i<this->len; i++) s+=*(this->vals+i);
        s = s/((T)this->len);
        return s;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
T Vec<T>::sum()
{
    T s = (T)0.0;
    if(this->len > 0)
    {
        for(int i=0; i<this->len; i++) s+=*(this->vals+i);
        return s;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
T Vec<T>::min()
{
    T s = *this->vals;
    if(this->len > 0)
    {
        for(int i=1; i<this->len; i++)
        {
            if(s>*(this->vals+i)) s=*(this->vals+i);
        }
        return s;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
T Vec<T>::max()
{
    T s = *this->vals;
    if(this->len > 0)
    {
        for(int i=1; i<this->len; i++)
        {
            if(s < *(this->vals+i)) s=*(this->vals+i);
        }
        return s;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}


template <typename T>
Vec<T>& Vec<T>::sort(const int& i_from, const int& i_end, const bool& rev) //default rev=false, from=0, end=-1
{
    int from = (this->len + i_from)%(this->len); 
    int end = (this->len + i_end)%(this->len);
    T* pTi, * pTj;
    T key, tmp;
    bool searchforward;
    int m=from; int step = 10;
    
    if(this->len > 0)
    {
        pTi=(this->vals+from); pTj=(this->vals+end);
        //std::cout << "(0,-1) "<< *pTi << "  " << *pTj << std::endl;
        key = *pTi; searchforward = false;
        
        if(!rev)
        {
            // case 1
            if(from >= end) return *this;
            // case 2
            if(from+1 == end)
            {
                if(*pTi > *pTj)
                {
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                }
                return *this;
            }
            // case 3
            while(pTi!=pTj)
            {
                if(!searchforward)
                {
                    for(;*(pTj) >= key && pTj>pTi; pTj--);
                    if(pTi==pTj) break;
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                    searchforward = true;
                }else
                {
                    for(;*(pTi) <= key && pTi<pTj; pTi++) m++;
                    if(pTi==pTj) break;
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                    searchforward = false;
                }
            }
            //std::cout << "############################## CALL LEFT ###########################"<< std::endl;
            if(m>from+1)
            {
                this->sort(from, m-1, rev);
            }
            //std::cout << "############################## CALL RIGHT###########################"<< std::endl;
            if(m<end-1)
            {
                this->sort(m+1, end, rev);
            }
            return *this;
        }else
        {
            // case 1
            if(from >= end) return *this;
            // case 2
            if(from+1 == end)
            {
                if(*pTi < *pTj)
                {
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                }
                return *this;
            }
            // case 3
            while(pTi!=pTj)
            {
                if(!searchforward)
                {
                    for(;*(pTj) <= key && pTj>pTi; pTj--);
                    if(pTi==pTj) break;
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                    searchforward = true;
                }else
                {
                    for(;*(pTi) >= key && pTi<pTj; pTi++) m++;
                    if(pTi==pTj) break;
                    tmp = *pTj; *pTj = *pTi; *pTi = tmp;
                    searchforward = false;
                }
            }
            //std::cout << "############################## CALL LEFT ###########################"<< std::endl;
            if(m>from+1)
            {
                this->sort(from, m-1, rev);
            }
            //std::cout << "############################## CALL RIGHT###########################"<< std::endl;
            if(m<end-1)
            {
                this->sort(m+1, end, rev);
            }
            return *this;
        }
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}
              
// Operators
template <typename T>
Vec<T>& Vec<T>::operator+(const T& a)
{
    if(this->len > 0)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<this->len; i++) *(pv->vals+i) += a ;
        return *pv;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator-(const T& a)
{
    if(this->len > 0)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<this->len; i++) *(pv->vals+i) -= a ;
        return *pv;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator*(const T& a)
{
    if(this->len > 0)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<this->len; i++) *(pv->vals+i) *= a ;
        return *pv;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator/(const T& a)
{
    // definition of !a
    if(a==0) { std::cerr << "Error division by zero" << std::endl; exit(-1); }
    if(this->len> 0 && a!=0)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<this->len; i++) *(pv->vals+i) /= a ;
        return *pv;
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
    return *this;
}

template <typename T>
Vec<T>& Vec<T>::operator=(const T& a)
{
    // definition of !a
    if(this->len > 0)
    {
        try
        {
            Vec<T>* pv = new Vec<T>(*this);
            for(int i=0; i<this->len; i++) *(pv->vals+i) = a ;
            return *pv;
        }
        catch(std::exception& e)
        {
            std::cout << "Standard exception: " << e.what() << std::endl; exit(-1);
        }
    }else
    {
        std::cerr << "Vec type without initization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator+(const Vec& v)
{
    unsigned int common_len = v.const_len(); 
    if(this->len > 0 && this->len == common_len)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<common_len; i++) *(pv->vals+i) += *(v.vals+i) ;
        return *pv;
    }else
    {
        std::cerr << "Vec without initialization or Vec size dismarch" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator-(const Vec& v)
{
    unsigned int common_len = v.const_len(); 
    if(this->len > 0 && this->len == common_len)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<common_len; i++) *(pv->vals+i) -= *(v.vals+i) ;
        return *pv;
    }else
    {
        std::cerr << "Vec without initialization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator-()
{
    if(this->len > 0 )
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<this->len; i++) *(pv->vals+i) = - *(this->vals+i) ;
        return *pv;
    }else
    {
        std::cerr << "Vec without initialization" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator*(const Vec& v)
{
    unsigned int common_len = v.const_len(); 
    if(this->len > 0 && this->len == common_len)
    {
        Vec<T>* pv = new Vec<T>(*this);
        for(int i=0; i<common_len; i++) *(pv->vals+i) *= *(v.vals+i) ;
        return *pv;
    }else
    {
        std::cerr << "Vec without initialization" << std::endl; exit(-1);
    }
}

// NOT SUGGESTED!
template <typename T>
Vec<T>& Vec<T>::operator/(const Vec& v)
{
    unsigned int common_len = v.const_len(); 
    if(this->len > 0 && this->len == common_len)
    {
        // might divide zero
        try
        {
            Vec<T>* pv = new Vec<T>(*this);
            for(int i=0; i<common_len; i++) *(pv->vals+i) /= *(v.vals+i) ;
            return *pv;
        }
        catch(std::exception& e)
        {
            std::cout << "Standard exception: " << e.what() << std::endl; exit(-1);
        }
    }else
    {
        std::cerr << "vec operation miss march or not initialed" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>& Vec<T>::operator=(const Vec& v)
{
    // definition of !a
    if(v.const_len()>0 && this->len == v.const_len() )
    {
        for(int i=0; i<this->len; i++) *(this->vals+i) = *(v.vals+i);
        return *this;
    }else
    {
        std::cerr << "Vec assignment size dismarch" << std::endl; exit(-1);
    }
}

template <typename T>
T&  Vec<T>::operator()(const int& d)
{
    if(d>=0 && d < this->len)
    {
        return *(this->vals+d);
    }else if (d<0 && -d <= this->len)
    {
        return *(this->vals+this->len+d);
    }else
    {
        std::cerr << "Vec boundary error" << std::endl; exit(-1);
    }
}

template <typename T>
Vec<T>&  Vec<T>::operator()(const int& d1, const int& d2, const bool& iscopy ) // NOTE NOT A COPY METHOD
{
    int from = d1; int end = (d2+this->len)%this->len;
    if(from >= 0 && d2 < this->len && from <= end)
    {
        if(iscopy)
        {
            Vec* pv = new Vec(end-from);
            for(int i=0; i< pv->len; i++) *(pv->vals+i) = *(this->vals+from+i) ;
            return *pv;
        }else
        {
            Vec* pv = new Vec();
            pv->len = end-from;
            pv->vals = this->vals+from;
            return *pv;
        }
    }else
    {
        std::cerr << "Vec boundary error" << std::endl; exit(-1);
    }
}


/*
template <typename T>
Vec<T>& Vec<T>::operator=(Vec&& v)
{
    // definition of !a
    delete[] this->vals;
    this->vals = v.vals;
    this->len = v.const_len();
    v.vals = NULL;
    return *this;
}*/

////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##### ####  #### #####
      #   #___  #___   #
      #   #        #   #
      #   ####  ####   #
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////
typedef Vec<int>        Vec_i;
typedef Vec<double>     Vec_d;
typedef Vec<Real>       Vec_r;
typedef Vec<Complex>    Vec_c;

/*
int main(int argc, char* argv[])
{
    std::cout << "DOING TEST--------------------------------------" << std::endl;
    Vec_d v1 = Vec_d();
    std::cout << "v1 = Vec(), then v1 = " << v1 << std::endl;
    Vec_d v2 = Vec_d(10);
    std::cout << "v2 = Vec(10), then v2 = " << v2 << std::endl;
    Vec_d v3 = Vec_d::range(20);
    std::cout << "v3 = Vec::range(20), then v3 = " << v3 << std::endl;
    Vec_d v4 = Vec_d::range(40);
    std::cout << "v4 = Vec::range(40), then v4 = " << v4 << std::endl;

    // test not change the values
    std::cout << "2-v3*3 = " << 2-v3*3 << std::endl;
    std::cout << "v3 = Vec::range(20), then v3 = " << v3 << std::endl;
    
    // assignment must equal size, note following will report error
    //v2 = v3;
    //std::cout << "v2 = v3, then v2 = " << v2 << std::endl;
    
    
    // test assignment with equal size
    v3 = v3*(v3-9.5);
    std::cout << "v3 * (v3-9.5), then v3 = " << v3 << std::endl;
    
    // test statistic function
    std::cout << "v3.mean() = " << v3.mean() << std::endl;
    std::cout << "v3.sum() = " << v3.sum() << std::endl;
    std::cout << "v3.max() = " << v3.max() << std::endl;
    std::cout << "v3.min() = " << v3.min() << std::endl;
    
    // test index fucntion
    std::cout << "now v3 = " << v3 << std::endl;
    std::cout << "v3(0) = " << v3(0) << "  " << "v3(-1) = " << v3(-1) << std::endl;
    v3(0) = 100; v3(-1) = -100;
    std::cout << "now let v3(0) = 100, v3(-1) = -100" << std::endl;
    std::cout << "now v3 = " << v3 << std::endl;
    std::cout << "v3(0) = " << v3(0) << "  " << "v3(-1) = " << v3(-1) << std::endl;
    std::cout << "v3(0,10) = " << v3(0,10) << std::endl;
    std::cout << "v3(3,9) = " << v3(3,9) << std::endl;
    v3(3,6) = Vec_d::range(3,1,0.5);
    std::cout << "v3(3,6) = Vec_d::rangle(3), now v3 = " << v3 << std::endl;
    std::cout << "v3(3,6) = Vec_d::rangle(3), now v3.len = " << v3.const_len() << std::endl;
    
    // test sort function
    v3.sort(0,-1);
    std::cout << "v3.sort(0,-1), then v3 = " << v3 << std::endl;
    v3.sort(0,-1,true);
    std::cout << "v3.sort(0,-1,true), then v3 = " << v3 << std::endl;
    std::cout << "DONE TEST!--------------------------------------" << std::endl;
    return 0;
}*/

#endif

