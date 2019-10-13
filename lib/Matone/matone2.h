////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DEFINITION OF MATRIX AND TENSOR TYPE

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////

template <typename T>
class mat: public Vec
{
    public:
        // constructor
        mat(int d1, int d2);
        // copy constructor
        mat(const mat& mat_copy);
        // move constructor
        mat(C&& mat_move);
        // copy assignment
        mat& operator= (const mat& mat_copy);
        // move assignment
        mat& operator= (mat&& mat_move);
        // destrcutor
        ~mat();
        
        // operations
        // refering
        T& operator()(int d1, int d2);
        // addition
        mat& operator+(const mat& m);
        // subtraction
        mat& operator+(const mat& m);
        // dot production
        mat& operator+(const mat& m);
        // trace
        // transpose
        // 
        
    private:
        static int err_level;
        int rows;
        int cols;
        T* vals;
};

template <typename T>

mat::mat<T>(int d1, int d2)
{
    this->rows = d1;
    this->cols = d2;
    vals = new T[d1*d2];
}

mat::mat<T>(int d1, int d2)
{
    this->rows = d1;
    this->cols = d2;
    vals = new T[d1*d2];
}

mat<T> mat::ones<T>(int d1, int d2)
{
    mat<T> example_ones = mat<T>(d1, d2);
    T** p = vals;
    for(int i=0; i<example_ones.rows; i++)
    {
        for(int j=0; j<example_ones.cols; j++)
        {
            *(p + i*this->rows + j)= (i==j)? T(1) : T(0) ;
        }
    }
    return example_ones;
}

mat::<T> & operator[](int &d1)
{
	return this->OrbC(idx);
}

mat::<T> & operator[](int &d1)
{
	return this->OrbC(idx);
}

////////////////////////////////////////////////////////////////////////
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DEFINITION OF REAL TYPE

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
////////////////////////////////////////////////////////////////////////


template <typename T>
class Tensor
{ 
}
