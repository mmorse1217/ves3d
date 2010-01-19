/**
 * @file   vesicle.h
 * @author Rahimian, Abtin <abtin@romario>
 * @date   Mon Jan 18 14:08:45 2010
 * 
 * @brief  The header function of the class vesicle.
 *
 */


/**
 * @class vesicle
 * @brief The main class
 *
 * vesicle is the class used to ...
 * 
 */
template <typename T> class vesicle
{
public:
    int nv; /**< The number of vesicles */
    int p; /**< The spherical harmonics basis degree.*/
    T* posVec; /**< The position vector for the grid points */

    /**
     * A constructor for the class vesicle.
     * 
     * @param nvIn Number of vesicles.
     * @param pIn Maximum degree of spherical harmonics.
     */
    vesicle(int nvIn, int pIn) : nv(nvIn), p(pIn)
    {
	posVec = new T[p];
	
	for(int ii=0;ii<p;++ii)
	    posVec[ii] = ii;
    }

    
    /**
     * The deconstructor for the class vesicle.
     * In the deconstructor the dynamic memory is freed. The position
     * vector is allocated dynamically, therefore it is deleted here. 
     */
    ~vesicle()
    {
	delete[] posVec;
    }

private:
    int np; /**< The number of points on each vesicle  */
    int vecLength; /**< The total length of the position vector */
    
};

