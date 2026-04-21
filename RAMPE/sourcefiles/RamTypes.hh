/**
 * @file
 * @brief Types used in C++ interface for RAM subroutines
 */

#include <complex>

/// Encapsulate the float/double precision choice here and in ramlib.f
namespace Ram {
    /// Type of RAM arguments that were originally real*4.
    /** To choose between float and double, define Real here, and
     * make compatible changes to the declaration section of module
     * ramlib near the beginning of file ramlib.h.
     */
    typedef double Real;
    
    /// Type of RAM arguments that were originally complex*8.
    /** Leave this alone; change only Real to choose precision. */
    typedef std::complex<Real> Complex;
    
    /// Type of RAM arguments that were originally complex*16.
    typedef std::complex<double> DoubleComplex;
    
    /// Type of RAM arguments that were originally real*8.
    typedef double Double;
}
