/**
 * @file
 * @brief C++ interface for RAM subroutines
 */

#include "RamTypes.hh"

#include <complex>

extern "C" {

/// Determine array sizes by reading some of file ram.in.
/** The file is closed before this function returns.
 */
void sizes( 
    int& mr,    ///< OUT: Maximum number of range steps
    int& mz,    ///< OUT: Maximum number of depth steps
    int& mp,    ///< OUT: Maximum number of terms in rational approximation
    const char* fname,  ///< IN: C string containing input filename (ram.in)
    const int& fsize    ///< IN: length of filename = strlen(fname)
);
    
/// Initialize the parameters, acoustic field, and matrices
/**
 * Open the input and output files on I/O Units 1-3.
 * Unit 1, 'ram.in', is the input: a text file as specified in the User's Guide.
 * Unit 2, 'tl.line', is a text output file that will get the dB transmission loss
 * versus range (2 columns) at the target depth.
 * Unit 3, 'tl.grid', is a Fortran unformatted file that will get the dB
 * transmission loss on a grid of ranges and depths.
 *
 * Read the scenario parameters, the bathymetry versus range, and
 * the sound speed profile at the initial range (zero).
 * Then compute the starting field by calling selfs,
 * and write it in the output files.
 * Then compute the Pade coefficients and initialize the matrix etc.
 *
 * All arguments are outputs except mr, mz, and mp, which are
 * used as array dimensions. The caller must allocate all memory.
 * No global memory is used; all inputs and outputs are via arguments.
 *
 * In the argument documentation,
 * [mr], [mz], and [mp] designate dimensions for arrays versus
 * range, depth, and Pade coefficient index, respectively.
 * [mp,mz] designates dimensions for matrices versus Pade index and depth,
 * where the C/C++ index order is used: depth varies fastest.
 */
void setup(
    const int& mr,  ///< Maximum number of range steps
    const int& mz,  ///< Maximum number of depth steps
    int& nz,        ///< Number of depth steps (zmax/dz-0.5)
    const int& mp,  ///< Maximum number of terms in rational approximation
    int& np,        ///< Number of terms in rational approximation
    int& ns,        ///< Number of stability constraints (1 or 2)
    int& mdr,       ///< Current range index (set to 0)
    int& ndr,       ///< Range decimation factor for tl.grid (1=no decimation)
    int& ndz,       ///< Depth decimation factor for tl.grid (1=no decimation)
    int& iz,        ///< Depth index of bottom
    int& nzplt,     ///< Maximum depth index in output to tl.grid (zmplt/dz-0.5)
    int& lz,        ///< Number of depths output to tl.grid (~nzplt/ndz)
    int& ib,        ///< Current index in bathymetry array (set to 1)
    int& ir,        ///< Depth index of receiver
    Ram::Real& dir,     ///< Depth of receiver below the ir'th depth step
    Ram::Real& dr,      ///< Range step (m)
    Ram::Real& dz,      ///< Depth grid spacing (m)
    Ram::Real& pi,      ///< Math constant 3.14... (4.0*atan(1.0))
    Ram::Real& eta,     ///< Attenuation conversion factor: 1.0/(40.0*pi*alog10(exp(1.0)))
    Ram::Real& eps,     ///< A very small number (1.0e-20)
    Ram::Real& omega,   ///< Angular frequency (radians per second): 2.0*pi*freq
    Ram::Real& rmax,    ///< Maximum range (m)
    Ram::Real& c0,      ///< Reference sound speed (m/s)
    Ram::Real& k0,      ///< Reference wave number (m^-1): omega/c0
    Ram::Real& r,       ///< Current range (m) (set to dr)
    Ram::Real& rp,      ///< Range of profile update (m)
    Ram::Real& rs,      ///< Maximum range of stability constraints (m)
    Ram::Real rb[],     ///< [mr] Range of bathymetry point (m) vs. range
    Ram::Real zb[],     ///< [mr] Depth of bathymetry point (m) vs. range
    Ram::Real cw[],     ///< [mz] Sound speed in water column (m/s) vs. depth
    Ram::Real rhow[],   ///< [mz] Density in water (g/cc) vs. depth
    Ram::Real attnw[],   ///< [mz] Attenuation in water (dB/wavelength) vs. depth
    Ram::Real cb[],     ///< [mz] Sound speed in sediment (m/s) vs. depth
    Ram::Real rhob[],   ///< [mz] Density in sediment (g/cc) vs. depth
    Ram::Real attn[],   ///< [mz] Attenuation in sediment (dB/wavelength) vs. depth
    Ram::Real alpw[],   ///< [mz] Energy conversion factor sqrt(cw/c0) vs. depth
    Ram::Real alpb[],   ///< [mz] Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    Ram::Complex ksq[],
        ///< [mz] Complex relative wave number squared (ksqw or ksqb) vs. depth
    Ram::Complex ksqw[],
        ///< [mz] Relative wave number squared (kw^2-k0^2) in water vs. depth
    Ram::Complex ksqb[], 
        ///< [mz] Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
    Ram::Real f1[], ///< [mz] Matrix input: rho/alpha vs. depth
    Ram::Real f2[], ///< [mz] Matrix input: rho vs. depth
    Ram::Real f3[], ///< [mz] Matrix input: ksq vs. depth
    Ram::Complex u[], ///< [mz] The answer: Complex pressure/ksq vs. depth
    Ram::Complex v[], ///< [mz] Right-hand side of matrix equation, vs. depth
    Ram::Complex r1[], ///< [mp,mz] LHS matrix elements below the diagonal?
    Ram::Complex r2[], ///< [mp,mz] LHS matrix elements on the diagonal?
    Ram::Complex r3[], ///< [mp,mz] LHS matrix elements above the diagonal?
    Ram::Complex s1[], ///< [mp,mz] RHS matrix elements below the diagonal?
    Ram::Complex s2[], ///< [mp,mz] RHS matrix elements on the diagonal?
    Ram::Complex s3[], ///< [mp,mz] RHS matrix elements above the diagonal?
    Ram::Complex pd1[], ///< [mp] Numerator Pade coefficients?
    Ram::Complex pd2[], ///< [mp] Denominator Pade coefficients?
    float tlg[], ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    float preal[], ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    float pimag[], ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    const char* incstr,   ///< IN: C string containing input param file (ram.in)
    const int& ilen,      ///< IN: length of input filename = strlen(incstr)
    const char* linecstr, ///< IN: C string containing TL at receiver file (tl.line)
    const int& llen,      ///< IN: length of line filename = strlen(linecstr)
    const char* gridcstr, ///< IN: C string containing TL field file (tl.grid)
    const int& glen,      ///< IN: length of grid filename = strlen(gridcstr)
    const int& isrc_type, ///< IN: source: 0 for self-starter, 1 for gaussian tapered point source
    Ram::Real& x0,        ///< Range of source (m)
    Ram::Real& x1,        ///< Range of starting field [for tapered pt src] (m)
    Ram::Real& theta0,    ///< For tapered point source: beam angle (deg)
    Ram::Real& dtheta     ///< For tapered point source: beam width (deg)
);


/// Read and interpolate sound speed, density, and attenuation (water & ground) vs altitude
/** The input data come from Fortran unit 1, which must be already open. */
void profl(
    const int& mz,      ///< Maximum number of depth steps
    const int& nz,      ///< Number of depth steps ((amax-amin)/dz-0.5)
    const Ram::Real& dz,    ///< Depth grid spacing (m)
    const Ram::Real& eta,   ///< Attenuation conversion factor: 1.0/(40.0*pi*alog10(exp(1.0)))
    const Ram::Real& omega, ///< Angular frequency (radians per second): 2.0*pi*freq
    const Ram::Real& rmax,  ///< Maximum range (m)
    const Ram::Real& c0,    ///< Reference sound speed (m/s)
    const Ram::Real& k0,    ///< Reference wave number (m^-1): omega/c0
    Ram::Real& rp,      ///< OUT: Range of profile update (m)
    Ram::Real cw[],     ///< [mz] OUT: Sound speed in water column (m/s) vs. depth
    Ram::Real rhow[],   ///< [mz] OUT: Density in water (g/cc) vs. depth
    Ram::Real attnw[],   ///< [mz] OUT: Attenuation in water (dB/wavelength) vs. depth
    Ram::Real cb[],     ///< [mz] OUT: Sound speed in sediment (m/s) vs. depth
    Ram::Real rhob[],   ///< [mz] OUT: Density in sediment (g/cc) vs. depth
    Ram::Real attn[],   ///< [mz] OUT: Attenuation in sediment (dB/wavelength) vs. depth
    Ram::Real alpw[],   ///< [mz] Energy conversion factor sqrt(cw/c0) vs. depth
    Ram::Real alpb[],   ///< [mz] OUT: Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    Ram::Complex ksqw[],
        ///< [mz] OUT: Relative wave number squared (ka^2-k0^2) in water vs. depth
    Ram::Complex ksqb[] 
        ///< [mz] OUT: Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
);


/// Read and interpolate a profile of any real value vs. altitude
/** The input data come from Fortran unit 1, which must be already open. */
void zread(
    const int& mz,       ///< Maximum number of depth steps
    const int& nz,       ///< Number of depth steps ((amax-amin)/dz-0.5)
    const Ram::Real& dz, ///< Depth grid spacing (m)
    Ram::Real prof[]     ///< [mz] OUT: The profile read and interpolated vs. z
);


/// Build the tridiagonal matrices
void matrc(
    const int& mz,      ///< Maximum number of depth steps
    const int& nz,      ///< Number of depth steps ((amax-amin)/dz-0.5)
    const int& mp,      ///< Maximum number of terms in rational approximation
    const int& np,      ///< Number of terms in rational approximation
    const int& iz,      ///< Depth index of bottom
    const int& jz,      ///< Previous depth index of bottom
    const Ram::Real& dz,    ///< Depth grid spacing (m)
    const Ram::Real& k0,    ///< Reference wave number (m^-1): omega/c0
    const Ram::Real rhow[], ///< [mz] Density in water (g/cc) vs. depth
    const Ram::Real rhob[], ///< [mz] Density in sediment (g/cc) vs. depth
    const Ram::Real alpw[], ///< [mz] Energy conversion factor sqrt(rhoa*ca/c0) vs. depth
    const Ram::Real alpb[], ///< [mz] Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    Ram::Complex ksq[],
        ///< [mz] OUT: Complex relative wave number squared (ksqw or ksqb) vs. depth
    const Ram::Complex ksqw[],
        ///< [mz] Relative wave number squared (kw^2-k0^2) in air vs. depth
    const Ram::Complex ksqb[], 
        ///< [mz] Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
    Ram::Real f1[],     ///< [mz] Matrix input: rho/alpha vs. depth
    Ram::Real f2[],     ///< [mz] Matrix input: rho vs. depth
    Ram::Real f3[],     ///< [mz] Matrix input: alpha vs. depth
    Ram::Complex r1[],  ///< [mp,mz] LHS matrix elements below the diagonal?
    Ram::Complex r2[],  ///< [mp,mz] LHS matrix elements on the diagonal?
    Ram::Complex r3[],  ///< [mp,mz] LHS matrix elements above the diagonal?
    Ram::Complex s1[],  ///< [mp,mz] RHS matrix elements below the diagonal?
    Ram::Complex s2[],  ///< [mp,mz] RHS matrix elements on the diagonal?
    Ram::Complex s3[],  ///< [mp,mz] RHS matrix elements above the diagonal?
    Ram::Complex pd1[], ///< [mp] Numerator Pade coefficients?
    Ram::Complex pd2[]  ///< [mp] Denominator Pade coefficients?
);


/// Solve the nested tridiagonal matrix equations
/**
 * For each set of Pade coefficients j in [1,np], do the following:
 *
 * v = s_j u
 *
 * u = r_j^{-1} v
 *
 * where s_j and r_j are tridiagonal matrices, and u and v are vectors,
 * all of which have depth as their 'hidden' index.
 */
void solve(
    const int& mz,      ///< Maximum number of depth steps
    const int& nz,      ///< Number of depth steps ((amax-amin)/dz-0.5)
    const int& mp,      ///< Maximum number of terms in rational approximation
    const int& np,      ///< Number of terms in rational approximation
    const int& iz,      ///< Depth index of bottom
    Ram::Complex u[],   ///< [mz] The answer: Complex pressure/ksq vs. depth
    Ram::Complex v[],   ///< [mz] Scratch: Right-hand side of matrix equation, vs. depth
    const Ram::Complex r1[],    ///< [mp,mz] LHS matrix elements below the diagonal?
    const Ram::Complex r2[],    ///< [mp,mz] LHS matrix elements on the diagonal?
    const Ram::Complex r3[],    ///< [mp,mz] LHS matrix elements above the diagonal?
    const Ram::Complex s1[],    ///< [mp,mz] RHS matrix elements below the diagonal?
    const Ram::Complex s2[],    ///< [mp,mz] RHS matrix elements on the diagonal?
    const Ram::Complex s3[]     ///< [mp,mz] RHS matrix elements above the diagonal?
);


/// Update the matrix for the current range step
/**
 * I think the main outputs are the elements of the two tridiagonal matrices
 * for each Pade index: r1, r2, r3, s1, s2, and s3.
 * This is done for each range step.
 * Other arguments may also be modified along the way.
 * I haven't combed through them to figure out which ones remain constant.
 * 
 * In the simplest case, with no change in ocean depth or sound speed,
 * nothing is done.
 *
 * If the ocean depth changes by enough to change the depth step of the bottom,
 * the work is done by calling matrc, which updates only the elements between the
 * old and new bottom depth steps.
 *
 * If the sound speed profile changes, profl is called to update the profile,
 * and then matrc is called to update the matrices for all depths.
 *
 * If the range exceeds rs, the Pade coefficients are re-computed without the
 * stability constraints, and matrc is called to update the matrices.
 * Then rs is set very large so it won't happen again.
 */
void updat(
    const int& mr,  ///< Maximum number of range steps
    const int& mz,  ///< Maximum number of depth steps
    const int& nz,  ///< Number of depth steps (zmax/dz-0.5)
    const int& mp,  ///< Maximum number of terms in rational approximation
    const int& np,  ///< Number of terms in rational approximation
    int& iz,        ///< Depth index of bottom
    int& ib,        ///< Current index in bathymetry array (set to 1)
    Ram::Real& dr,      ///< Range step (m)
    Ram::Real& dz,      ///< Depth grid spacing (m)
    Ram::Real& eta,     ///< Attenuation conversion factor: 1.0/(40.0*pi*alog10(exp(1.0)))
    Ram::Real& omega,   ///< Angular frequency (radians per second): 2.0*pi*freq
    Ram::Real& rmax,    ///< Maximum range (m)
    Ram::Real& c0,      ///< Reference sound speed (m/s)
    Ram::Real& k0,      ///< Reference wave number (m^-1): omega/c0
    Ram::Real& r,       ///< Current range (m) (set to dr)
    Ram::Real& rp,      ///< Range of profile update (m)
    Ram::Real& rs,      ///< Maximum range of stability constraints (m)
    Ram::Real rb[],     ///< [mr] Range of bathymetry point (m) vs. range
    Ram::Real zb[],     ///< [mr] Depth of bathymetry point (m) vs. range
    Ram::Real cw[],     ///< [mz] Sound speed in water column (m/s) vs. depth
    Ram::Real rhow[],   ///< [mz] Density in water (g/cc) vs. depth
    Ram::Real attnw[],   ///< [mz] Attenuation in water (dB/wavelength) vs. depth
    Ram::Real cb[],     ///< [mz] Sound speed in sediment (m/s) vs. depth
    Ram::Real rhob[],   ///< [mz] Density in sediment (g/cc) vs. depth
    Ram::Real attn[],   ///< [mz] Attenuation in sediment (dB/wavelength) vs. depth
    Ram::Real alpw[],   ///< [mz] Energy conversion factor sqrt(rhow*cw/c0) vs. depth
    Ram::Real alpb[],   ///< [mz] Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    Ram::Complex ksq[],
        ///< [mz] Complex relative wave number squared (ksqw or ksqb) vs. depth
    Ram::Complex ksqw[],
        ///< [mz] Relative wave number squared (kw^2-k0^2) in water vs. depth
    Ram::Complex ksqb[], 
        ///< [mz] Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
    Ram::Real f1[], ///< [mz] Matrix input: rho/alpha vs. depth
    Ram::Real f2[], ///< [mz] Matrix input: rho vs. depth
    Ram::Real f3[], ///< [mz] Matrix input: ksq vs. depth
    Ram::Complex r1[], ///< [mp,mz] LHS matrix elements below the diagonal?
    Ram::Complex r2[], ///< [mp,mz] LHS matrix elements on the diagonal?
    Ram::Complex r3[], ///< [mp,mz] LHS matrix elements above the diagonal?
    Ram::Complex s1[], ///< [mp,mz] RHS matrix elements below the diagonal?
    Ram::Complex s2[], ///< [mp,mz] RHS matrix elements on the diagonal?
    Ram::Complex s3[], ///< [mp,mz] RHS matrix elements above the diagonal?
    Ram::Complex pd1[], ///< [mp] Numerator Pade coefficients?
    Ram::Complex pd2[]  ///< [mp] Denominator Pade coefficients?
);


/// The self-starter: Compute the starting field
void selfs(
    const int& mz,      ///< Maximum number of depth steps
    const int& nz,      ///< Number of depth steps ((amax-amin)/dz-0.5)
    const int& mp,      ///< Maximum number of terms in rational approximation
    const int& np,      ///< Number of terms in rational approximation
    const int& ns,      ///< Number of stability constraints (1 or 2)
    const int& iz,      ///< Depth index of bottom
    const Ram::Real& zs,   ///< Source depth (m)
    const Ram::Real& dr,    ///< Range step (m)
    const Ram::Real& dz,    ///< Depth grid spacing (m)
    const Ram::Real& pi,    ///< Math constant 3.14... (4.0*atan(1.0))
    const Ram::Real& c0,    ///< Reference sound speed (m/s)
    const Ram::Real& k0,    ///< Reference wave number (m^-1): omega/c0
    const Ram::Real rhow[], ///< [mz] Density in water (g/cc) vs. depth
    const Ram::Real rhob[], ///< [mz] Density in sediment (g/cc) vs. depth
    const Ram::Real alpw[], ///< [mz] Energy conversion factor sqrt(rhoa*ca/c0) vs. depth
    const Ram::Real alpb[], ///< [mz] Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    Ram::Complex ksq[],
        ///< [mz] OUT: Complex relative wave number squared (ksqw or ksqb) vs. depth
    const Ram::Complex ksqw[],
        ///< [mz] Relative wave number squared (kw^2-k0^2) in water vs. depth
    const Ram::Complex ksqb[], 
        ///< [mz] Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
    Ram::Real f1[],     ///< [mz] Matrix input: rho/alpha vs. depth
    Ram::Real f2[],     ///< [mz] Matrix input: rho vs. depth
    Ram::Real f3[],     ///< [mz] Matrix input: alpha vs. depth
    Ram::Complex u[],   ///< [mz] The answer: Complex pressure/ksq vs. depth
    Ram::Complex v[],   ///< [mz] Scratch: Right-hand side of matrix equation, vs. depth
    Ram::Complex r1[],  ///< [mp,mz] LHS matrix elements below the diagonal
    Ram::Complex r2[],  ///< [mp,mz] LHS matrix elements on the diagonal
    Ram::Complex r3[],  ///< [mp,mz] LHS matrix elements above the diagonal
    Ram::Complex s1[],  ///< [mp,mz] RHS matrix elements below the diagonal
    Ram::Complex s2[],  ///< [mp,mz] RHS matrix elements on the diagonal
    Ram::Complex s3[],  ///< [mp,mz] RHS matrix elements above the diagonal
    Ram::Complex pd1[], ///< [mp] Numerator Pade coefficients
    Ram::Complex pd2[]  ///< [mp] Denominator Pade coefficients
);

/// Compute the starting field (2D) for a  Gaussian tapered point source 
void pt_src_gtaper(
    const int& mz,      ///< Maximum number of depth steps
    const int& nz,      ///< Number of depth steps ((amax-amin)/dz-0.5)
    Ram::Real& dz,      ///< Depth grid spacing (m)
    const Ram::Real& freq,  ///< Source freq (Hz)
    const Ram::Real& c0,    ///< Reference sound speed (m/s)
    const Ram::Real& zs,    ///< Source depth (m)
    const Ram::Real& x0,    ///< Source range (m)
    const Ram::Real& x1,    ///< Range of computed vertical field (m)
    const Ram::Real& theta0,    ///< Beam angle (deg)  (-ve is towards surface)
    const Ram::Real& dtheta,    ///< Angular step (deg)
    Ram::Complex p[]   ///< [mz] OUT: Complex pressure
);


/// Output the transmission loss
void outpt(
    const int& mz,  ///< Maximum number of depth steps
    int& mdr,       ///< Current range index (incremented by 1)
    int& ndr,       ///< Range decimation factor for tl.grid (1=no decimation)
    int& ndz,       ///< Depth decimation factor for tl.grid (1=no decimation)
    int& iz,        ///< Depth index of bottom
    int& nzplt,     ///< Maximum depth index in output to tl.grid (zmplt/dz-0.5)
    int& lz,        ///< Number of depths output to tl.grid (~nzplt/ndz)
    int& ir,        ///< Depth index of receiver
    Ram::Real& dir,     ///< Depth of receiver below the ir'th depth step
    Ram::Real& eps,     ///< A very small number (1.0e-20)
    Ram::Real& r,       ///< Current range (m)
    Ram::Real f3[],     ///< [mz] Matrix input: ksq vs. depth
    Ram::Complex u[], ///< [mz] The answer: Complex pressure/ksq vs. depth
    float tlg[],     ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    float preal[],     ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    float pimag[]     ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
);


/// Compute the coefficients of the rational approximation
void epade(
    const int& mp,       ///< Maximum number of terms in rational approximation
    const int& np,       ///< Number of terms in rational approximation
    const int& ns,       ///< Number of stability constraints (1 or 2)
    const int& ip,       ///< =1 if called from setup or updat, 2 if called from selfs
    const Ram::Real& k0, ///< Reference wave number (m^-1): omega/c0
    const Ram::Real& c0, ///< Reference sound speed (m/s)
    const Ram::Real& dr, ///< Range step (m)
    Ram::Complex pd1[],  ///< [mp] Numerator Pade coefficients?
    Ram::Complex pd2[]   ///< [mp] Denominator Pade coefficients?
);


/// The operator function
/** Not used; replaced by Pade series
Ram::DoubleComplex g(
    const Ram::Double& sig, ///< k0*dr
    const Ram::Double& x,   ///< The operand
    const Ram::Double& alp, ///< 0 or -0.25
    const Ram::Double& nu   ///< 0 or 1.0d0
);
*/

/// The derivatives of the operator function at x=0.
void deriv(
    const int& m,           ///< Dimension of in, dg, dh1, dh2, dh3
    const int& n,           ///< Number of derivatives required (2*np)
    const Ram::Double& sig, ///< k0*dr
    const Ram::Double& alp, ///< 0 or -0.25
    Ram::Double dg[],       ///< [m]
    Ram::Double dh1[],      ///< [m]
    Ram::Double dh2[],      ///< [m]
    Ram::Double dh3[],      ///< [m]
    const Ram::Double bin[], ///< [m*m] Binomial coefficients
    const Ram::Double& nu   ///< 0 or 1.0d0
);


/// Solve a complex linear system by Gaussian elimination
void gauss(
    const int& m,           ///< Dimension of a, b
    const int& n,           ///< Index range of a, b
    Ram::DoubleComplex a[], ///< [m*m] The matrix to be inverted; n output, its inverse?
    Ram::DoubleComplex b[]  ///< [m] On input, the RHS vector; on output, the result
);


/// In Gaussian elimination, rows are interchanged for stability
void pivot(
    const int& m,            ///< Dimension of a, b
    const int& n,            ///< Index range of a, b
    const int& i,            ///< Index of pivot element
    Ram::DoubleComplex a[],  ///< [m*m] The matrix to be inverted; on output, its inverse?
    Ram::DoubleComplex b[]   ///< [m] On input, the RHS vector; on output, the result
);


/// Find the roots of a complex polynomial
void fndrt(
    Ram::DoubleComplex a[], ///< [m] Coefficients of the polynomial (over-written)
    const int& n,           ///< Order of the polynomial
    Ram::DoubleComplex z[], ///< OUT: Roots of the polynomial
    const int& m            ///< Dimension of a, z
);


/// Find a root of a polynomial of degree n > 2 by Laguerre's method.
void guerre(
    const Ram::DoubleComplex a[], ///< [m] Coefficients of the polynomial
    const int& n,           ///< Order of the polynomial
    const int& m,           ///< Dimension of a, z
    Ram::DoubleComplex& z,  ///< OUT: A root of the polynomial
    const Ram::Double& err, ///< Tolerance for root z
    const int& nter         ///< Maximum number of iterations
);


/// Close the input and output files used by RAMair (units 1, 2, and 3)
void close_ram_files();

}   // end extern "C"

