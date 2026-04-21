/**
 * @file
 * @brief RamRunner: C++ class to encapsulate a run of RAM
 */

#ifndef RamRunner_h
#define RamRunner_h
// $Id: RamRunner.hh 7495 2015-05-13 19:12:14Z rpg $

#include "RamTypes.hh"

#include <string>
#include <vector>

/// One run of RAM, along one vertical slice, at one frequency
/**
 * This wrapper holds all of the state that is passed from one RAM
 * Fortran subroutine to another.
 * It greatly simplifies the caller's interface by limiting the input and
 * output variables to those that the caller is interested in.
 *
 * The intent is that multiple instances of this class could be run in parallel.
 * This won't work yet because it uses fixed Fortran unit numbers for I/O.
 */
class RamRunner {

public:
    /// Full Constructor taking inputs from constructor arguments, not a file
    /**
     * This constructor does NOT read any input files, nor does it open
     * any output files.
     * 
     * The input parameters are in the same order as in the RAM input file.
     */
    RamRunner(
        Ram::Real freq, ///< Acoustic frequency [Hz]
        Ram::Real zs,   ///< Source depth (m)
        Ram::Real zr,   ///< Receiver depth (m) for output
        Ram::Real rmax, ///< Maximum range (m)
        Ram::Real dr,   ///< Range step (m)
        int ndr,        ///< Range decimation factor for output (1=no decimation)
        Ram::Real zmax, ///< Maximum depth (m) of computational domain
        Ram::Real dz,   ///< Depth grid spacing (m)
        int ndz,        ///< Depth decimation factor for output (1=no decimation)
        Ram::Real zmplt,///< Maximum depth of output (m)
        Ram::Real c0,   ///< Reference sound speed [m/s]
        int np,         ///< Number of terms in rational approximation
        int ns,         ///< Number of stability constraints (1 or 2)
        Ram::Real rs,   ///< Maximum range of stability constraints (m)
        const std::vector<Ram::Real>& rb,   ///< Range of bathymetry point (m)
        const std::vector<Ram::Real>& zb,   ///< Depth of bathymetry point (m) vs. range
        const std::vector<Ram::Real>& zcw,  ///< Depth points for cw
        const std::vector<Ram::Real>& cw,   ///< Sound speed in water column (m/s) vs. depth
        const std::vector<Ram::Real>& zrhow,///< Depth points for rhow
        const std::vector<Ram::Real>& rhow, ///< Density in water (rel. water) vs. depth
        const std::vector<Ram::Real>& zattnw,///< Depth points for attnw
        const std::vector<Ram::Real>& attnw, ///< Attenuation in water (dB/wavelength) vs. depth
        const std::vector<Ram::Real>& zcb,  ///< Depth points for cb
        const std::vector<Ram::Real>& cb,   ///< Sound speed in sediment (m/s) vs. depth
        const std::vector<Ram::Real>& zrhob,///< Depth points for rhob
        const std::vector<Ram::Real>& rhob, ///< Density in sediment (rel. water) vs. depth
        const std::vector<Ram::Real>& zattn,///< Depth points for attn
        const std::vector<Ram::Real>& attn, ///< Attenuation in sediment (dB/wavelength) vs. depth
        int isrc_type,       ///< const int& isrc_type, /// 0 for self-starter, 1 for gaussian tapered point source
        Ram::Real x0,       ///< source range (m) [for tapered point source, else equals 0]
        Ram::Real x1,       ///< starting field range (m)
        Ram::Real theta0,   ///< Beam angle (deg)
        Ram::Real dtheta    ///< Beam width (deg)
    );

    /// Constructor needs the names of the input and output files
    /**
     * This constructor calls the original RAM's 'setup' program,
     * which reads all of its inputs from a file whose name is in inFile
     * and opens output files whose names are in liuneFile and gridFile.
     * Then it initializes the calculation using the 'self-starter'
     * in the Fortran subroutine 'selfs', and writes the first parts of
     * the output files by calling 'outpt'.
     * This is the initialization behavior of the original Fortran
     * 'ram' main program.
     */
    RamRunner(
        const std::string& inFile,      ///< Name of the input text file
        const std::string& lineFile,    ///< Name of the output line text file
        const std::string& gridFile     ///< Name of the output grid binary file
    );
    
    /// Destructor closes the file and deallocates memory
    ~RamRunner();
  
    /// Set up all of the profiles versus z
    /**
     * The input data for sound speed in water and bottom, bottom density,
     * and bottom attenuation are given at user-specified depth points.
     * The output profiles (RamRunner member data) are tabulated on a 
     * common fixed depth grid:
     * nz_ samples starting at depth 0 with spacing dz_.
     * 
     * This function interpolates the input profiles onto the output grid,
     * and uses the results to compute several more profiles (member data)
     * on the same grid. These include k^2 - k0^2 in the water and bottom,
     * and alpha (density*soundspeed/c0) in the water and bottom.
     */
    void SetProfiles(
        const std::vector<Ram::Real>& zcw,    ///< Depth points for cw
        const std::vector<Ram::Real>& cw,     ///< Sound speed in water column (m/s) vs. depth
        const std::vector<Ram::Real>& zrhow,  ///< Depth points for rhow
        const std::vector<Ram::Real>& rhow,   ///< [mz] Density in water (rel. water) vs. depth
        const std::vector<Ram::Real>& zattnw,  ///< Depth points for attnw
        const std::vector<Ram::Real>& attnw,    ///< [mz] Attenuation in water (dB/wavelength) vs. depth
        const std::vector<Ram::Real>& zcb,    ///< Depth points for cb
        const std::vector<Ram::Real>& cb,     ///< [mz] Sound speed in sediment (m/s) vs. depth
        const std::vector<Ram::Real>& zrhob,  ///< Depth points for rhob
        const std::vector<Ram::Real>& rhob,   ///< [mz] Density in sediment (rel. water) vs. depth
        const std::vector<Ram::Real>& zattn,  ///< Depth points for attn
        const std::vector<Ram::Real>& attn    ///< [mz] Attenuation in sediment (dB/wavelength) vs. depth
    );

    /// Interpolate a profile vs. z (depth)
    /**
     * The input data are given at user-specified depth points.
     * The output data are tabulated on the fixed depth grid common to
     * all depth profiles: nz_ samples starting at 0 with spacing dz_.
     */
    void InterpVsZ(
        const std::vector<Ram::Real>& zin,  ///< Input depth points
        const std::vector<Ram::Real>& pin,  ///< INput values vs. depth
        std::vector<Ram::Real>& pout        ///< [mz] Interpolated values vs. depth
    );
    
    /// Compute the tridiagonal matrices.
    /** This should be called after SetProfiles changes any profile.
     * All inputs and outputs are RamRunner member data.
     */
    void ComputeMatrices();
    
    /// Set the sizes of the various vectors
    void setSizes(
        int mr,     ///< Maximum number of range steps
        int mz,     ///< Maximum number of depth steps
        int mp      ///< Maximum number of terms in rational approximation
    );
    
    /// Propagate the field to the next range.
    /** @return the range to which the field has been propagated */
    Ram::Real PropagateToNextRange();
    
    /// Get the current range [m]
    Ram::Real range() const
    {
        return r_;
    }
    
    /// Get the relative pressure field vs. depth at the current range
    /** 
     * The phase does not include the fast-changing phase exp(+i*k*r_).
     * The amplitude does not include the spherical spreading factor.
     *
     * The output dz specifies the depths at which pressures are returned:
     * pressure[j] is at depth j*dz.
     * The maximum depth is dz*(pressure.size()-1).
     * They are determined by the input variables auplt, alplt, dz, and ndz.
     * 
     * @pre  nzplt_ / ndz_ <= pressure.size()
     */
    void getPressureField(
        std::vector<Ram::Complex>& pressure, ///< OUT: Pressure vs. depth
        Ram::Real& dz           ///< OUT: Increment in depth
    ) const;
    
    /// Get the complex pressure at the current range, at a given depth
    /**
     * The result linearly interpolated between the 2 nearest z samples.
     *
     * The phase does not include the fast-changing exp(+i*k*r_).
     * The amplitude DOES include the spherical spreading factor.
     */
    Ram::Complex getPressureAtDepth(
        Ram::Real depth    ///< Depth [m] below the water surface
    ) const;
    
    /// Write the transmission loss files in the original RAM format
    void writeLossFiles();
    
    /// Is there more to come?
    bool moreToCome() const
    {
        return r_ < rmax_;
    }
    
    
private:
    // ------------------------------------------------------------------------
    // Object State
    int mr_;        ///< Maximum number of range steps
    int mz_;        ///< Maximum number of depth steps
    int nz_;        ///< Number of depth steps ((amax-amin)/dz-0.5)
    int mp_;        ///< Maximum number of terms in rational approximation
    int np_;        ///< Number of terms in rational approximation
    int ns_;        ///< Number of stability constraints (1 or 2)
    int mdr_;       ///< Current range index (set to 0)
    int ndr_;       ///< Range decimation factor for tl.grid (1=no decimation)
    int ndz_;       ///< Depth decimation factor for tl.grid (1=no decimation)
    int iz_;        ///< Depth index of bottom
    int nzplt_;     ///< Maximum depth index in output to tl.grid (zmplt/dz-0.5)
    int lz_;        ///< Number of depths output to tl.grid (~nzplt/ndz)
    int ib_;        ///< Current index in bathymetry array
    int ir_;        ///< Depth index of receiver
    Ram::Real dir_;     ///< Depth of receiver below the (iz+ir)'th depth step
    Ram::Real dr_;      ///< Range step (m)
    Ram::Real dz_;      ///< Depth grid spacing (m)
    Ram::Real pi_;      ///< Math constant 3.14... (4.0*atan(1.0))
    Ram::Real eta_;     ///< Attenuation conversion factor: 1.0/(40.0*pi*alog10(exp(1.0)))
    Ram::Real eps_;     ///< A very small number (1.0e-20)
    Ram::Real omega_;   ///< Angular frequency (radians per second): 2.0*pi*freq
    Ram::Real rmax_;    ///< Maximum range (m)
    Ram::Real c0_;      ///< Reference sound speed (m/s)
    Ram::Real k0_;      ///< Reference wave number (m^-1): omega/c0
    Ram::Real r_;       ///< Current range (m) (set to dr)
    Ram::Real rp_;      ///< Range of profile update (m)
    Ram::Real rs_;      ///< Maximum range of stability constraints (m)
    std::vector<Ram::Real> rb_;     ///< [mr] Range of bathymetry point (m) vs. range
    std::vector<Ram::Real> zb_;     ///< [mr] Depth of bathymetry point (m) vs. range
    std::vector<Ram::Real> cw_;     ///< [mz] Sound speed in water column (m/s) vs. depth
    std::vector<Ram::Real> rhow_;   ///< [mz] Density in water (rel. water) vs. depth
    std::vector<Ram::Real> attnw_;   ///< [mz] Attenuation in water (dB/wavelength) vs. depth
    std::vector<Ram::Real> cb_;     ///< [mz] Sound speed in sediment (m/s) vs. depth
    std::vector<Ram::Real> rhob_;   ///< [mz] Density in sediment (rel. water) vs. depth
    std::vector<Ram::Real> attn_;   ///< [mz] Attenuation in sediment (dB/wavelength) vs. depth
    std::vector<Ram::Real> alpw_;   ///< [mz] Energy conversion factor sqrt(cw/c0) vs. depth
    std::vector<Ram::Real> alpb_;
        ///< [mz] OUT: Energy conversion factor sqrt(rhob*cb/c0) vs. depth
    std::vector<Ram::Complex> ksq_;
        ///< [mz] Complex relative wave number squared (ksqa or ksqb) vs. depth
    std::vector<Ram::Complex> ksqw_;
        ///< [mz] Relative wave number squared (kw^2-k0^2) in water vs. depth
    std::vector<Ram::Complex> ksqb_; 
        ///< [mz] Complex relative wave number squared (kb^2-k0^2) in bottom vs. depth
    std::vector<Ram::Real> f1_;     ///< [mz] Matrix input: rho/alpha vs. depth
    std::vector<Ram::Real> f2_;     ///< [mz] Matrix input: rho vs. depth
    std::vector<Ram::Real> f3_;     ///< [mz] Matrix input: alpha vs. depth
    std::vector<Ram::Complex> u_;   ///< [mz] The answer: Complex pressure/ksq vs. depth
    std::vector<Ram::Complex> v_;
        ///< [mz] Right-hand side of matrix equation; vs. depth (scratch used by solve)
    std::vector<Ram::Complex> r1_;  ///< [mp,mz] LHS matrix elements below the diagonal?
    std::vector<Ram::Complex> r2_;  ///< [mp,mz] LHS matrix elements on the diagonal?
    std::vector<Ram::Complex> r3_;  ///< [mp,mz] LHS matrix elements above the diagonal?
    std::vector<Ram::Complex> s1_;  ///< [mp,mz] RHS matrix elements below the diagonal?
    std::vector<Ram::Complex> s2_;  ///< [mp,mz] RHS matrix elements on the diagonal?
    std::vector<Ram::Complex> s3_;  ///< [mp,mz] RHS matrix elements above the diagonal?
    std::vector<Ram::Complex> pd1_; ///< [mp] Numerator Pade coefficients?
    std::vector<Ram::Complex> pd2_; ///< [mp] Denominator Pade coefficients?
    std::vector<float> tlg_;    ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    std::vector<float> preal_;    ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    std::vector<float> pimag_;    ///< [mz] Transmission loss in dB, vs. depth (scratch used by outpt)
    bool filesAreOpen_;     ///< Were the Fortran I/O files opened at construction time?
    int isrc_type_;      ///< source: 0 for self-starter, 1 for gaussian tapered point source
    Ram::Real x0_;       ///< source range (m)
    Ram::Real x1_;       ///< starting field range (m)
    Ram::Real theta0_;   ///< Beam angle (deg)
    Ram::Real dtheta_;   ///< Beam width (deg)
};

#endif
/* DON'T ADD STUFF AFTER THIS #endif */
