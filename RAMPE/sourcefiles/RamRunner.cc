/**
 * @file
 * @brief RamRunner: C++ class to encapsulate a run of RAM
 */

#include "RamRunner.hh"

#include "ramlib.hh" // C interfaces for the Fortran subroutines

#include <cassert>

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
RamRunner::RamRunner(
    const std::string& inFile,      ///< Name of the input text file
    const std::string& lineFile,    ///< Name of the output line text file
    const std::string& gridFile     ///< Name of the output grid binary file
)
:   mr_( 0 ),
    mz_( 0 ),
    nz_( 0 ),
    mp_( 0 ),
    np_( 0 ),
    ns_( 0 ),
    mdr_( 0 ),
    ndr_( 0 ),
    ndz_( 0 ),
    iz_( 0 ),
    nzplt_( 0 ),
    lz_( 0 ),
    ib_( 0 ),
    ir_( 0 ),
    dir_( 0.0 ),
    dr_( 0.0 ),
    dz_( 0.0 ),
    pi_( 0.0 ),
    eta_( 0.0 ),
    eps_( 0.0 ),
    omega_( 0.0 ),
    rmax_( 0.0 ),
    c0_( 0.0 ),
    k0_( 0.0 ),
    r_( 0.0 ),
    rp_( 0.0 ),
    rs_( 0.0 ),
    filesAreOpen_( true ),
    isrc_type_( 0 ),
    x0_( 0.0 ),
    x1_( 0.0 ),
    theta0_( 0.0 ),
    dtheta_( 0.0 )
{
    int clen = inFile.size();
    int llen = lineFile.size();
    int glen = gridFile.size();
    int mr, mz, mp;
    
    // Determine the sizes of the various vectors
    // by reading part of the input file
    sizes( mr, mz, mp, inFile.c_str(), clen );
    
    // Set the sizes of the various vectors
    setSizes( mr, mz, mp );
    
    //Initialize the parameters, acoustic field, and matrices.
    setup( mr_, mz_, nz_, mp_, np_, ns_, mdr_, ndr_, ndz_, iz_, nzplt_, lz_, ib_, ir_,
          dir_,  dr_, dz_, pi_, eta_, eps_, omega_, rmax_, c0_, k0_, r_, rp_, rs_,
          &rb_[0], &zb_[0], &cw_[0],&rhow_[0], &attnw_[0], &cb_[0],  &rhob_[0],  &attn_[0],
          &alpw_[0], &alpb_[0], &ksq_[0], &ksqw_[0], &ksqb_[0], 
          &f1_[0], &f2_[0], &f3_[0], &u_[0], &v_[0], 
          &r1_[0], &r2_[0], &r3_[0], &s1_[0], &s2_[0], &s3_[0], 
          &pd1_[0], &pd2_[0], &tlg_[0], &preal_[0], &pimag_[0],
          inFile.c_str(), clen, lineFile.c_str(), llen, gridFile.c_str(), glen,
          isrc_type_, x0_, x1_, theta0_, dtheta_);
    
}   // end RamRunner constructor taking file names


/// Full Constructor taking inputs from constructor arguments, not a file
/**
 * This constructor does NOT read any input files, nor does it open
 * any output files.
 * 
 * The input parameters are in the same order as in the RAM input file.
 */
RamRunner::RamRunner(
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
    int isrc_type,      ///< IN: source: 0 for self-starter, 1 for gaussian tapered point source
    Ram::Real x0,       ///< source range (m)
    Ram::Real x1,       ///< starting field range (m)
    Ram::Real theta0,   ///< beam angle (deg)
    Ram::Real dtheta    ///< beam width (deg)
)
:   mr_( 0 ),       // set in body
    mz_( 0 ),       // set in body
    nz_( zmax/dz - 0.5 ),
    mp_( np ),
    np_( np ),
    ns_( ns ),
    mdr_( 1 ),
    ndr_( ndr ),
    ndz_( ndz ),
    iz_( 0 ),       // initialized in body
    nzplt_( int( zmplt/dz-0.5 ) ),
    lz_( 0 ),       // set in body
    ib_( 1 ),
    ir_( 0 ),       // initialized in body
    dir_( 0.0 ),    // set in body
    dr_( dr ),
    dz_( dz ),
    pi_( 4.0*atan( 1.0 ) ),
    eta_( 1.0/(40.0*pi_*log10(exp(1.0))) ),
    eps_( 1.0e-20 ),
    omega_( 2.0*pi_*freq ),
    rmax_( rmax ),
    c0_( c0 ),
    k0_( omega_/c0 ),
    r_( dr ),
    rp_( 2.0*rmax ), // Profile changes come from caller, not from Fortran
    rs_( rs ),
    filesAreOpen_( false ),
    isrc_type_( isrc_type ),
    x0_( x0 ),
    x1_( x1 ),
    theta0_( theta0 ),
    dtheta_( dtheta )
{
    // Check inputs for sanity
    assert( freq > 0.0 );
    assert( zs >= 0.0 && zs <= zmax );
    assert( zr >= 0.0 && zr <= zmax );
    assert( rmax > 0.0 );
    assert( dr > 0.0 );
    assert( ndr > 0 );
    assert( zmax > 0.0 );
    assert( dz > 0.0 );
    ndz = std::max( 1, ndz );
    assert( zmplt >= 0.0 && zmplt <= zmax );
    assert( c0 > 0.0 );
    assert( np > 0 && np <= 12 );
    assert( ns >= 0 && ns <= 2 );
    assert( rb.size() == zb.size() );
    assert ( zb.size() > 0 );
    assert( cw.size() == zcw.size() );
    assert ( cw.size() > 0 );
    assert( rhow.size() == zrhow.size() );
    assert ( rhow.size() > 0 );
    assert( attnw.size() == zattnw.size() );
    assert ( attnw.size() > 0 );
    assert( cb.size() == zcb.size() );
    assert ( cb.size() > 0 );
    assert( rhob.size() == zrhob.size() );
    assert ( rhob.size() > 0 );
    assert( attn.size() == zattn.size() );
    assert ( attn.size() > 0 );
    assert( isrc_type >= 0 && isrc_type <= 1 );
    if ( isrc_type == 1 ) {
        assert( x1 >= x0 );
        assert( theta0 > -90.0 && theta0 < 90.0);
        assert( dtheta > 0.0 && dtheta < 90.0 );
    }
    
    // How many depth steps do we need?
    nz_ = zmax/dz - 0.5;
    mz_ = nz_ + 2;
    
    // Set the sizes of the various vectors
    setSizes( zb.size()+1, mz_, np );
    
    // Copy the bathymetry arrays
    rb_ = rb;
    zb_ = zb;
    rb_.push_back( 2.0*rmax );
    zb_.push_back( zb.back() );
    assert( rb_.size() == zb_.size() );
    assert( rb_.size() == mr_ );
    
    // Receiver depth as step number and fraction
    Ram::Real ri = 1.0 + zr/dz;
    ir_ = int( ri );    // Fortran-style 1-base index
    ir_ = std::max( 1, std::min( nz_, ir_ ) );
    dir_ = ri - ir_;
    
    // Number of output depths may not exceed number of depth samples
    nzplt_ = std::min( nzplt_, nz_ );
    
    // Initialize the parameters that will change from step to step
    Ram::Real z = zb[0];
    iz_ = 1.0 + z/dz;   // 1-based Fortran index
    iz_ = std::max( 2, iz_ );
    iz_ = std::min( nz_, iz_ );
    if ( rs < dr )
        rs_ = 2.0*rmax;
    
    // Number of depths output to tl.grid (~nzplt/ndz)
    lz_ = 0;
    for ( int i = ndz; i <= nzplt_; i += ndz )
        ++lz_;
    
    // Set up all of the profiles versus z
    SetProfiles( zcw, cw, zrhow, rhow, zattnw, attnw, zcb, cb, zrhob, rhob, zattn, attn );
    
    // Compute the starting field
    if ( isrc_type == 0 ) {
        // Use the "self starter" (point source)
        selfs( mz_, nz_, mp_, np_, ns_, iz_, zs, dr_, dz_, pi_, c0_, k0_,
            &rhow_[0], &rhob_[0], &alpw_[0], &alpb_[0], &ksq_[0], &ksqw_[0], &ksqb_[0],
            &f1_[0], &f2_[0], &f3_[0], &u_[0], &v_[0],
            &r1_[0], &r2_[0], &r3_[0], &s1_[0], &s2_[0], &s3_[0],
            &pd1_[0], &pd2_[0] );
    } else { // isrc_type == 1
        
        // Gaussian tapered point source for starting field
        pt_src_gtaper( mz_, nz_, dz_, freq, c0_, zs, x0, x1, theta0,
            dtheta, &u_[0] );
        matrc( mz_, nz_, mp_, np, iz_, iz_, dz_, k0_,
            &rhow_[0] ,&rhob_[0], &alpw_[0], &alpb_[0], &ksq_[0],
            &ksqw_[0], &ksqb_[0], &f1_[0], &f2_[0], &f3_[0], &r1_[0],
            &r2_[0], &r3_[0], &s1_[0], &s2_[0], &s3_[0], &pd1_[0], &pd2_[0]);
    }
    
    // Compute the Pade coefficients
    epade( mp_, np_, ns_, 1, k0_, c0_, dr_, &pd1_[0], &pd2_[0] );
    
    // Compute the propagation matrices.
    ComputeMatrices();
   
}   // end RamRunner constructor taking parameter values


/// Compute the tridiagonal matrices.
/** This should be called after SetProfiles changes any profile.
 * All inputs and outputs are RamRunner member data.
 */
void RamRunner::ComputeMatrices()
{
    matrc( mz_, nz_, mp_, np_, iz_, iz_, dz_, k0_,
        &rhow_[0], &rhob_[0], &alpw_[0], &alpb_[0], &ksq_[0],
        &ksqw_[0], &ksqb_[0], &f1_[0], &f2_[0], &f3_[0],
        &r1_[0], &r2_[0], &r3_[0], &s1_[0], &s2_[0], &s3_[0], 
        &pd1_[0], &pd2_[0] );
}


/// Set the sizes of the various vectors
void RamRunner::setSizes(
    int mr,     ///< Maximum number of range steps
    int mz,     ///< Maximum number of depth steps
    int mp      ///< Maximum number of terms in rational approximation
)
{
    assert( mr >= 0 );  // Can be zero; use push_back to fill bathymetry
    assert( mz > 0 );
    assert( mp > 0 && mp <= 12 );
    
    mr_ = mr;
    mz_ = mz;
    mp_ = mp;
    rb_.resize( mr_ );
    zb_.resize( mr_ );
    cw_.resize( mz_ );
    rhow_.resize( mz_ );
    attnw_.resize( mz_ );
    cb_.resize( mz_ );
    rhob_.resize( mz_ );
    attn_.resize( mz_ );
    alpw_.resize( mz_ );
    alpb_.resize( mz_ );
    ksq_.resize( mz_ );
    ksqw_.resize( mz_ );
    ksqb_.resize( mz_ );
    f1_.resize( mz_ );
    f2_.resize( mz_ );
    f3_.resize( mz_ );
    u_.resize( mz_ );
    v_.resize( mz_ );
    r1_.resize( mz_*mp_ );
    r2_.resize( mz_*mp_ );
    r3_.resize( mz_*mp_ );
    s1_.resize( mz_*mp_ );
    s2_.resize( mz_*mp_ );
    s3_.resize( mz_*mp_ );
    pd1_.resize( mp_ );
    pd2_.resize( mp_ );
    tlg_.resize( mz_ );
    preal_.resize( mz_ );
    pimag_.resize( mz_ );
}       // end RamRunner::setSizes


/// Destructor closes the files and deallocates memory
RamRunner::~RamRunner()
{
    if ( filesAreOpen_ )
        close_ram_files();
}   // end RamRunner destructor


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
void RamRunner::SetProfiles(
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
)
{
    // Interpolate the profiles onto the common depth grid
    InterpVsZ( zcw, cw, cw_ );
    InterpVsZ( zrhow, rhow, rhow_ );
    InterpVsZ( zattnw, attnw, attnw_ );
    InterpVsZ( zcb, cb, cb_ );
    InterpVsZ( zrhob, rhob, rhob_ );
    InterpVsZ( zattn, attn, attn_ );

    // Compute the derived profiles from the interpolated ones
    for ( int i = 0; i < nz_+2; ++i ) {
        // ksqw_[i]=(omega/cw_[i])**2-k0_**2
        Ram::Complex tempkw = (omega_/cw_[i])
                * Ram::Complex( 1.0, eta_*attnw_[i] );
        ksqw_[i]=(tempkw-k0_)*(tempkw+k0_);
        // ksqb_[i]=((omega/cb_[i])*(1.0+ci*eta*attn_[i]))**2-k0**2
        Ram::Complex tempkb = (omega_/cb_[i])
                * Ram::Complex( 1.0, eta_*attn_[i] );
        ksqb_[i] = (tempkb-k0_)*(tempkb+k0_);
        alpw_[i] = std::sqrt(rhow_[i]*cw_[i]/c0_);
        alpb_[i] = std::sqrt(rhob_[i]*cb_[i]/c0_);
    }
}       // end RamRunner::SetProfiles


/// Interpolate a profile vs. z (depth)
/**
 * The input data are given at user-specified depth points.
 * The output data are tabulated on the fixed depth grid common to
 * all depth profiles: nz_ samples starting at 0 with spacing dz_.
 */
void RamRunner::InterpVsZ(
    const std::vector<Ram::Real>& zin,  ///< Input depth points
    const std::vector<Ram::Real>& pin,  ///< Input values vs. depth
    std::vector<Ram::Real>& pout        ///< [mz] Interpolated values vs. depth
)
{
    assert( zin.size() == pin.size() );
    int nin = zin.size();
    assert( nin > 0 );
    assert( pout.size() > nz_ );    // Need one extra value
    
    // Output samples before the first input sample get the first value
    int i = zin[0]/dz_ + 0.5;
    i = std::max( 0, std::min( nz_, i ) );
    std::fill( pout.begin(), pout.begin() + i, pin[0] );
    
    // Output samples between input samples are interpolated linearly
    for ( int j = 1; j < nin; ++j ) {
        int iold = i;   // next sample to set
        i = zin[j]/dz_ + 0.5;
        i = std::max( iold, std::min( nz_, i ) );
        int nramp = i - iold;
        if ( nramp > 0 ) {
            Ram::Real dp = (pin[j] - pin[j-1])/nramp;
            for ( int k = 0; k < nramp; ++k )
                pout[iold+k] = pin[j-1] + k*dp;
        }
    }
    
    // Output samples after the last input sample get the last value
    std::fill( pout.begin() + i, pout.end(), pin.back() );
}       // end RamRunner::InterpVsZ

/// Propagate the field to the next output range.
/** @return the range to which the field has been propagated */
Ram::Real RamRunner::PropagateToNextRange()
{
    if ( mdr_ >= ndr_ )
        mdr_ = 0;
    for ( int mdr = mdr_; mdr < ndr_; ++mdr ) {
        // Update the tridiagonal matrix
        updat(mr_, mz_, nz_, mp_, np_, iz_, ib_, dr_, dz_, eta_, omega_, rmax_, c0_, k0_, r_, 
              rp_, rs_,&rb_[0], &zb_[0], &cw_[0], &rhow_[0], &attnw_[0],  &cb_[0], &rhob_[0], &attn_[0],
              &alpw_[0], &alpb_[0], &ksq_[0], &ksqw_[0], &ksqb_[0], &f1_[0], &f2_[0], &f3_[0], 
              &r1_[0], &r2_[0], &r3_[0], &s1_[0], &s2_[0], &s3_[0], &pd1_[0], &pd2_[0]);
        
        //Solve the tridiagonal matrix equation
        solve(mz_, nz_, mp_, np_, iz_, &u_[0], &v_[0], &r1_[0], &r2_[0], &r3_[0],
              &s1_[0], &s2_[0], &s3_[0]);
        r_ += dr_;
    }
    mdr_ = ndr_;    // to force outpt to write its output
    return r_;
}   // end RamRunner::propagateToNextRange


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
void RamRunner::getPressureField(
    std::vector<Ram::Complex>& pressure, ///< OUT: Pressure vs. depth
    Ram::Real& dz           ///< OUT: Increment in depth
) const
{
    // Copy pressure from internal data
    int lz = nzplt_ / ndz_;
    assert( lz <= pressure.size() );
    for ( int j = 0; j < lz; ++j ) {
        int i = j*ndz_;
        pressure[j] = u_[i] * f3_[i];
    }
    
    // Fill it out with zeros
    size_t nz = pressure.size();
    for ( int j = lz; j < nz; ++j )
        pressure[j] = Ram::Complex();
    
    // Set the output depth increment
    dz = dz_*ndz_;
}   // end RamRunner::getPressureField


/// Get the complex pressure at the current range, at a given depth
/**
 * The result linearly interpolated between the 2 nearest z samples.
 *
 * The phase does not include the fast-changing exp(+i*k*r_).
 * The amplitude DOES include the spherical spreading factor.
 */
Ram::Complex RamRunner::getPressureAtDepth(
    Ram::Real depth    ///< Depth (m) below the water surface
) const
{
    Ram::Real ri = depth/dz_;
    int irz = int( ri );
    irz = std::max( 0, std::min( mz_-2, irz ) );
    Ram::Real dir = ri - Ram::Real( irz );
    Ram::Complex p = (1.0-dir)*f3_[irz]*u_[irz] + dir*f3_[irz+1]*u_[irz+1];
    return p / std::sqrt( r_ );
}   // end RamRunner::getPressureAtDepth


/// Write the transmission loss files in the original RAM format
void RamRunner::writeLossFiles()
{
    outpt( mz_, mdr_, ndr_, ndz_, iz_, nzplt_, lz_, ir_, dir_, eps_, r_,
           &f3_[0], &u_[0], &tlg_[0] , &preal_[0], &pimag_[0]);
}
