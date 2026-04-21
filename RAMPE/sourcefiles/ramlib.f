c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 1.5, 13-Sep-00 *****
c     ******************************************************************
c     
c     This code was developed by Michael D. Collins at the Naval
c     Research Laboratory in Washington, DC. It solves range-dependent 
c     ocean acoustics problems with the split-step Pade algorithm 
c     [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 (1993)]. A 
c     user's guide and updates of the code are available via anonymous 
c     ftp from ram.nrl.navy.mil.
c
c     Library for RAM 1.5
c
c     This file contains all of Michael Collins's RAM 1.5
c     except for the main program, packaged as a module called ramlib.
c     The functions are now portably C-callable, with C-compatible
c     types. Also, most computation is now double precision, but
c     easily changed back to single precision for comparison.
c     Those changes were made in 2012 by Robert P. Goddard.
c
c     Version 1.5 contains a correction to a bug in the dimension of
c     quantities passed to subroutines fndrt and guerre that Laurie
c     Fialkowski noticed. 
c
c     Version 1.4 contains a correction to a minor bug in subroutine
c     guerre that Dave King noticed (amp1 and amp2 were declared
c     twice) and a few other minor improvements. 
c
c     Version 1.3 contains a new root-finding subroutine.
c
c     Version 1.2 contains a minor modification. The output to tl.grid
c     is no longer zeroed out along the ocean bottom. This was done in
c     previous versions so that the ocean bottom would be highlighted
c     in graphical displays. The graphics codes ramclr, ramctr, and
c     ramcc read in the bathymetry from ram.in and plot the ocean
c     bottom directly. 
c
c     Version 1.1 contains two improvements:
c
c     (1) An improved self starter. Stability is improved by using the 
c     factor (1-X)**2 instead of (1+X)**2 to smooth the delta function. 
c     The factor (1+X)**2 is nearly singular for some problems involving
c     deep water and/or weak attenuation. Numerical problems associated 
c     with this singularity were detected by Eddie Scheer of Woods Hole 
c     Oceanographic Institute. 
c
c     (2) Elimination of underflow problems. A very small number is 
c     added to the solution in subroutine solve to prevent underflow,
c     which can adversely affect run time on some computers. This
c     improvement was suggested by Ed McDonald of the SACLANT Undersea
c     Research Centre. 
c

      module ramlib
c       These 'kind' type parameters define the types used in the interfaces.
c       To choose single precision, use c_float and c_float_complex here.
c       To choose double precision, use c_double and c_double_complex here.
        use, intrinsic :: iso_c_binding
        implicit real(c_double) (a-h, o-z)
        integer, parameter :: ram_int = c_int
        integer, parameter :: ram_real = c_double
        integer, parameter :: ram_complex = c_double_complex
c       Note: implicit is deprecated, other than 'implicit none'.
      contains
c
c     The module body
c
c     ----------------------------------------------------------------
c
c     helper subroutine to convert c string to FORTRAN string
c
      subroutine c2fstring (cstring, ilen, fstring)
      character(c_char), intent(in) :: cstring(ilen)
      integer(ram_int) :: ilen
      character*128 fstring
      
      do 2, ii=1,ilen
        fstring(ii:ii) = cstring(ii)
    2 continue
      do 3, ii=ilen+1,128
        fstring(ii:ii) = ' '
    3 continue

      return
      end subroutine c2fstring
      
c     ----------------------------------------------------------------
c
c     Determine array sizes by reading some of file ram.in.
c     The file is closed before this subroutine returns.
c
      subroutine sizes( mr, mz, mp, cfname, len )
     > bind(c)
      integer(ram_int) :: mr,mz,mp
      integer(ram_int) :: len
      character(c_char), intent(in) :: cfname(len)
c
      character*128 :: fname
c
      call c2fstring(cfname,len,fname)
c
c     Open the input and output files on I/O Units 1-3.
c     Unit 1, 'ram.in', is the input: a text file as specified in the
c     User's Guide. open(unit=1,status='old',file='ram.in')

      open(unit=1,status='old',file=fname)
      read(1,*)
      read(1,*)isrc_type, x0, x1, theta0, dtheta
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs
c
c     How many ranges do we need in the bathymetry array?
      i=1
    1 continue
        read(1,*)rb,zb
        if(rb.lt.0.0)go to 2
        i=i+1
        go to 1
    2 mr = i
c
c     How many depth steps do we need?
      nz=zmax/dz-0.5
      mz = nz + 2
c
c     How many Pade terms do we need?
      mp = np
c
c     Close the file. Subroutine setup will re-open it.
      close( 1 )
      
      return
      end subroutine sizes
      
c     ----------------------------------------------------------------
c
c     Initialize the parameters, acoustic field, and matrices.
c
      subroutine setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,
     >   ir,dir,dr,dz,pi,eta,eps,omega,rmax,c0,k0,r,rp,rs,rb,zb,cw,
     >   rhow,attnw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,
     >   u,v,r1,r2,r3,s1,s2,s3,pd1,pd2,tlg,preal,pimag,
     >   incstr, ilen, linecstr, llen, gridcstr, glen,
     >   isrc_type, x0, x1, theta0, dtheta)
     > bind(c)
      integer(ram_int) :: 
     >   mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,isrc_type
      complex(ram_complex) ::
     >   u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),ksqw(mz),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
      real(ram_real) :: 
     >   k0,rb(mr),zb(mr),cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),
     >   alpb(mz),f1(mz),f2(mz),f3(mz),rhow(mz),attnw(mz),
     >   dir,dr,dz,pi,eta,eps,omega,rmax,c0,r,rp,rs,
     >   x0, x1, theta0, dtheta
      real(c_float) ::
     >   tlg(mz), preal(mz), pimag(mz)
      integer(ram_int) :: ilen, llen, glen
      character(c_char), intent(in) :: 
     >   incstr(ilen), linecstr(llen), gridcstr(glen)
      
c     Locals
      character*128 :: infile, linefile, gridfile
c
c     Translate C strings to Fortran strings for filenames
      call c2fstring(incstr,ilen,infile)
      call c2fstring(linecstr,llen,linefile)
      call c2fstring(gridcstr,glen,gridfile)
c
c     Open the input and output files on I/O Units 1-3.
c     Unit 1, 'ram.in', is the input: a text file as specified in
c      the User's Guide.
      open(unit=1,status='old',file=infile)
c
c     Unit 2, 'tl.line', is a text output file that will get the dB
c     transmission loss versus range (2 columns) at the target depth.
      open(unit=2,status='unknown',file=linefile)
c
c     Unit 3, 'tl.grid', is a Fortran unformatted file that will
c      get the dB transmission loss on a grid of ranges and depths.
      open(unit=3,status='unknown',file=gridfile,form='unformatted')
c
      read(1,*)
      read(1,*) isrc_type, x0, x1, theta0, dtheta
      read(1,*)freq,zs,zr
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs
c
      i=1
    1 continue
        read(1,*)rb(i),zb(i)
        if(rb(i).lt.0.0)go to 2
        i=i+1
        go to 1
    2 rb(i)=2.0*rmax
      zb(i)=zb(i-1)
c
      pi=4.0*atan(1.0)
      eta=1.0/(40.0*pi*log10(exp(1.0)))
      eps=1.0e-20
      ib=1

      if (isrc_type .eq. 0) then
         mdr=1
         r=dr
      else if (isrc_type .eq. 1) then
C  ( x1/dr+1 ?? - need to check)
         mdr= nint(x1/dr)
         r=x1
      else
         write(*,*)' '
         write(*,*)' ERROR in setup: isrc_type must b 0 or 1'
         write(*,*)' '
         stop
      end if

      omega=2.0*pi*freq
      ri=1.0+zr/dz
      ir=int(ri)
      dir=ri-ir
      k0=omega/c0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5
      z=zb(1)
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      if(rs.lt.dr)rs=2.0*rmax
c
      if(nz+2.gt.mz)then
        write(*,*)'   Need to increase parameter mz to ',nz+2
        stop
      end if
      if(np.gt.mp)then
        write(*,*)'   Need to increase parameter mp to ',np
        stop
      end if
      if(i.gt.mr)then
        write(*,*)'   Need to increase parameter mr to ',i
        stop
      end if
c
      do 3 j=1,mp
        r3(1,j)=0.0
        r1(nz+2,j)=0.0
    3 continue
      do 4 i=1,nz+2
      u(i)=0.0
      v(i)=0.0
    4 continue
      lz=0
      do 5 i=ndz,nzplt,ndz
      lz=lz+1
    5 continue
      write(3)lz
c
c     The initial profiles and starting field.
c
      call profl(mz,nz,dz,eta,omega,rmax,c0,k0,rp,cw,rhow,attnw,
     >   cb,rhob,attn,alpw,alpb,ksqw,ksqb)
c
      if (isrc_type .eq. 0) then
c
c use "self starter"
         x0 = 0.0
         call selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,pi,c0,k0,rhow,rhob,
     >      alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
     >      pd1,pd2)
      else if (isrc_type .eq. 1) then
c
c use "a Gaussian tapered point source"
c 
   
         call pt_src_gtaper(mz, nz, dz, freq, c0, zs, x0, x1, 
     >                    theta0, dtheta, u)
         call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >       alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      else
         write(*,*)' '
         write(*,*)'   Unrecognized source type'
         write(*,*)'   0: self starter'
         write(*,*)'   1: Gaussian tapered point source'
         write(*,*)' '
         stop
      end if

      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg,
     >     preal,pimag)
c
c     The propagation matrices.
c
      call epade(mp,np,ns,1,k0,c0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
      return
      end subroutine setup

c     ----------------------------------------------------------------
c
c     Set up the profiles.
c
      subroutine profl(mz,nz,dz,eta,omega,rmax,c0,k0,rp,cw,rhow,attnw,
     >   cb,rhob,attn,alpw,alpb,ksqw,ksqb)
     > bind(c)
      integer(ram_int) :: mz,nz
      complex(ram_complex) :: ksqb(mz),ksqw(mz)
      real(ram_real) :: dz,eta,omega,rmax,c0,rp,
     >   k0,cw(mz),rhow(mz),attnw(mz),cb(mz),rhob(mz),attn(mz),
     >   alpw(mz),alpb(mz)
c     Locals
      complex(ram_complex) :: tempkb,tempkw

      call zread(mz,nz,dz,cw)
      call zread(mz,nz,dz,rhow)
	  call zread(mz,nz,dz,attnw)
      call zread(mz,nz,dz,cb)
      call zread(mz,nz,dz,rhob)
      call zread(mz,nz,dz,attn)
      rp=2.0*rmax
      read(1,*,end=1)rp
c
    1 do 2 i=1,nz+2
c       ksqw(i)=(omega/cw(i))**2-k0**2
        tempkw = (omega/cw(i))* cmplx(1.0,eta*attnw(i),ram_complex)
        ksqw(i)=(tempkw-k0)*(tempkw+k0)
c       ksqb(i)=((omega/cb(i))*(1.0+ci*eta*attn(i)))**2-k0**2
        tempkb = (omega/cb(i)) * cmplx(1.0,eta*attn(i),ram_complex)
        ksqb(i)=(tempkb-k0)*(tempkb+k0)
        alpw(i)=sqrt(rhow(i)*cw(i)/c0)
        alpb(i)=sqrt(rhob(i)*cb(i)/c0)

    2 continue
c
      return
      end subroutine profl

c     ----------------------------------------------------------------
c
c     Profile reader and interpolator.
c
      subroutine zread(mz,nz,dz,prof)
     > bind(c)
      integer(ram_int) :: mz,nz
      real(ram_real) :: dz, prof(mz), dprof
c
      do 1 i=1,nz+2
        prof(i)=-1.0
    1 continue
      read(1,*)zi,profi
      prof(1)=profi
      i=1.5+zi/dz
      prof(i)=profi
      iold=i
    2 continue
        read(1,*)zi,profi
        if(zi.lt.0.0)go to 3
        i=1.5+zi/dz
        if(i.eq.iold)i=i+1
        prof(i)=profi
        iold=i
        go to 2
    3 prof(nz+2)=prof(i)
      i=1
      j=1
    4 continue
        i=i+1
        if(prof(i).lt.0.0)go to 4
        if(i-j.eq.1)go to 6
c        dprof = (prof(i)-prof(j))/(i-j)
        do 5 k=j+1,i-1
c          prof(k)=prof(j) + (k-j)*dprof
        prof(k)=prof(j)+float(k-j)*(prof(i)-prof(j))/float(i-j)

    5   continue
    6   j=i
        if(j.lt.nz+2)go to 4
c
      return
      end subroutine zread

c     ----------------------------------------------------------------
c
c     The tridiagonal matrices.
c
      subroutine matrc(mz,nz,mp,np,iz,jz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
     > bind(c)
      integer(ram_int) :: mz,nz,mp,np,iz,jz
      complex(ram_complex) :: 
     >   d1,d2,d3,rfact,ksq(mz),ksqw(mz),ksqb(mz),
     >   r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),
     >   pd1(mp),pd2(mp)
      real(ram_real) :: 
     >   dz,k0,rhow(mz),rhob(mz),f1(mz),f2(mz),f3(mz),
     >   alpw(mz),alpb(mz)
c
      a1=k0**2/6.0
      a2=2.0*k0**2/3.0
      a3=k0**2/6.0
      cfact=0.5/dz**2
      dfact = 12.0
      dfact=1.0/dfact
c
c     New matrices when iz.eq.jz.
      if(iz.eq.jz)then
      i1=2
      i2=nz+1
      do 1 i=1,iz
        f1(i)=rhow(i)/alpw(i)
        f2(i)=1.0/rhow(i)
        f3(i)=alpw(i)
        ksq(i)=ksqw(i)


    1 continue
      do 2 i=iz+1,nz+2
        f1(i)=rhob(i)/alpb(i)
        f2(i)=1.0/rhob(i)
        f3(i)=alpb(i)
        ksq(i)=ksqb(i)
    2 continue
      end if
c
c     Updated matrices when iz.ne.jz.
c
      if(iz.gt.jz)then
      i1=jz
      i2=iz+1
      do 3 i=jz+1,iz
        f1(i)=rhow(i)/alpw(i)
        f2(i)=1.0/rhow(i)
        f3(i)=alpw(i)
        ksq(i)=ksqw(i)
    3 continue
      end if
c
      if(iz.lt.jz)then
        i1=iz
        i2=jz+1
        do 4 i=iz+1,jz
            f1(i)=rhob(i)/alpb(i)
            f2(i)=1.0/rhob(i)
            f3(i)=alpb(i)
            ksq(i)=ksqb(i)

    4   continue
      end if
c
      do 6 i=i1,i2
c
c       Discretization by Galerkin's method.
c
        c1=cfact*f1(i)*(f2(i-1)+f2(i))*f3(i-1)
        c2=-cfact*f1(i)*(f2(i-1)+2.0*f2(i)+f2(i+1))*f3(i)
        c3=cfact*f1(i)*(f2(i)+f2(i+1))*f3(i+1)
        d1=c1+dfact*(ksq(i-1)+ksq(i))
        d2=c2+dfact*(ksq(i-1)+6.0*ksq(i)+ksq(i+1))
        d3=c3+dfact*(ksq(i)+ksq(i+1))
c
        do 5 j=1,np
            r1(i,j)=a1+pd2(j)*d1
            r2(i,j)=a2+pd2(j)*d2
            r3(i,j)=a3+pd2(j)*d3
            s1(i,j)=a1+pd1(j)*d1
            s2(i,j)=a2+pd1(j)*d2
            s3(i,j)=a3+pd1(j)*d3
    5   continue
    6 continue
c
c     The matrix decomposition.
c
      do 9 j=1,np
        do 7 i=i1,iz
            rfact=1.0/(r2(i,j)-r1(i,j)*r3(i-1,j))
             r1(i,j)=r1(i,j)*rfact
            r3(i,j)=r3(i,j)*rfact
            s1(i,j)=s1(i,j)*rfact
            s2(i,j)=s2(i,j)*rfact
            s3(i,j)=s3(i,j)*rfact
    7   continue
c
        do 8 i=i2,iz+2,-1
            rfact=1.0/(r2(i,j)-r3(i,j)*r1(i+1,j))
            r1(i,j)=r1(i,j)*rfact
            r3(i,j)=r3(i,j)*rfact
            s1(i,j)=s1(i,j)*rfact
            s2(i,j)=s2(i,j)*rfact
            s3(i,j)=s3(i,j)*rfact
    8   continue
c
        r2(iz+1,j)=r2(iz+1,j)-r1(iz+1,j)*r3(iz,j)
        r2(iz+1,j)=r2(iz+1,j)-r3(iz+1,j)*r1(iz+2,j)
        r2(iz+1,j)=1.0/r2(iz+1,j)
c
    9 continue
c
      return
      end subroutine matrc

c     ----------------------------------------------------------------
c
c     The tridiagonal solver.
c
      subroutine solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
     > bind(c)
      integer(ram_int) :: mz,nz,mp,np,iz
      complex(ram_complex) ::
     >   u(mz),v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp)
      
      eps=1.0e-30
c
c     The loop through Pade coefficients
      do 6 j=1,np
c
c     The right side.
c
      do 1 i=2,nz+1
      v(i)=s1(i,j)*u(i-1)+s2(i,j)*u(i)+s3(i,j)*u(i+1)+eps
    1 continue
c
c     The elimination steps.
c
      do 2 i=3,iz
      v(i)=v(i)-r1(i,j)*v(i-1)+eps
    2 continue
      do 3 i=nz,iz+2,-1
      v(i)=v(i)-r3(i,j)*v(i+1)+eps
    3 continue
c
      u(iz+1)=(v(iz+1)-r1(iz+1,j)*v(iz)-r3(iz+1,j)*v(iz+2))*
     >   r2(iz+1,j)+eps
c
c     The back substitution steps.
c
      do 4 i=iz,2,-1
      u(i)=v(i)-r3(i,j)*u(i+1)+eps
    4 continue
      do 5 i=iz+2,nz+1
      u(i)=v(i)-r1(i,j)*u(i-1)+eps
    5 continue
    6 continue
c
      return
      end subroutine solve

c     ----------------------------------------------------------------
c
c     Matrix updates.
c
      subroutine updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,rmax,c0,
     >   k0,r,rp,rs,rb,zb,cw,rhow,attnw,cb,rhob,attn,alpw,alpb,ksq,
     >   ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
     > bind(c)
      integer(ram_int) :: mr,mz,nz,mp,np,iz,ib
      complex(ram_complex) ::
     >   ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),
     >   s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp),ksqw(mz)
      real(ram_real) :: 
     >   dr,dz,eta,omega,rmax,c0,k0,r,rp,rs,
     >   rb(mr),zb(mr),cw(mz),rhow(mz),attnw(mz),cb(mz),rhob(mz),
     >   attn(mz),f1(mz),f2(mz),f3(mz),alpw(mz),alpb(mz)
c
c     Varying bathymetry.
c
      if(r.ge.rb(ib+1))ib=ib+1
      jz=iz
      z=zb(ib)+(r+0.5*dr-rb(ib))*(zb(ib+1)-zb(ib))/(rb(ib+1)-rb(ib))
      iz=1.0+z/dz
      iz=max(2,iz)
      iz=min(nz,iz)
      if(iz.ne.jz)call matrc(mz,nz,mp,np,iz,jz,dz,k0,rhow,rhob,alpw,
     >    alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c     Varying profiles.
c
      if(r.ge.rp)then
        call profl(mz,nz,dz,eta,omega,rmax,c0,k0,rp,cw,rhow,attnw,
     >      cb,rhob,attn,alpw,alpb,ksqw,ksqb)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      end if
c
c     Turn off the stability constraints.
c
      if(r.ge.rs)then
        ns=0
        rs=2.0*rmax
        call epade(mp,np,ns,1,k0,c0,dr,pd1,pd2)
        call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      end if
c
      return
      end subroutine updat

c     ----------------------------------------------------------------
c
c     The self-starter.
c
      subroutine selfs(mz,nz,mp,np,ns,iz,zs,dr,dz,pi,c0,k0,rhow,rhob,
     >   alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
     >   pd1,pd2)
     > bind(c)
      integer(ram_int) :: mz,nz,mp,np,ns,iz
      complex(ram_complex) ::
     >   u(mz),v(mz),ksq(mz),ksqb(mz),r1(mz,mp),r2(mz,mp),ksqw(mz),
     >   r3(mz,mp),s1(mz,mp),s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
      real(ram_real) :: zs,dr,dz,pi,c0,
     >   k0,rhob(mz),alpw(mz),alpb(mz),f1(mz),f2(mz),f3(mz),
     >   rhow(mz)
c
c     Conditions for the delta function.
c
      si=1.0+zs/dz
      is=int(si)
      dis=si-is
      u(is)=(1.0-dis)*sqrt(2.0*pi/k0)/(dz*alpw(is))
      u(is+1)=dis*sqrt(2.0*pi/k0)/(dz*alpw(is))
c
c     Divide the delta function by (1-X)**2 to get a smooth rhs.
c
      pd1(1)=0.0
      pd2(1)=-1.0
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
      call solve(mz,nz,mp,1,iz,u,v,r1,r2,r3,s1,s2,s3)
c
c    Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
c
      call epade(mp,np,ns,2,k0,c0,dr,pd1,pd2)
      call matrc(mz,nz,mp,np,iz,iz,dz,k0,rhow,rhob,alpw,
     >   alpb,ksq,ksqw,ksqb,f1,f2,f3,r1,r2,r3,s1,s2,s3,pd1,pd2)
      call solve(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
c
      return
      end subroutine selfs

c     ----------------------------------------------------------------
c
c     Gaussian tapered point source
c
      subroutine pt_src_gtaper(mz,nz,dz,freq,c0,zs,x0,x1,
     >                         theta0,dtheta,p )
     > bind(c)
c
c A beamed point source in 3-dimensions, with a Gaussian beam pattern,
c is given by
c
c psi3D(z) = exp( -(theta-theta0)^2/(2*dtheta^2) ) * exp( ikr ) / r
c 
c where theta = atan( (z-zs)/ (x1-x0) )
c   and     r = sqrt( (x1-x0)^2 + (z-zs)^2 ) . k = 2pi*freq/c0 .
c
c For a PE source, this must be converted to one that is 2D, and
c cylindrically symmetric. This conversion involves the
c asymptotic form of the H_0 Hankel funtion. We output
c
c   p(z)   = psi3D(z) * sqrt(x-x0) *  exp( -ik(x-x0) ) .
c
c       where psi3D = psi3D_1 - psi3D_2
c      (difference of source and image_source)
c
c p(z) will be evaluated along the z_grid at
c     a range x1 for a point source at x0,zs.
c
c   4Mar15, modified for inclusion in RAM
c --
      integer(ram_int) :: mz, nz
      real(ram_real) ::
     >    dz, freq, c0, zs, x0, x1, theta0d, dthetad

      complex(ram_complex) ::  p(mz)

c     Locals
      complex(ram_complex) :: psi3d_1
      complex(ram_complex) :: psi3d_2
      complex(ram_complex) :: factr
      real(ram_real) ::  pi, k0, dx, dx_sq, two_dtheta_sq, 
     >    z, theta_1, theta_2, r_1, r_2, theta0, dtheta
      integer (ram_int) :: nd

      pi       = 4.d0 * atan(1.0d0)
      k0       = 2.d0 * pi * freq/c0
      dx       = x1 - x0
      dx_sq    = dx * dx  
      dtheta = pi * dthetad/180.d0		! convert to radians
      theta0 = pi * theta0d/180.d0
      two_dtheta_sq = 2.0 * dtheta*dtheta

      factr    =  sqrt(dx) * exp( cmplx(0.,-k0*dx) )

      do nd=1, nz

         z = (nd-1)*dz

         r_1 = sqrt( dx_sq + (z-zs)**2 )
         r_2 = sqrt( dx_sq + (z+zs)**2 )

         theta_1 = atan( (z-zs) / dx )
         theta_2 = atan( (z+zs) / dx )

         psi3d_1 = exp( -(theta_1-theta0)**2 / (two_dtheta_sq) )
     +                    * exp( cmplx(0.,k0*r_1) )/ r_1 

         psi3d_2 = exp( -(theta_2+theta0)**2 / (two_dtheta_sq) )
     +                    * exp( cmplx(0.,k0*r_2) )/ r_2 
	  

         p(nd) = (psi3d_1 - psi3d_2) * factr

      end do

      return
      end subroutine pt_src_gtaper

c     ----------------------------------------------------------------
c
c     Output transmission loss.
c
      subroutine outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,
     >                 f3,u,tlg,preal,pimag)
     > bind(c)
      integer(ram_int) :: mz,mdr,ndr,ndz,iz,nzplt,lz,ir
      complex(ram_complex) :: u(mz)
      real(ram_real) :: f3(mz),dir,eps,r
      real(c_float) :: tlg(mz),preal(mz),pimag(mz)
      
c     Locals
      complex(ram_complex) :: ur
c
      ur=(1.0-dir)*f3(ir)*u(ir)+dir*f3(ir+1)*u(ir+1)
      tl=-20.0*log10(abs(ur)+eps)+10.0*log10(r+eps)
      write(2,*)r,tl
c
      if(mdr.eq.ndr)then
        mdr=0
c
        j=0
        do 1 i=ndz,nzplt,ndz
            ur=u(i)*f3(i)
            j=j+1
            tlg(j)=-20.0*log10(abs(ur)+eps)+10.0*log10(r+eps)
            preal(j) = real((ur+eps)/sqrt(r+eps))
c                                               *exp(cmplx(0.,k0*r)))
            pimag(j) = aimag((ur+eps)/sqrt(r+eps))
c                                               *exp(cmplx(0.,k0*r)))
    1   continue
c        write(3)(tlg(j),j=1,lz)
        write(3)(preal(j),j=1,lz)
        write(3)(pimag(j),j=1,lz)

      end if
      mdr=mdr+1
c
      return
      end subroutine outpt

c     ----------------------------------------------------------------
c
c     The coefficients of the rational approximation.
c
      subroutine epade(mp,np,ns,ip,k0,c0,dr,pd1,pd2)
     > bind(c)
      implicit none
      integer(ram_int) :: mp,np,ns,ip
      real(ram_real) :: k0,c0,dr
      complex(ram_complex) :: pd1(mp),pd2(mp)
c
      real*8 pi,sig,alp
      complex*16 z1,z2,g
      real*8 nu
      integer m
      parameter (m=40)
      real*8 bin(m,m),fact(m)
      complex*16 a(m,m),b(m),dg(m),dh1(m),dh2(m),dh3(m)
      integer i,j,n
      pi=4.0d0*datan(1.0d0)
      sig=k0*dr
      n=2*np
c
      if(ip.eq.1)then
c       Called from setup or updat
        nu=0.0d0
        alp=0.0d0
      else
c       Called from selfs (self-starter)
        nu=1.0d0
        alp=-0.25d0
      end if
c
c     The factorials.
c
      fact(1)=1.0d0
      do 1 i=2,n
        fact(i)=dfloat(i)*fact(i-1)
    1 continue
c
c     The binomial coefficients.
c
      do 2 i=1,n+1
        bin(i,1)=1.0d0
        bin(i,i)=1.0d0
    2 continue
      do 4 i=3,n+1
        do 3 j=2,i-1
            bin(i,j)=bin(i-1,j-1)+bin(i-1,j)
    3   continue
    4 continue
c
      do 6 i=1,n
        do 5 j=1,n
            a(i,j)=0.0d0
    5   continue
    6 continue
c
c     The accuracy constraints.
c
      call deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
c
      do 7 i=1,n
        b(i)=dg(i+1)
    7 continue
      do 9 i=1,n
        if(2*i-1.le.n)a(i,2*i-1)=fact(i)
        do 8 j=1,i
            if(2*j.le.n)a(i,2*j)=-bin(i+1,j+1)*fact(j)*dg(i-j+1)
    8   continue
    9 continue
c
c     The stability constraints.
c
      if(ns.ge.1)then
      z1=-3.0d0
      b(n)=-1.0d0
      do 10 j=1,np
        a(n,2*j-1)=z1**j
        a(n,2*j)=0.0d0
   10 continue
      end if
c
      if(ns.ge.2)then
        z1=-1.5d0
        b(n-1)=-1.0d0
        do 11 j=1,np
            a(n-1,2*j-1)=z1**j
            a(n-1,2*j)=0.0d0
   11   continue
      end if
c
      call gauss(m,n,a,b)
c
      dh1(1)=1.0d0
      do 12 j=1,np
        dh1(j+1)=b(2*j-1)
   12 continue
      call fndrt(dh1,np,dh2,m)
      do 13 j=1,np
        pd1(j)=-1.0d0/dh2(j)
   13 continue
c
      dh1(1)=1.0d0
      do 14 j=1,np
        dh1(j+1)=b(2*j)
   14 continue
      call fndrt(dh1,np,dh2,m)
      do 15 j=1,np
        pd2(j)=-1.0d0/dh2(j)
   15 continue
c
      return
      end subroutine epade

c     ----------------------------------------------------------------
c
c     The operator function.
c
      function g(sig,x,alp,nu)
     > bind(c)
      complex(c_double_complex) :: g,arg
      real(c_double) :: alp,sig,x,nu
      
      arg = cmplx( alp*log(1.0+x), sig*(-1.0d0+sqrt(1.0+x)),
     >  ram_complex )
      g=(1.0d0-nu*x)**2*exp( arg )
      return
      end function g

c     ----------------------------------------------------------------
c
c     The derivatives of the operator function at x=0.
c
      subroutine deriv(m,n,sig,alp,dg,dh1,dh2,dh3,bin,nu)
     > bind(c)
      implicit none
      integer(ram_int) :: m,n
      real(c_double) :: sig,alp,bin(m,m),nu
      complex(c_double_complex) :: dg(m),dh1(m),dh2(m),dh3(m)

      integer i,j
      real*8 exp1,exp2,exp3
      complex*16 ci
      ci=dcmplx(0.0d0,1.0d0)
c
      dh1(1)=0.5d0*ci*sig
      exp1=-0.5d0
      dh2(1)=alp
      exp2=-1.0d0
      dh3(1)=-2.0d0*nu
      exp3=-1.0d0
      do 1 i=2,n
        dh1(i)=dh1(i-1)*exp1
        exp1=exp1-1.0d0
        dh2(i)=dh2(i-1)*exp2
        exp2=exp2-1.0d0
        dh3(i)=-nu*dh3(i-1)*exp3
        exp3=exp3-1.0d0
    1 continue
c
      dg(1)=1.0d0
      dg(2)=dh1(1)+dh2(1)+dh3(1)
      do 3 i=2,n
        dg(i+1)=dh1(i)+dh2(i)+dh3(i)
        do 2 j=1,i-1
            dg(i+1)=dg(i+1)+bin(i,j)*(dh1(j)+dh2(j)+dh3(j))*dg(i-j+1)
    2   continue
    3 continue
c
      return
      end subroutine deriv

c     ----------------------------------------------------------------
c
c     Gaussian elimination.
c
      subroutine gauss(m,n,a,b)
     > bind(c)
      implicit none
      integer(ram_int) :: m,n
      complex(c_double_complex) :: a(m,m),b(m)
c
c     Downward elimination.
c
      integer i,j,k
      do 4 i=1,n
        if(i.lt.n)call pivot(m,n,i,a,b)
        a(i,i)=1.0d0/a(i,i)
        b(i)=b(i)*a(i,i)
        if(i.lt.n)then
            do 1 j=i+1,n
                a(i,j)=a(i,j)*a(i,i)
    1       continue
            do 3 k=i+1,n
                b(k)=b(k)-a(k,i)*b(i)
                do 2 j=i+1,n
                    a(k,j)=a(k,j)-a(k,i)*a(i,j)
    2           continue
    3       continue
        end if
    4 continue
c
c     Back substitution.
c
      do 6 i=n-1,1,-1
        do 5 j=i+1,n
            b(i)=b(i)-a(i,j)*b(j)
    5   continue
    6 continue
c
      return
      end subroutine gauss

c     ----------------------------------------------------------------
c
c     Rows are interchanged for stability.
c
      subroutine pivot(m,n,i,a,b)
     > bind(c)
      implicit none
      integer(ram_int) :: m,n,i
      complex(c_double_complex) :: a(m,m),b(m)
      
      integer i0,j
      real*8 amp0,amp
      complex*16 temp
c
      i0=i
      amp0=cdabs(a(i,i))
      do 1 j=i+1,n
        amp=cdabs(a(j,i))
        if(amp.gt.amp0)then
            i0=j
            amp0=amp
        end if
    1 continue
      if(i0.eq.i)return
c
      temp=b(i)
      b(i)=b(i0)
      b(i0)=temp
      do 2 j=i,n
        temp=a(i,j)
        a(i,j)=a(i0,j)
        a(i0,j)=temp
    2 continue
c
      return
      end subroutine pivot

c     ----------------------------------------------------------------
c
c     The root-finding subroutine. 
c
      subroutine fndrt(a,n,z,m)
     > bind(c)
      integer(ram_int) :: n,m
      complex(c_double_complex) :: a(m),z(m)

      complex*16 root
      real*8 err
c
      if(n.eq.1)then
        z(1)=-a(1)/a(2)
        return
      end if
      if(n.eq.2)go to 4
c
        do 3 k=n,3,-1
c
c           Obtain an approximate root.
c
            root=0.0d0
            err=1.0d-12
            call guerre(a,k,m,root,err,1000)
c
c           Refine the root by iterating five more times.
c
            err=0.0d0
            call guerre(a,k,m,root,err,5)
            z(k)=root
c
c           Divide out the factor (z-root).
c
            do 1 i=k,1,-1
                a(i)=a(i)+root*a(i+1)
    1       continue
            do 2 i=1,k
                a(i)=a(i+1)
    2       continue
c
    3   continue
c
c     Solve the quadratic equation.
c
    4 continue
      z(2)=0.5*(-a(2)+sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
      z(1)=0.5*(-a(2)-sqrt(a(2)**2-4.0*a(1)*a(3)))/a(3)
c
      return
      end subroutine fndrt

c     ----------------------------------------------------------------
c
c     This subroutine finds a root of a polynomial of degree n > 2
c     by Laguerre's method.
c
      subroutine guerre(a,n,m,z,err,nter)
     > bind(c)
      integer(ram_int) :: n,m,nter
      real(c_double) :: err
      complex(c_double_complex) :: a(m),z

      complex*16 az(50),azz(50),dz,p,pz,pzz,f,g,h,ci
      real*8 amp1,amp2,rn,eps
      ci=cmplx(0.0d0,1.0d0)
      eps=1.0d-20
      rn=real(n)
c
c     The coefficients of p'(z) and p''(z).
c
      do 1 i=1,n
        az(i)=i*a(i+1)
    1 continue
      do 2 i=1,n-1
        azz(i)=i*az(i+1)
    2 continue
c
      iter=0
    3 continue
        p=a(n)+a(n+1)*z
        do 4 i=n-1,1,-1
            p=a(i)+z*p
    4   continue
        if(abs(p).lt.eps)return
c
        pz=az(n-1)+az(n)*z
        do 5 i=n-2,1,-1
            pz=az(i)+z*pz
    5   continue
c
        pzz=azz(n-2)+azz(n-1)*z
        do 6 i=n-3,1,-1
            pzz=azz(i)+z*pzz
    6   continue
c
c       The Laguerre perturbation.
c
        f=pz/p
        g=f**2-pzz/p
        h=sqrt((rn-1.0d0)*(rn*g-f**2))
        amp1=abs(f+h)
        amp2=abs(f-h)
        if(amp1.gt.amp2)then
            dz=-rn/(f+h)
        else
            dz=-rn/(f-h)
        end if
c
        iter=iter+1
c
c       Rotate by 90 degrees to avoid limit cycles. 
c
        jter=jter+1
        if(jter.eq.10)then
            jter=1
            dz=dz*ci
        end if
        z=z+dz
c
        if(iter.eq.200)then
            write(*,*)' '
            write(*,*)'   Laguerre method not converging.'
            write(*,*)'   Try a different combination of DR and NP.'
            write(*,*)' '
            stop
        end if
c
        if((abs(dz).gt.err).and.(iter.lt.nter))go to 3
c
      return
      end subroutine guerre
      
c     ----------------------------------------------------------------
c
c     Close the input and output files used by RAMair (units 1, 2, and 3)
c
      subroutine close_ram_files()
     > bind(c)
      
      close( unit=1 )
      close( unit=2 )
      close( unit=3 )
      return
      
      end subroutine close_ram_files
      
      end module ramlib
