
!------------------------------------------------------------------------------!
! Module for SZ cluster counts likelihood
!
! Module prepared for cosmomcplanck by A. Bonaldi 2015
! This module have been used to produce the results in
! Planck 2015 results XXIV: Cosmology from Sunyaev-Zeldovich cluster counts
! for questions and problems contact anna.bonaldi@manchester.ac.uk
!
! Original code to compute cluster counts written from 2001
!
!                   v 0.1
!       Jochen Weller and Richard Battye
!                  4/2/2007
!    all the necessary modules and subroutines
!
!                   v 1.0
!               Jochen Weller
!                 10/1/2008
!       modified for optical cluster counts
!
!                   v 2.0
!               Jochen Weller
!                 30/12/2009
!         modified for SZ cluster counts
!
!                   v 3.0
!             Anna Bonaldi 2011-2013
!    Include realistic selection function for Planck
!       missing redshifts, errors on redshifts
!             Matthieu Roman 2012-2013
!               QA completeness
!
!                   v 4.0
!             Anna Bonaldi 2014-2015
!           updates for 2015 analysis
!            dN and dzdq likelihood
!------------------------------------------------------------------------------!

module constants_sz
   use PRECISION
   implicit none
   real(dl), PARAMETER :: Mpc = 3.08568025d22 ! Mpc in metres
   real(dl), PARAMETER :: G = 6.67300d-11 ! Newton's constant in m^3kg^-1s^-2
   real(dl), PARAMETER :: c = 3.0d8 ! speed of light in m/s
   real(dl), PARAMETER :: msun = 1.98892d30 ! mass of sun in kg
   real(dl), PARAMETER :: rhocrit0 = 2.7751973751261264d11 ! rhocrit in units of h^-1 Msun/ h^-3 Mpc^3
   real(dl), PARAMETER :: PI=3.141592653589793238463D0
end module constants_sz

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module cosmology
   ! module computing cosmological functions
   USE precision
   USE CONSTANTS_sz
   use Interpolation
   implicit none
   ! note that we switch to the generalized Linder parametrization now
   public
   TYPE cospar
   REAL(dl) :: H0,w0,w1,omegam,omegav,n,sig8,omegak,omegabh2,gamma, s8omegamp3
   real(dl) :: ystar,alpha,sigmaM,bias,biasinv,logystar,beta
   real(dl) :: bias_act, a_act, b_act, d_act, ylim_act  ! :D ! for ACT
   real(dl) :: bias_spt, a_spt, b_spt, c_spt, d_spt     ! :D ! for SPT
   class(TCubicSpline), pointer :: sigmaR
   END TYPE cospar
   Type (cospar), SAVE :: cosmopar

   contains

   !---------------------------------------------------------------------------!

   function Eh(z)
      ! E=H(z)/H0
      real(dl) :: Eh
      real(dl), intent(in) :: z
      real(dl) ::  a
      a = 1.0/(1.0+z)
      Eh=sqrt(cosmopar%omegam*a**(-3)+cosmopar%omegav*a**(-3.0*(1.0+cosmopar%w0+cosmopar%w1))*exp(3.0*cosmopar%w1*(a-1.0))+cosmopar%omegak*a**(-2))
      return
   end function Eh

   !---------------------------------------------------------------------------!

   function Omegam(z)
      real(dl), intent(in) :: z
      real(dl) :: Omegam
      Omegam=cosmopar%omegam*(1.0+z)**3/(Eh(z))**2
      return
   end function Omegam

   !---------------------------------------------------------------------------!

   function Omegade(z)
      real(dl), intent(in) :: z
      real(dl) :: Omegade
      Omegade=cosmopar%omegav*(1.0+z)**(3.0*(1.0+cosmopar%w0+cosmopar%w1))*exp(-3.0*cosmopar%w1*z/(1.0+z))/(Eh(z))**2
      return
   end function Omegade

   !---------------------------------------------------------------------------!

   function Omegak(z)
      real(dl), intent(in) :: z
      real(dl) :: Omegak
      Omegak=1.0_dl-Omegam(z)-Omegade(z)
      return
   end function Omegak

   !---------------------------------------------------------------------------!

   function r(z)
      ! coordinate distance , in units of h^-1 Mpc
      real(dl), INTENT(IN) :: z
      real(dl) :: r
      real(dl), PARAMETER :: tol=1.0d-5
      real(dl) :: integral,rombint
      external rombint
      integral=rombint(rint,0._dl,z,tol)
      ! maybe soften this to 10^-5
      if (cosmopar%omegak == 0._dl) then
         r= c*integral/1.0d5
         elseif (cosmopar%omegak > 0._dl) then
            r=c/sqrt(cosmopar%omegak)*dsinh(sqrt(cosmopar%omegak)*integral)/1.0d5
         else
            r=c/sqrt(-cosmopar%omegak)*dsin(sqrt(-cosmopar%omegak)*integral)/1.0d5
      end if
      return
   end function r

   !---------------------------------------------------------------------------!

   function rint(z)
      ! integrand for coordinate distance
      real(dl), INTENT(IN) :: z
      real(dl) :: rint
      rint = 1._dl/Eh(z)
      return
   end function rint

   !---------------------------------------------------------------------------!

   function dA(z)
      ! angular diameter distance in units of h^-1 Mpc
      real(dl), INTENT(IN) :: z
      real(dl) :: dA
      dA = r(z)/(1._dl+z)
      return
   end function dA

   !---------------------------------------------------------------------------!

   function dVdzdO(z)
      ! volume element in units of h^-3 Mpc^3
      real(dl), INTENT(IN) :: z
      real(dl) :: dVdzdO
      dVdzdO = c/1.0d5*r(z)**2/Eh(z)
      return
   end function dVdzdO

   !---------------------------------------------------------------------------!

end module cosmology

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module survey
   use PRECISION
   implicit none
   TYPE survpar
   ! Note that in the SZ version mlimin is only used for the number of mass bins
   REAL(dl) :: ylimin,deg2,ymaxin, deg2_a, deg2_s, deg2_SOfull, deg2_SOtest
   INTEGER :: ybin,Nscat
   REAL(dl) :: sfid,ab,nb
   REAL(dl), allocatable :: b(:)
   ! REAL(dl), allocatable :: mbias(:),mscatter(:)
   END TYPE survpar
   TYPE (survpar), SAVE :: surveypar

end module survey

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module massobservable
   !scaling relations
   USE PRECISION
   USE SURVEY
   USE COSMOLOGY
   implicit none
   real(dl), parameter :: thetastar = 6.997 ! dic 2012 new scalings
   real(dl), parameter :: alpha_theta = 1./3.
   REAL(dl), parameter :: q = 6. ! signal-to-noise threshold for planck
   real(dl), allocatable :: Mlim_arr_a(:,:), Mlim_arr_s(:), y0_act_arr(:,:) ! :D ! pre-calculation of Mlim of ACT and SPT
   real(dl), allocatable :: y0_SO_arr(:,:), surveydeg2_SO(:), surveydeg2_act(:)
   real(dl), allocatable :: y0_SO_tiles_arr(:,:,:), surveydeg2_SO_tiles(:), precalc_y0(:,:,:)

   contains

   !---------------------------------------------------------------------------!

   ! Planck scaling relation parameter
   function theta500(m,z)
      !size-mass/redshift relation
      real(dl), INTENT(IN) :: z,m
      real(dl) :: theta500,thetastar2,m2
      m2=m*cosmopar%bias !hydro bias
      thetastar2=thetastar*(cosmopar%H0/70.)**(-2./3.)
      theta500=thetastar2*(m2/3.e14*(100/cosmopar%H0))**alpha_theta*Eh(z)**(-2./3)*(100.0*da(z)/500.0/cosmopar%H0)**(-1.)
      RETURN
   end function theta500

   !---------------------------------------------------------------------------!

   ! Planck scaling relation parameter
   function y500(m,z)
      !y-mass/redshift relation
      real(dl), INTENT(IN) :: z,m
      real(dl) :: y500,ystar2,m2,alpha
      m2=m*cosmopar%bias !hydro bias
      alpha=cosmopar%alpha
      ystar2=cosmopar%ystar
      ystar2=ystar2*(cosmopar%H0/70.)**(-2.+alpha)
      y500=ystar2*(m2/3.0e14*(100./cosmopar%H0))**alpha*Eh(z)**(cosmopar%beta)*(100.0*da(z)/500.0/cosmopar%H0)**(-2.)
      RETURN
   end function y500

   !---------------------------------------------------------------------------!

   ! :D ! ACT scaling paramter : theta_500(m,z)
   ! analytical cluster size theta as a function of mass and redshift
   !
   ! (slightly different definitions between here and test code in python)
   ! H0  (.py) = cosmopar%H0*1.e3/Mpc
   ! D_a (.py) = dA(z)/cosmopar%H0*1.e2*Mpc
   ! 1 rad     = 3437.75 arcmin
   function theta_act(m, z)
      real(dl), intent(in) :: z
      real(dl) :: m, m2, theta_act, hh

      !introducing mass bias
      m2 = m*cosmopar%bias_act
      !m2 = 1.54*(m2**0.72)*((1.e14*msun)**(1. - 0.72)) !fitting
      theta_act = (G*m2/250._dl)**(1./3)*(Eh(z)*cosmopar%H0*1.e3/Mpc)**(-2./3)*(100._dl*dA(z)*Mpc/cosmopar%H0)**(-1.)*3437.75_dl

      return
   end function theta_act

   !---------------------------------------------------------------------------!

   ! :D ! ACT scaling parameter : mismatch function Q(theta(m,z))
   ! mismatch between the true cluster size and the size encoded in the filter
   ! response function to reconstruct the cluster central decrement
   ! function of cluster angular size theta_act(m, z)
   ! when Q = 1, filter is perfectly matched
   ! sadly, interpolated Q goes quite crazy when both ends of theta values

   function Qfunc(m, z)
      real(dl), intent(in) :: z
      real(dl) :: m, Qfunc, col1, col2, splQ, tt
      real(dl), allocatable :: theta(:), Q0(:), Q1(:)
      integer :: iostat, i, nQrow, io

      ! reading Qfit file

      open(unit=33, file='data/Qfit.txt')                             ! ACTPol map (99 -> 103)
      !open(unit=33, file='data/MFMF_SOSim_3freq_small_Qfit.txt')     ! SO testmap with a single tile (51)
      !open(unit=33, file='data/SOSim_3freq_tiles_Qfit_mean.txt')     ! SO fullmap with one average Q func (51)

      nQrow = 0
      do
        read(33, *, iostat=io)
        if (io/=0) exit
        nQrow = nQrow + 1
      enddo
      nQrow = nQrow - 1

      if (allocated(theta)) deallocate(theta)
      if (allocated(Q0)) deallocate(Q0)
      if (allocated(Q1)) deallocate(Q1)
      allocate(theta(nQrow), Q0(nQrow), Q1(nQrow), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      endif

      rewind(33)
      do i = 1, nQrow
         read(33,*) col1, col2
         theta(i) = col1
         Q0(i) = col2
      end do
      close(unit=33)

      ! subroutine spline in Interpolation.f90
      ! calculates array of the second derivatives used by cubic spline interpolation
      ! input: xdata, ydata, number, 1st derivatives at end points and output
      ! need to check about end points of 1st derivative
      ! it is very not okay there

      ! need Q1 in order to obtain splQ
      call spline(theta, Q0, nQrow, 0.3_dl, -0.017_dl, Q1)
      !call spline(theta, Q0, nQrow, 10.0_dl, 0.1_dl, Q1)             ! for SO Qfunc

      tt = theta_act(m, z)
      !print*, z, log(m/msun), tt
      ! min tt = 0.86533721
      ! max tt = 25.3004629

      ! check if theta is valid - do I not need this anymore?
      !if (tt > 30.) tt = 30.
      !if (tt < 0.5) tt = 0.5

      ! function splintp gives new Q for a given theta_act
      Qfunc = splintp(theta, Q0, Q1, nQrow, tt, splQ)
      !print*, tt, Qfunc
      ! min Q = 0.3065542339
      ! max Q = 1.3497906666
      !if (Qfunc < 0) Qfunc = 1.e-3

      return
   end function Qfunc

   !---------------------------------------------------------------------------!

   function Qfunc_SO(m, z, theta, Q0)
      real(dl), intent(in) :: z, theta(:), Q0(:)
      !real(dl) :: m, Qfunc_SO, col1, col2, splQ, tt
      real(dl) :: m, Qfunc_SO, splQ, tt
      real(dl), allocatable :: Q1(:)
      integer :: iostat
      integer, parameter :: nQrow = 51

      !open(unit=44, file='data/MFMF_SOSim_3freq_small_Qfit.txt')     ! SO testmap with a single tile (51)
      !open(unit=44, file='data/SOSim_3freq_tiles_Qfit_mean.txt')    ! SO fullmap with one average Q func (51)

      !nQrow = 0
      !do
      !  read(44, *, iostat=io)
      !  if (io/=0) exit
      !  nQrow = nQrow + 1
      !enddo
      !print*, nQrow

      !if (allocated(theta)) deallocate(theta)
      !if (allocated(Q0)) deallocate(Q0)
      if (allocated(Q1)) deallocate(Q1)
      allocate(Q1(nQrow), stat=iostat)
      if (iostat/=0) then
         print*,'allocation error'
      endif

      !rewind(44)
      !do i = 1, nQrow
      !   read(44,*) col1, col2
      !   theta(i) = col1
      !   Q0(i) = col2
      !end do
      !close(unit=44)

      ! subroutine spline in Interpolation.f90
      ! calculates array of the second derivatives used by cubic spline interpolation
      ! input: xdata, ydata, number, 1st derivatives at end points and output

      ! need Q1 in order to obtain splQ
      call spline(theta, Q0, nQrow, 10.0_dl, 0.1_dl, Q1)

      tt = theta_act(m, z)

      ! function splintp gives new Q for a given theta_act
      Qfunc_SO = splintp(theta, Q0, Q1, nQrow, tt, splQ)

      deallocate(Q1)

      return
   end function Qfunc_SO

   !---------------------------------------------------------------------------!

   function Qfunc_SO_tiles(m, z, tile, theta, Q00)
      real(dl), intent(in) :: z, theta(:), Q00(:,:)
      real(dl) :: m, Qfunc_SO_tiles, splQ, tt
      real(dl), allocatable :: Q1(:)
      integer, intent(in) :: tile
      integer :: iostat
      integer, parameter :: nQrow = 51
      integer, parameter :: ncol = 10  ! 2  346 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (allocated(Q1)) deallocate(Q1)
      allocate(Q1(nQrow), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      endif

      ! need Q1 in order to obtain splQ
      call spline(theta, Q00(:,tile), nQrow, 10.0_dl, 0.1_dl, Q1)

      tt = theta_act(m, z)

      ! function splintp gives new Q for a given theta_act
      Qfunc_SO_tiles = splintp(theta, Q00(:,tile), Q1, nQrow, tt, splQ)

      deallocate(Q1)

      return
   end function Qfunc_SO_tiles

   !---------------------------------------------------------------------------!

   ! :D ! ACT scaling parameter : relativistic correction f_rel(m,z)
   ! from Hasselfield (2013) page 5

   function relfn(m, z)
      real(dl), intent(in) :: z
      real(dl) :: m, m2, t, relfn, hh

      m2    = m*cosmopar%bias_act

      t     = -0.00848_dl*(m2/(3.e14*msun*70./cosmopar%H0)*Eh(z))**(-0.585_dl)
      !t     = -0.00848_dl*(((1.54_dl/(3.*70./cosmopar%H0))*(m2/(1.e14*msun))**0.72_dl)*Eh(z))**(-0.585_dl) ! fitting
      relfn = 1.0_dl + (3.79_dl*t) - 28.2_dl*(t**2.)

      return
   end function relfn

   !---------------------------------------------------------------------------!

   ! :D ! to calculate a limiting mass from ACT scaling relation
   ! need to find a root of this function
   ! a root will be a mass of a given redshift z and limiting y
   ! ylim is fixed here (constant y0lim)
   ! be careful with the unit of y0lim_const

   function func(m, z)
      real(dl), intent(in) :: z
      real(dl) :: m, m2, A_a, B_a, func, y0lim_const

      A_a  = cosmopar%a_act    ! 4.95e-5
      B_a  = cosmopar%b_act    ! 0.08        from Hilton (2017)

      !when using Planck scaling relation
      !A_a  = (10.**cosmopar%logystar)/((4.38e3)*(2.**cosmopar%alpha))
      !B_a  = cosmopar%alpha - 5./3.

      m2     = m*cosmopar%bias_act
      y0lim_const = cosmopar%ylim_act ! from ini file (constant ylim) : ylim_SZact

      func = A_a*Eh(z)**2.*(m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a)*Qfunc(m, z)*relfn(m, z) - y0lim_const*1.e-4

      return
   end function func

   !---------------------------------------------------------------------------!

   ! :D ! same as above, but y0lim is from the map now
   function func_ymap(m, z, y0lim)
      real(dl), intent(in) :: z, y0lim
      real(dl) :: m, m2, A_a, B_a, func_ymap

      A_a  = cosmopar%a_act
      B_a  = cosmopar%b_act

      !when using Planck scaling relation
      !A_a  = (10.**cosmopar%logystar)/((4.38e3)*(2.**cosmopar%alpha))
      !B_a  = cosmopar%alpha - 5./3.

      m2   = m*cosmopar%bias_act

      func_ymap = A_a*(Eh(z)**2.)*((m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a))*Qfunc(m, z)*relfn(m, z) - y0lim*1.e-4

      !print*, m, relfn(m,z)

      return
   end function func_ymap

   !---------------------------------------------------------------------------!

   ! below is just copied and pasted
   ! for the constant ylim case
   ! NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING p1186
   ! find the root of a function func thought to lie between x1 and x2.
   ! The root, returned as secant, is refined until its accuracy is ±xacc.

   function secant(z, x1, x2, xacc)
      implicit none
      real(dl), intent(in) :: z, x1, x2, xacc
      real(dl) :: secant
      integer, parameter :: maxit = 1000 ! maximum number of allowed iterations
      integer :: j
      real(dl) :: dx, f, fl, xl, swap
      fl = func(x1, z)
      f  = func(x2, z)
      if (abs(fl) < abs(f)) then  ! pick the bound with the smaller function value
         secant = x1              ! as the most recent guess
         xl = x2
         swap = fl
         fl = f
         f = swap
      else
         xl = x1
         secant = x2
      end if
      do j = 1, maxit                 ! secant loop
         dx = (xl - secant)*f/(f-fl)  ! increment with respect to latest value
         xl = secant
         fl = f
         secant = secant + dx
         f = func(secant, z)
         if (abs(dx) < xacc .or. f == 0.0) return ! convergence
      end do
   end function secant

   !---------------------------------------------------------------------------!

   ! same as above, but for the case of including limiting y0 map values
   ! additional y0lim parameter, just passing

   function secant_y(z, x1, x2, xacc, y0lim)
      implicit none
      real(dl), intent(in) :: z, x1, x2, xacc, y0lim
      real(dl) :: secant_y
      integer, parameter :: maxit = 1000 ! maximum number of allowed iterations
      integer :: j
      real(dl) :: dx, f, fl, xl, swap
      fl = func_ymap(x1, z, y0lim)
      f  = func_ymap(x2, z, y0lim)
      if (abs(fl) < abs(f)) then  ! pick the bound with the smaller function value
         secant_y = x1            ! as the most recent guess
         xl = x2
         swap = fl
         fl = f
         f = swap
      else
         xl = x1
         secant_y = x2
      end if
      do j = 1, maxit                   ! secant loop
         dx = (xl - secant_y)*f/(f-fl)  ! increment with respect to latest value
         xl = secant_y
         fl = f
         secant_y = secant_y + dx
         f = func_ymap(secant_y, z, y0lim)
         if (abs(dx) < xacc .or. f == 0.0) return ! convergence
      end do
   end function secant_y

   !---------------------------------------------------------------------------!

   ! :D ! finding ACT limiting mass with ylim0_const
   ! ACT limiting mass is a function of redshift

   function Mlim_act(z)
      real(dl), intent(in) :: z
      real(dl) :: Mlim_act, mm1, mm2, macc

      mm1   = 2.e14*msun      ! initial guess
      mm2   = 4.e14*msun      ! based on M_pivot = 3.e14*msun
      macc  = 1.e11*msun      ! accuracy of root m

      Mlim_act = secant(z, mm1, mm2, macc)

      return
   end function Mlim_act

   !---------------------------------------------------------------------------!

   ! :D ! finding ACT limiting mass with limiting y0 map values
   ! ACT limiting mass is a function of redshift and limiting y0

   function Mlim_act_ymap(z, y0lim)
      real(dl), intent(in) :: z, y0lim
      real(dl) :: Mlim_act_ymap, mm1, mm2, macc

      mm1   = 2.e14*msun      ! initial guess
      mm2   = 4.e14*msun      ! based on M_pivot = 3.e14*msun
      macc  = 1.e11*msun      ! accuracy of root m

      Mlim_act_ymap = secant_y(z, mm1, mm2, macc, y0lim)

      return
   end function Mlim_act_ymap

   !---------------------------------------------------------------------------!

   ! :D ! ACT limiting mass is pre-calculated and allocated in an array : Mlim_arr_a
   ! redshift is newly assigned just for this subroutine

   subroutine Mlim_act_allocate(y0lim_arr)
      real(dl), allocatable :: zz_a(:)
      real(dl), intent(in) :: y0lim_arr(:)
      integer :: Nzz_a, i, j, iostat
      real(dl) :: zz0_a, zzmax_a, dzz_a
      integer, parameter :: Ny0lim = 48

      zz0_a = 0.0d0
      zzmax_a = 1.5
      dzz_a = 0.1d0
      Nzz_a = int((zzmax_a - zz0_a)/dzz_a) + 2

      if (allocated(zz_a)) deallocate(zz_a)
      if (allocated(Mlim_arr_a)) deallocate(Mlim_arr_a)
      allocate(zz_a(Nzz_a), Mlim_arr_a(Nzz_a, Ny0lim), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      endif

      ! allocate redshift for ACT
      do i = 0, Nzz_a-1
         zz_a(i+1) = zz0_a + i*dzz_a !+ 0.5_dl*dzz_a
      end do
      if (zz0_a == 0._dl) zz_a(1) = zz_a(1) + 1.e-8

      ! allocate ACT limiting mass for given redshift and ylim
      do i = 1, Nzz_a
         do j = 1, Ny0lim
            Mlim_arr_a(i, j) = Mlim_act_ymap(zz_a(i), y0lim_arr(j))
            !print*, i, j, log(Mlim_arr_a(i, j)/msun)
            !if (i==6 .and. j==48) then
            ! print*, Mlim_act_ymap(zz_a(i), y0lim_arr(j))/msun
            !end if
         end do
      end do

      return
   end subroutine Mlim_act_allocate

   !---------------------------------------------------------------------------!

   ! :D ! manually deallocate Mlim_arr_a
   subroutine Mlim_act_deallocate
      if (allocated(Mlim_arr_a)) deallocate(Mlim_arr_a)
   end subroutine Mlim_act_deallocate

   !---------------------------------------------------------------------------!

   ! :D ! find a matched index in y0lim_arr
   ! only called by scatter_act_ymap

   function ylim_index(y0lim, y0lim_arr)
      real(dl), intent(in) :: y0lim, y0lim_arr(:)
      integer, parameter :: Narr = 48
      integer :: i, ylim_index

      do i = 1, Narr
         if (y0lim_arr(i) .eq. y0lim) then
            ylim_index = i
         end if
      end do

      return
   end function ylim_index

   !---------------------------------------------------------------------------!

   ! :D ! scatter around ACT limiting mass in the case of including limiting y0 map
   ! mlim_scatter_act_ymap = [0, 1]
   ! this is just a simple version of completeness

   function mlim_scatter_act_ymap(m, z, y0lim, y0lim_arr)
      real(dl), intent(in) :: m, z, y0lim, y0lim_arr(:)
      real(dl) :: D_a, arg, mlim_scatter_act_ymap
      integer :: index_z, index_y

      D_a = 0.1 ! given by me :D

      index_z = int(10.*z) + 1
      index_y = ylim_index(y0lim, y0lim_arr) ! need to find an index that's pointing the same one

      if (allocated(Mlim_arr_a)) then
         arg = (log(m) - log(Mlim_arr_a(index_z, index_y)/msun))/sqrt(2.)/D_a
         mlim_scatter_act_ymap = (1. + erf(arg))/2.
      else
         print*, 'error in scatter_act'
         call Mpistop()
      end if

      return
   end function mlim_scatter_act_ymap

   !---------------------------------------------------------------------------!

   ! :D ! scatter around ACT limiting mass with constant y0lim
   function mlim_scatter_act(m, z)
      real(dl), intent(in) :: m, z
      real(dl) :: D_a, arg, mlim_scatter_act

      D_a = 0.1 ! given by me :D

      arg = (log(m) - log(Mlim_act(z)/msun))/sqrt(2.)/D_a
      mlim_scatter_act = (1. + erf(arg))/2.

      return
   end function mlim_scatter_act

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   function y0_act(m, z)
      implicit none
      real(dl), intent(in) :: m, z
      real(dl) :: m2, A_a, B_a, y0_act

      A_a  = cosmopar%a_act
      B_a  = cosmopar%b_act
      m2   = m*msun*cosmopar%bias_act

      y0_act = A_a*(Eh(z)**2.)*((m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a))*Qfunc(m*msun, z)*relfn(m*msun, z)
      !y0_act = A_a*(Eh(z)**2.)*(((1.54_dl/(3.*70./cosmopar%H0))*(m*cosmopar%bias_act/1.e14)**0.72_dl)**(1. + B_a))*Qfunc(m*msun, z)*relfn(m*msun, z) !fitting

      if (y0_act < 0) y0_act = -1.*y0_act

      return
   end function y0_act

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   ! to calculate ACT mass estimate
   ! :D ! this is used for calculating mass estimates of ACT
   function scatter_y0(m, z, SZcat, y0_cat, y0err_cat)
      real(dl), intent(in) :: m, z, SZcat(:,:), y0_cat, y0err_cat
      real(dl) :: fac, sigma_y0, mu, scatter_y0
      integer :: i

      sigma_y0 = 0.2_dl
      fac = 1./sqrt(2.*pi*(sigma_y0**2.))
      mu = log(y0_act(m, z))
      scatter_y0 = fac*(exp(-((log(y0_cat) - mu)**2.)/(2.*(sigma_y0**2. + y0err_cat**2.))))
      !scatter_y0 = fac*(exp(-((log(y0_cat) - mu)**2.)/(2.*(sigma_y0**2.))))

      return
   end function scatter_y0

   !---------------------------------------------------------------------------!

   subroutine y0_act_allocate(y0_act_arr)
      real(dl), allocatable :: y0_act_arr(:,:), zz(:), mm(:)
      integer :: Nzz, Nmm, iostat, i, j
      real(dl) :: zz_min, zz_max, dzz, mm_min, mm_max, dmm

      zz_min = 0.0d0
      zz_max = 1.5
      dzz = 0.1d0
      Nzz = int((zz_max - zz_min)/dzz) + 2

      mm_min = log(1.e13)
      mm_max = log(1.e16)
      dmm = 0.05_dl
      Nmm = int((mm_max - mm_min)/dmm)

      if (allocated(y0_act_arr)) deallocate(y0_act_arr)
      if (allocated(zz)) deallocate(zz)
      if (allocated(mm)) deallocate(mm)
      allocate(y0_act_arr(Nmm, Nzz), zz(Nzz), mm(Nmm), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      end if

      do i = 0, Nzz-1
         zz(i+1) = zz_min + i*dzz
      end do
      if (zz_min == 0._dl) zz(1) = zz(1) + 1.e-8

      do j = 0, Nmm-1
         mm(j+1) = mm_min + j*dmm
      end do

      do i = 1, Nzz
         do j = 1, Nmm
            y0_act_arr(j, i) = y0_act(exp(mm(j)), zz(i))
            !print*, i, j, y0_act_arr(j, i)
         end do
         !stop
      end do

      return
   end subroutine y0_act_allocate

   !---------------------------------------------------------------------------!

   subroutine y0_act_deallocate
      if (allocated(y0_act_arr)) deallocate(y0_act_arr)
   end subroutine y0_act_deallocate

   !---------------------------------------------------------------------------!

   !function completeness_act(m, z, y0lim, y0lim_arr)
   function completeness_act(m, z, noise_act, noise_act_arr, surveydeg2_act)
      !real(dl), intent(in) :: m, z, y0lim, y0lim_arr(:)
      real(dl), intent(in) :: m, z, noise_act, noise_act_arr(:), surveydeg2_act(:)
      real(dl) :: sel, sels, intsc, completeness_act, totalArea, erfn, sigma, summ
      real(dl) :: D_a, fac, mu, arg, arg0, arg1
      real(dl) :: lnymin, lnymax, dlny, lny, y0, y1, dy, py, yy
      integer :: Ny, Npatch, i, j, iostat, index_mm, index_zz
      real(dl), allocatable :: self(:)

      D_a = cosmopar%d_act
      !D_a = cosmopar%sigmaM     ! when using Planck YM relation to ACT SZ cluster data

      Npatch = size(noise_act_arr)
      totalArea = sum(surveydeg2_act)
      !print*, totalArea/3.046174198d-4, Npatch
      !stop

      index_zz = int(10.*z) + 1
      index_mm = nint((log(m)- log(1.e13))/0.05_dl) + 1
      ! nint : rounds its argument to the nearest whole number

      !print*, y0_act(m,z), y0_act_arr(index_mm, index_zz)

      ! without intrinsic scatter
      if (D_a == 0.) then
         !arg = (y0_act(m, z) - y0lim*1.e-4)/(sqrt(2.)*(y0lim*(1.e-4)/5.)) ! this works but slow
         !arg = (y0_act_arr(index_mm, index_zz) - y0lim*1.e-4)/(sqrt(2.)*(1.e-5)) ! test
!!!         arg = (y0_act_arr(index_mm, index_zz) - y0lim*1.e-4)/(sqrt(2.)*(y0lim*(1.e-4)/5.))
!!!         completeness_act = (1. + erf(arg))/2.

        do i = 1, Npatch
            arg = (y0_act_arr(index_mm, index_zz) - noise_act*5.)/(sqrt(2.)*noise_act)
            erfn = (1. + erf(arg))/2.
            completeness_act = completeness_act + erfn*surveydeg2_act(i)/totalArea
        end do

         ! with intrinsic scatter in y0 : log-normal
      else
         fac = 1./sqrt(2.*pi*D_a**2.)
         lnymin = -25                ! ln(1.e-10) = -23
         lnymax = 0.					! ln(1.e-2) = -4.6
         dlny = 0.1_dl
         Ny = int((lnymax - lnymin)/dlny) + 1  ! Ny = 250
!!!         Npatch = size(y0lim_arr)

         if (allocated(self)) deallocate(self)
         allocate(self(Ny),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         lny = lnymin
         do i = 1, Ny
            y0 = exp(lny)
            lny = lny + dlny
            sels = 0.
            do j = 1, Npatch
!!!               yy = y0lim_arr(j)*1.e-4
!!!               arg = (y0 - yy)/(sqrt(2.)*(yy/5.))
!!!               sel = (1. + erf(arg))/2.
!!!               sels = sels + sel/Npatch
               sigma = noise_act_arr(j)
               arg = (y0 - sigma*5.)/(sqrt(2.)*sigma)
               sel = (1. + erf(arg))/2.
               sels = sels + sel*surveydeg2_act(j)/totalArea
            end do
            self(i) = sels
         end do

         !mu = log(y0_act(m, z))
         mu = log(y0_act_arr(index_mm, index_zz))

         intsc = 0.
         lny = lnymin
         do i = 1, Ny-1
            y0 = exp(lny)
            y1 = exp(lny + dlny)
            dy = y1 - y0
            arg0 = (lny - mu)/(sqrt(2.)*D_a)
            lny = lny + dlny
            arg1 = (lny - mu)/(sqrt(2.)*D_a)
            py = 0.5_dl*fac*(self(i)/y0*exp(-arg0**2.) + self(i+1)/y1*exp(-arg1**2.))
            intsc = intsc + py*dy
            !print*, intsc
         end do
         !stop
         completeness_act = intsc
         !print*, completeness_act
         deallocate(self)
      end if

      return
   end function completeness_act

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   function y0_SO(m, z, theta, Q0)
      implicit none
      real(dl), intent(in) :: m, z, theta(:), Q0(:)
      real(dl) :: m2, A_a, B_a, y0_SO

      A_a  = cosmopar%a_act
      B_a  = cosmopar%b_act
      m2   = m*msun*cosmopar%bias_act

      y0_SO = A_a*(Eh(z)**2.)*((m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a))*Qfunc_SO(m*msun, z, theta, Q0)*relfn(m*msun, z)

      if (y0_SO < 0) y0_SO = -1.*y0_SO

      return
   end function y0_SO

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine y0_SO_allocate(y0_SO_arr, theta, Q0)
      real(dl), intent(in) :: theta(:), Q0(:)
      real(dl), allocatable :: y0_SO_arr(:,:), zz(:), mm(:)
      integer :: Nzz, Nmm, iostat, i, j
      real(dl) :: zz_min, zz_max, dzz, mm_min, mm_max, dmm

      zz_min = 0.0d0
      zz_max = 2.2
      dzz = 0.1d0
      Nzz = int((zz_max - zz_min)/dzz) + 2

      mm_min = log(1.e13)
      mm_max = log(1.e16)
      dmm = 0.05_dl
      Nmm = int((mm_max - mm_min)/dmm)

      if (allocated(y0_SO_arr)) deallocate(y0_SO_arr)
      if (allocated(zz)) deallocate(zz)
      if (allocated(mm)) deallocate(mm)
      allocate(y0_SO_arr(Nmm, Nzz), zz(Nzz), mm(Nmm), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      end if

      do i = 0, Nzz-1
         zz(i+1) = zz_min + i*dzz
      end do
      if (zz_min == 0._dl) zz(1) = zz(1) + 1.e-8

      do j = 0, Nmm-1
         mm(j+1) = mm_min + j*dmm
      end do

      do i = 1, Nzz
         do j = 1, Nmm
            y0_SO_arr(j, i) = y0_SO(exp(mm(j)), zz(i), theta, Q0)
            !print*, i, j, y0_act_arr(j, i)
         end do
         !stop
      end do

      return
   end subroutine y0_SO_allocate

   !---------------------------------------------------------------------------!

   subroutine y0_SO_deallocate
      if (allocated(y0_SO_arr)) deallocate(y0_SO_arr)
   end subroutine y0_SO_deallocate

   !---------------------------------------------------------------------------!

   function completeness_SO(m, z, noise_SO, noise_SOtest, surveydeg2_SO)
      real(dl), intent(in) :: m, z, noise_SO, noise_SOtest(:), surveydeg2_SO(:)
      real(dl) :: sel, sels, intsc, completeness_SO, erfn, totalArea
      real(dl) :: D_a, fac, mu, arg, arg0, arg1
      real(dl) :: lnymin, lnymax, dlny, lny, y0, y1, dy, py, sigma
      integer :: Ny, Npatch, i, j, iostat, index_mm, index_zz
      real(dl), allocatable :: self(:)

      D_a = cosmopar%d_act

      Npatch = size(noise_SOtest)
      totalArea = sum(surveydeg2_SO)

      index_zz = int(10.*z) + 1
      index_mm = nint((log(m)- log(1.e13))/0.05_dl) + 1
      ! nint : rounds its argument to the nearest whole number

      !print*, y0_SO(m,z), y0_SO_arr(index_mm, index_zz)

      ! without intrinsic scatter
      if (D_a == 0.) then

         completeness_SO = 0.
         do i = 1, Npatch
            arg = (y0_SO_arr(index_mm, index_zz) - noise_SO*5.)/(sqrt(2.)*noise_SO)
            erfn = (1. + erf(arg))/2.
            completeness_SO = completeness_SO + erfn*surveydeg2_SO(i)/totalArea
         end do
         !print*, y0_SO_arr(index_mm, index_zz), noise_SO

         ! with intrinsic scatter in y0 : log-normal
      else
         fac = 1./sqrt(2.*pi*D_a**2.)
         lnymin = -25                ! ln(1.e-10) = -23
         lnymax = 0.					! ln(1.e-2) = -4.6
         dlny = 0.1_dl
         Ny = int((lnymax - lnymin)/dlny) + 1  ! Ny = 250

         if (allocated(self)) deallocate(self)
         allocate(self(Ny),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         lny = lnymin
         do i = 1, Ny
            y0 = exp(lny)
            lny = lny + dlny
            sels = 0.
            do j = 1, Npatch
               sigma = noise_SOtest(j)
               arg = (y0 - sigma*5.)/(sqrt(2.)*sigma)
               sel = (1. + erf(arg))/2.
               sels = sels + sel*surveydeg2_SO(j)/totalArea
            end do
            !print*, sels
            self(i) = sels
         end do
         !stop

         !mu = log(y0_SO(m, z))
         mu = log(y0_SO_arr(index_mm, index_zz))

         intsc = 0.
         lny = lnymin
         do i = 1, Ny-1
            y0 = exp(lny)
            y1 = exp(lny + dlny)
            dy = y1 - y0
            arg0 = (lny - mu)/(sqrt(2.)*D_a)
            lny = lny + dlny
            arg1 = (lny - mu)/(sqrt(2.)*D_a)
            py = 0.5_dl*fac*(self(i)/y0*exp(-arg0**2.) + self(i+1)/y1*exp(-arg1**2.))
            intsc = intsc + py*dy
            !print*, intsc
         end do
         !stop
         completeness_SO = intsc
         !print*, completeness_act
         deallocate(self)
      end if

      return
   end function completeness_SO

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   function y0_SO_tiles(m, z, tile, theta, Q00)
      implicit none
      real(dl), intent(in) :: m, z, theta(:), Q00(:,:)
      integer, intent(in) :: tile
      real(dl) :: m2, A_a, B_a, y0_SO_tiles

      A_a  = cosmopar%a_act
      B_a  = cosmopar%b_act
      m2   = m*msun*cosmopar%bias_act

      y0_SO_tiles = A_a*(Eh(z)**2.)*((m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a))*Qfunc_SO_tiles(m*msun, z, tile, theta, Q00)*relfn(m*msun, z)

      if (y0_SO_tiles < 0) y0_SO_tiles = -1.*y0_SO_tiles

      return
   end function y0_SO_tiles

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine y0_SO_tiles_allocate(y0_SO_tiles_arr, theta, Q0)
      real(dl), allocatable :: y0_SO_tiles_arr(:,:,:), zz(:), mm(:)
      integer :: Nzz, Nmm, iostat, i, j, k, Ntt
      real(dl) :: zz_min, zz_max, dzz, mm_min, mm_max, dmm
      real(dl), intent(in) :: theta(:), Q0(:,:)

      Ntt = 346

      zz_min = 0.0d0
      zz_max = 2.2
      dzz = 0.1d0
      Nzz = int((zz_max - zz_min)/dzz) + 2

      mm_min = log(1.e13)
      mm_max = log(1.e16)
      dmm = 0.05_dl
      Nmm = int((mm_max - mm_min)/dmm)

      if (allocated(y0_SO_tiles_arr)) deallocate(y0_SO_tiles_arr)
      if (allocated(zz)) deallocate(zz)
      if (allocated(mm)) deallocate(mm)
      allocate(y0_SO_tiles_arr(Nmm, Nzz, Ntt), zz(Nzz), mm(Nmm), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      end if

      do i = 0, Nzz-1
         zz(i+1) = zz_min + i*dzz
      end do
      if (zz_min == 0._dl) zz(1) = zz(1) + 1.e-8

      do j = 0, Nmm-1
         mm(j+1) = mm_min + j*dmm
      end do

      do i = 1, Nzz
         do j = 1, Nmm
            do k = 1, Ntt
               y0_SO_tiles_arr(j, i, k) = y0_SO_tiles(exp(mm(j)), zz(i), k, theta, Q0)
               !print*, i, j, y0_act_arr(j, i)
            end do
         end do
         !stop
      end do
      return
   end subroutine y0_SO_tiles_allocate

   subroutine y0_SO_tiles_deallocate
      if (allocated(y0_SO_tiles_arr)) deallocate(y0_SO_tiles_arr)
   end subroutine y0_SO_tiles_deallocate

   !---------------------------------------------------------------------------!

   function completeness_SO_tiles(m, z, noise_SO_tiles, noise_SOfull, surveydeg2_SO_tiles, tile)
      real(dl), intent(in) :: m, z, noise_SO_tiles, noise_SOfull(:), surveydeg2_SO_tiles(:)
      integer, intent(in) :: tile
      real(dl) :: sel, sels, intsc, completeness_SO_tiles, erfn, totalArea
      real(dl) :: D_a, fac, mu, arg, arg0, arg1
      real(dl) :: lnymin, lnymax, dlny, lny, y0, y1, dy, py, sigma
      integer :: Ny, Npatch, i, j, iostat, index_mm, index_zz
      real(dl), allocatable :: self(:)

      D_a = cosmopar%d_act

      Npatch = size(noise_SOfull)
      totalArea = sum(surveydeg2_SO_tiles)

      ! index for pre-calculated array
      index_zz = int(10.*z) + 1
      index_mm = nint((log(m)- log(1.e13))/0.05_dl) + 1       ! nint : rounds its argument to the nearest whole number

      !print*, y0_SO_tiles(m,z), y0_SO_tiles_arr(index_mm, index_zz)

      ! without intrinsic scatter
      ! which means y0 is equal to the mean value
      ! only requires the selection function
      !
      if (D_a == 0.) then

         completeness_SO_tiles = 0.
         do i = 1, Npatch
            arg = (y0_SO_tiles_arr(index_mm, index_zz, tile) - noise_SO_tiles*5.)/(sqrt(2.)*noise_SO_tiles)
            erfn = (1. + erf(arg))/2.
            completeness_SO_tiles = completeness_SO_tiles + erfn*surveydeg2_SO_tiles(i)/totalArea
         end do

         !print*, y0_SO_arr(index_mm, index_zz), noise_SO

         ! with intrinsic scatter in y0 : log-normal
      else

         fac = 1./sqrt(2.*pi*D_a**2.)
         lnymin = -25                ! ln(1.e-10) = -23
         lnymax = 0.					! ln(1.e-2) = -4.6
         dlny = 0.1_dl
         Ny = int((lnymax - lnymin)/dlny) + 1  ! Ny = 250

         if (allocated(self)) deallocate(self)
         allocate(self(Ny),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         lny = lnymin
         do i = 1, Ny
            y0 = exp(lny)
            lny = lny + dlny
            sels = 0.
            do j = 1, Npatch
               sigma = noise_SOfull(j)
               arg = (y0 - sigma*5.)/(sqrt(2.)*sigma)
               sel = (1. + erf(arg))/2.
               sels = sels + sel*surveydeg2_SO_tiles(j)/totalArea
            end do
            self(i) = sels
         end do

         !mu = log(y0_act(m, z))
         mu = log(y0_SO_tiles_arr(index_mm, index_zz, tile))

         intsc = 0.
         lny = lnymin
         do i = 1, Ny-1
            y0 = exp(lny)
            y1 = exp(lny + dlny)
            dy = y1 - y0
            arg0 = (lny - mu)/(sqrt(2.)*D_a)
            lny = lny + dlny
            arg1 = (lny - mu)/(sqrt(2.)*D_a)
            py = 0.5_dl*fac*(self(i)/y0*exp(-arg0**2.) + self(i+1)/y1*exp(-arg1**2.))
            intsc = intsc + py*dy
            !print*, lny, py
         end do
         !stop
         completeness_SO_tiles = intsc
         !print*, completeness_act
         deallocate(self)
      end if

      return
   end function completeness_SO_tiles

   !---------------------------------------------------------------------------!

   ! SPT limiting mass with constant xi
   ! limiting mass as a function of redshift

   function Mlim_spt(z)
      real(dl), intent(in) :: z
      real(dl) :: xi, A_s, B_s, C_s, Mlim_spt

      A_s = cosmopar%a_spt     !5.38_dl
      B_s = cosmopar%b_spt     !1.34_dl
      C_s = cosmopar%c_spt     !0.49_dl   from DeHaan (2016)

      xi  = 5._dl              !          detection significance : 4.5 (full catalogue) and 5.0 (cosmology)
      !zeta = sqrt(xi**2. -3.) ! unbiased detection significance : 4.1 (full)               4.6 (cosmology)

      Mlim_spt = (3.e14*msun/((cosmopar%H0/100.)*cosmopar%bias_spt))*((sqrt(xi**2. - 3.)/A_s)*(Eh(0.6_dl)/Eh(z))**C_s)**(1./B_s)

      return
   end function Mlim_spt

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   ! SPT observable : unbiased significance zeta
   ! related to the detection significance xi
   ! this function is called by completeness_spt

   function zeta(m, z)
      real(dl), intent(in) :: m, z
      real(dl) :: A_s, B_s, C_s, m2, zeta, gamma_field

      A_s = cosmopar%a_spt   ! normalisation
      B_s = cosmopar%b_spt   ! slope
      C_s = cosmopar%c_spt   ! redshift evolution

      !gamma_field = 1.198
      !gamma_field = 1.504
      !gamma_field = 1.573
      gamma_field = 1.664

      m2 = m*cosmopar%bias_spt
      zeta = gamma_field*A_s*(m2/(3.e14*100./cosmopar%H0))**B_s*(Eh(z)/Eh(0.6_dl))**C_s   ! now with gamma

      return
   end function zeta

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   ! :D ! intrinsic scatter + selection function = completeness

   function completeness_spt(m, z)
      real(dl), intent(in) :: m, z
      real(dl) :: D_s, arg, completeness_spt
      !integer :: index    ! index starts from 1
      real(dl) :: xi_cut, zeta_cut
      real(dl) :: fac, lnzmin, lnzmax, dlnz, lnzeta, tempf, mu, intsc, zeta0, zeta1, dzeta, arg0, arg1, pzeta, zz
      integer :: i, Nzeta, iostat
      real(dl), allocatable :: sel(:)

      D_s = cosmopar%d_spt

      ! log() = ln and log10() = log10

      ! option 1 : pre-calculated array
      !!index = int(10.*z) + 1.
      !!if (allocated(Mlim_arr_s)) then
      !!   arg = (log(m) - log(Mlim_arr_s(index)/msun))/(sqrt(2.)*D_s)
      !!   completeness_spt = (1. + erf(arg))/2. ![0, 1]
      !!else
      !!   print*, 'error in completeness_spt'
      !!   call Mpistop()
      !!end if

      ! option 2 : simply
      !arg = (log(m) - log(Mlim_spt(z)/msun))/(sqrt(2.)*D_s)
      !completeness_spt = (1. + erf(arg))/2. ![0, 1]

      ! option3 : as in Planck
      xi_cut = 5.
      zeta_cut = sqrt(xi_cut**2. - 3.)       ! = 4.69    ln zeta_cut = 1.54

      ! without intrinsic scatter
      if (D_s == 0.) then
         arg = (zeta(m, z) - zeta_cut)/sqrt(2.)
         completeness_spt = (1. + erf(arg))/2.

      ! with intrinsic scatter
      else
         fac = 1./sqrt(2.*pi*D_s**2.)
         lnzmin = -10.              ! zeta_min = 4.53e-5
         lnzmax = 10.               ! zeta_max = 2.20e+4
         dlnz = 0.1_dl
         Nzeta = int((lnzmax - lnzmin)/dlnz)

         if (allocated(sel)) deallocate(sel)
         allocate(sel(Nzeta), stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         lnzeta = lnzmin
         do i = 1, Nzeta
            zz = exp(lnzeta)
            lnzeta = lnzeta + dlnz
            arg = (zz - zeta_cut)/sqrt(2.)
            sel(i) = (1. + erf(arg))/2.
            !print*, sel(i)
         end do
         !stop

         mu = log(zeta(m, z))
         intsc = 0.
         lnzeta = lnzmin
         do i = 1, Nzeta-1
            zeta0 = exp(lnzeta)
            zeta1 = exp(lnzeta + dlnz)
            dzeta = zeta1 - zeta0
            arg0 = (lnzeta - mu)/(sqrt(2.)*D_s)
            lnzeta = lnzeta + dlnz
            arg1 = (lnzeta - mu)/(sqrt(2.)*D_s)
            pzeta = 0.5_dl*fac*(sel(i)/zeta0*exp(-arg0**2.) + sel(i+1)/zeta1*exp(-arg1**2.))
            intsc = intsc + pzeta*dzeta
         end do

         completeness_spt = intsc
         !print*, completeness_spt

         deallocate(sel)
      endif

      return
   end function completeness_spt

   !---------------------------------------------------------------------------!

   ! for SPT mass estimate calculation

   function scatter_zeta(m, z, SZcat, xi_cat, gamma_cat)
      real(dl), intent(in) :: m, z, SZcat(:,:), xi_cat, gamma_cat
      real(dl) :: fac, D_s, mu, arg0, scatter_zeta, zeta_cat
      integer :: i

      D_s = cosmopar%d_spt
      fac = 1./sqrt(2.*pi*(D_s**2.))
      mu = log(gamma_cat*zeta(m, z))
      zeta_cat = sqrt(xi_cat**2. - 3)
      arg0 = (log(zeta_cat) - mu)/(sqrt(2.)*D_s)
      scatter_zeta = fac*exp(-arg0**2.)

      return
   end function scatter_zeta

   !---------------------------------------------------------------------------!

   ! SPT limiting mass is pre-calculated and allocated in an array : Mlim_arr_s
   ! redshift is newly assigned just for this subroutine

   subroutine Mlim_spt_allocate
      real(dl), allocatable :: zz_s(:)
      integer :: Nzz_s, i, iostat
      real(dl) :: xii, A_ss, B_ss, C_ss, zz0_s, zzmax_s, dzz_s, hh

      !hh = 0.7

      zz0_s = 0.0d0
      zzmax_s = 1.8
      dzz_s = 0.1d0
      Nzz_s = int((zzmax_s - zz0_s)/dzz_s) + 2

      if (allocated(zz_s)) deallocate(zz_s)
      if (allocated(Mlim_arr_s)) deallocate(Mlim_arr_s)
      allocate(zz_s(Nzz_s), Mlim_arr_s(Nzz_s), stat=iostat)

      if (iostat/=0) then
         print*,'allocation error'
      endif

      A_ss = cosmopar%a_spt  !5.38_dl
      B_ss = cosmopar%b_spt  !1.34_dl
      C_ss = cosmopar%c_spt  !0.49_dl
      xii  = 5._dl

      do i = 0, Nzz_s-1
         zz_s(i+1) = zz0_s + i*dzz_s ! + 0.5_dl*dzz_s
      end do
      if (zz0_s == 0._dl) zz_s(1) = zz_s(1) + 1.e-8

      ! need to check number of redshift bins !!!!!!!!!!!!!
      do i = 1, Nzz_s
         Mlim_arr_s(i) = (3.0e14*msun/((cosmopar%H0*1.e-2)*cosmopar%bias_spt))*((sqrt(xii**2. - 3.)/A_ss)*(Eh(0.6_dl)/Eh(zz_s(i)))**C_ss)**(1./B_ss)
         !Mlim_arr_s(i) = (3.0e14*msun/(hh*cosmopar%bias_spt))*((sqrt(xii**2. - 3.)/A_ss)*(Eh(0.6_dl)/Eh(zz_s(i)))**C_ss)**(1./B_ss)      ! massmass
         !print*, i, log(Mlim_arr_s(i)/msun)
      end do

      return
   end subroutine Mlim_spt_allocate

   !---------------------------------------------------------------------------!

   ! :D ! manually deallocate Mlim_arr_s
   subroutine Mlim_spt_deallocate
      if (allocated(Mlim_arr_s)) deallocate(Mlim_arr_s)
   end subroutine Mlim_spt_deallocate

   !---------------------------------------------------------------------------!

   function splintp(XA,YA,Y2A,N,X,Y)
      ! just copied and pasted : spline + interpolation
      ! given the arrays xa(1:n) and ya(1:n) of length n,
      ! which tabulate a function with the xa(i) in order,
      ! and given the array y2a(1:n),
      ! which is the output from the subroutine spline,
      ! and given a value of x,
      ! this routine returns a cubic spline interpolated value y
      INTEGER :: N
      REAL(DL) :: XA(N),YA(N),Y2A(N),X,Y,H,A,B, splintp
      INTEGER :: KLO,KHI,K
      KLO=1
      KHI=N
      1   IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
         GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) stop 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      splintp = y
      RETURN
   END function splintp

   !---------------------------------------------------------------------------!

   function erf_compl(y,sn,q)
      !completeness with error function
      REAL(dl):: y,sn,arg,erf_compl,q
      arg=(y-q*sn)/sqrt(2.)/sn
      erf_compl=(erf(arg)+1.)/2.
   end function erf_compl

   !---------------------------------------------------------------------------!

end module massobservable

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module power_sz
   use PRECISION
   use cosmology
   use Calculator_Cosmology
   implicit none
   public
   REAL(dl) :: normsig8
   REAL(dl) :: normgrowth

   contains

   !---------------------------------------------------------------------------!

   subroutine INIGROWTH
      ! normalize growth factor to 1 today
      normgrowth = 1.0_dl ! need to define this first because it is used in delta(z)
      normgrowth = 1.0/delta(0.0_dl)

   end subroutine INIGROWTH

   !---------------------------------------------------------------------------!

   function delta(z)
      ! growth factor
      real(dl), INTENT(IN) :: z
      real(dl) :: delta
      real(dl), PARAMETER :: zmax=1000.0_dl ! infinity for growth integration
      real(dl), PARAMETER :: tol=1.0d-8
      INTEGER, PARAMETER :: n=2
      INTEGER, PARAMETER :: nw = n
      real(dl) :: c(32),w(nw,9)
      real(dl) :: y(n)
      real(dl) :: a,aend
      real :: dummyr
      real(dl) :: integ,gdum,rombint
      INTEGER :: ind
      external dverk,rombint

      y(1) = 1.0_dl
      a = 1.0_dl/(1.0_dl+zmax)
      y(2) = 0.0_dl
      aend=1.0_dl/(1.0_dl+z)
      ind = 1
      ! dverk is defined in CAMB subroutines.f90
      ! dummyr is not assigned any value

      if (cosmopar%gamma.eq.-1.0) then
         call dverk(dummyr,n,ddelta,a,y,aend,tol,ind,c,nw,w)
         if (ind.ne.3) write(*,*) 'Problem in dverk',ind
         delta = normgrowth*y(2)
         !       write(*,*) a,normgrowth*(y(1)/a-y(2)/a/a)
      else
         integ=rombint(growint,a,aend,tol)
         gdum=aend*exp(integ)
         delta=gdum*normgrowth
      endif
      RETURN
   end function delta

   !---------------------------------------------------------------------------!

   subroutine ddelta(dummyr,n,a,y,yprime)
      ! note dummyr is completely irrelevant but required for CAMB dverk
      real(dl), INTENT(IN) :: a
      real, INTENT(IN) :: dummyr
      INTEGER, INTENT(IN) :: n
      real(dl) :: y(n),yprime(n)
      real(dl) :: z,w,u,del

      z=1.0_dl/a-1.0_dl
      w=cosmopar%w0+cosmopar%w1*(1.0_dl-a)
      u = y(1)
      del = y(2)

      yprime(1) = -1.5_dl*(1.0_dl+omegak(z)/3.0_dl-w*omegade(z))*u/a+1.5_dl*omegam(z)*del/a/a
      !yprime(1) = -1.5_dl*(1.0_dl+omegak(z)/3.0_dl-w*omegade(z))*u/a+1.5_dl*0.3*del/a/a !!!!!!!!!!!!!! massmass
      yprime(2) = u
      RETURN
   end subroutine ddelta

   !---------------------------------------------------------------------------!

   function growint(a)
      real(dl), INTENT(IN) :: a
      real(dl) :: growint,z
      z=1.0_dl/a-1.0_dl
      growint=(omegam(z)**cosmopar%gamma-1.0_dl)/a
      !growint=(0.3**cosmopar%gamma-1.0_dl)/a !!!!!!!!!!!!1 massmass!!!!!!!!!!!
      return
   end function growint

   !---------------------------------------------------------------------------!

end module power_sz

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module massfunction
   USE PRECISION
   USE CONSTANTS_sz
   USE POWER_SZ
   USE massobservable

   implicit none
   TYPE MASSPAR
   REAL(dl) :: Amf,Bmf,epsmf,dso,pmf
   INTEGER :: psind
   END TYPE MASSPAR
   TYPE (masspar) :: massfnpar

   integer :: scatter_mlim_switch, completeness_act_switch
   ! switches for ACT are added

   contains

   !---------------------------------------------------------------------------!

   function dndlnM(z,M,g)
      !mass function: Jenkins, Tinker (default), Watson
      !mass in units of h^-1 M_sun
      real(dl), INTENT(in) :: z,M,g
      real(dl) :: dndlnM
      real(dl) :: R,rhom0,rhom
      real(dl) :: dMdR,sR,fJen
      real(dl) :: fTink
      real(dl) :: fWat,p,q,CDelta,Gamma,ddz,z0
      real(dl) :: alpha,dsoz
      real(dl) :: del(1:9),par_aa(1:9),par_a(1:9),par_b(1:9),par_c(1:9)
      real(dl) :: der_aa(1:9),der_a(1:9),der_b(1:9),der_c(1:9)
      real(dl) :: par1,par2,par3,par4
      integer  :: total,i

      !rhom0=0.3*rhocrit0 !!!!!!!!!!!!!!!!!!!!!!!! massmass
      rhom0=cosmopar%omegam*rhocrit0
      rhom=rhom0*(1.0+z)**3

      total=9
      del(1)=200
      del(2)=300
      del(3)=400
      del(4)=600
      del(5)=800
      del(6)=1200
      del(7)=1600
      del(8)=2400
      del(9)=3200
      par_aa(1)=0.186
      par_aa(2)=0.200
      par_aa(3)=0.212
      par_aa(4)=0.218
      par_aa(5)=0.248
      par_aa(6)=0.255
      par_aa(7)=0.260
      par_aa(8)=0.260
      par_aa(9)=0.260
      par_a(1)=1.47
      par_a(2)=1.52
      par_a(3)=1.56
      par_a(4)=1.61
      par_a(5)=1.87
      par_a(6)=2.13
      par_a(7)=2.30
      par_a(8)=2.53
      par_a(9)=2.66
      par_b(1)=2.57
      par_b(2)=2.25
      par_b(3)=2.05
      par_b(4)=1.87
      par_b(5)=1.59
      par_b(6)=1.51
      par_b(7)=1.46
      par_b(8)=1.44
      par_b(9)=1.41
      par_c(1)=1.19
      par_c(2)=1.27
      par_c(3)=1.34
      par_c(4)=1.45
      par_c(5)=1.58
      par_c(6)=1.80
      par_c(7)=1.97
      par_c(8)=2.24
      par_c(9)=2.44
      der_aa(1)=0.00
      der_aa(2)=0.50
      der_aa(3)=-1.56
      der_aa(4)=3.05
      der_aa(5)=-2.95
      der_aa(6)=1.07
      der_aa(7)=-0.71
      der_aa(8)=0.21
      der_aa(9)=0.00
      der_a(1)=0.00
      der_a(2)=1.19
      der_a(3)=-6.34
      der_a(4)=21.36
      der_a(5)=-10.95
      der_a(6)=2.59
      der_a(7)=-0.85
      der_a(8)=-2.07
      der_a(9)=0.00
      der_b(1)=0.00
      der_b(2)=-1.08
      der_b(3)=12.61
      der_b(4)=-20.96
      der_b(5)=24.08
      der_b(6)=-6.64
      der_b(7)=3.84
      der_b(8)=-2.09
      der_b(9)=0.00
      der_c(1)=0.00
      der_c(2)=0.94
      der_c(3)=-0.43
      der_c(4)=4.61
      der_c(5)=0.01
      der_c(6)=1.21
      der_c(7)=1.43
      der_c(8)=0.33
      der_c(9)=0.00

      do i=1,9
         del(i)=log10(del(i))
      enddo

      ! radius of shell of mass M with density rhom0
      R = (0.75_dl*M/Pi/rhom0)**(1._dl/3._dl) ! R in units of h^-1 Mpc
      ! check if powerspectrum for sigma(R) is actually calculated to large
      ! enough k
      ! maybe put this later outside dndlnM, for lowest mass limit and test ahead

      ! dM/dR
      dMdR = 3*M/R

      sR = cosmopar%sigmaR%Value(R)
      !    write(*,*) cosmopar%omegam,sR
      ! g = delta(z)
      ! Parameters from Jenkins et al. for LCDM (Appendix) and SO(324)
      ! Amf = 0.316 ; Bmf = 0.67 ; epsmf = 3.82
      ! Note SO(324) is mass in mass density units or 0.3*324=97.2

      dsoz = massfnpar%dso/Omegam(z)
      !dsoz = massfnpar%dso/0.3_dl !!!!!!!!!!!!!!!!!!!!!!! massmass

      if (massfnpar%psind==1) then  !Jenkins et al.
         fJen = massfnpar%Amf*exp(-abs(-dlog(g*sR)+massfnpar%Bmf)**massfnpar%epsmf)
         dndlnM = -rhom0*fJen*cosmopar%sigmaR%Derivative(R)/dMdR/sR
         elseif (massfnpar%psind==2) then  !Tinker et al.
            call SPLINTNR(del,par_aa,der_aa,total,log10(dsoz),par1)
            call SPLINTNR(del,par_a,der_a,total,log10(dsoz),par2)
            call SPLINTNR(del,par_b,der_b,total,log10(dsoz),par3)
            call SPLINTNR(del,par_c,der_c,total,log10(dsoz),par4)

            alpha=10**(-((0.75_dl/dlog10(dsoz/75.0_dl))**1.2_dl))
            massfnpar%Amf = par1*((1.0_dl+z)**(-0.14_dl))
            massfnpar%epsmf = par2*((1.0_dl+z)**(-0.06_dl))
            massfnpar%Bmf = par3*((1.0_dl+z)**(-alpha))
            massfnpar%pmf = par4

            fTink = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)

            dndlnM = -rhom0*fTink*cosmopar%sigmaR%Derivative(R)/dMdR/sR

         elseif (massfnpar%psind==3) then !Watson et al. 2013

            !FOF

            par1=0.282_dl !A
            par2=2.163_dl !alpha
            par3=1.406_dl !beta
            par4=1.210_dl !gamma

            !CPMSO
            !           par1=0.287_dl !A
            !           par2=2.234_dl !alpha
            !           par3=1.478_dl !beta
            !           par4=1.318_dl !gamma

            !!$           !redshift dependence:
            !par1=Omegam(z)*(0.990_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !A
            !           par1=Omegam(z)*(1.097_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !gaz
            !           par2=Omegam(z)*(5.907_dl*(1.0_dl+z)**(-3.058_dl)+2.349_dl) !bz
            !           par3=Omegam(z)*(3.136_dl*(1.0_dl+z)**(-3.599_dl)+2.344_dl) !az

            massfnpar%Amf = par1
            massfnpar%epsmf = par2
            massfnpar%Bmf = par3
            massfnpar%pmf = par4

            fWat = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)

            !Delta dependence:
            dsoz=massfnpar%dso/Omegam(z)
            ddz=-0.456_dl*Omegam(z)-0.139
            CDelta=exp(0.023_dl*(dsoz/178._dl-1.0_dl))*0.947
            p=0.072_dl
            q=2.130_dl

            Gamma = CDelta*(dsoz/178._dl)**ddz*exp(p*(1.0_dl-dsoz/178._dl)/(sR*g)**q)

            fwat = fwat*Gamma

            dndlnM = -rhom0*fWat*cosmopar%sigmaR%Derivative(R)/dMdR/sR
         else
            write(*,*) 'Invalid mass function indicator: ',massfnpar%psind
            stop
         end if

      RETURN
   end function dndlnM

   !---------------------------------------------------------------------------!

   ! :D ! y0lim is passed (just copied from above)
   function dndlnMy(z, y0lim, M, g)
      !mass function: Jenkins, Tinker (default), Watson
      !mass in units of h^-1 M_sun
      real(dl), INTENT(in) :: z, y0lim, M, g
      real(dl) :: dndlnMy
      real(dl) :: R,rhom0,rhom
      real(dl) :: dMdR,sR,fJen
      real(dl) :: fTink
      real(dl) :: fWat,p,q,CDelta,Gamma,ddz,z0
      real(dl) :: alpha,dsoz
      real(dl) :: del(1:9),par_aa(1:9),par_a(1:9),par_b(1:9),par_c(1:9)
      real(dl) :: der_aa(1:9),der_a(1:9),der_b(1:9),der_c(1:9)
      real(dl) :: par1,par2,par3,par4
      integer  :: total,i

      rhom0=cosmopar%omegam*rhocrit0
      rhom=rhom0*(1.0+z)**3

      total=9
      del(1)=200
      del(2)=300
      del(3)=400
      del(4)=600
      del(5)=800
      del(6)=1200
      del(7)=1600
      del(8)=2400
      del(9)=3200
      par_aa(1)=0.186
      par_aa(2)=0.200
      par_aa(3)=0.212
      par_aa(4)=0.218
      par_aa(5)=0.248
      par_aa(6)=0.255
      par_aa(7)=0.260
      par_aa(8)=0.260
      par_aa(9)=0.260
      par_a(1)=1.47
      par_a(2)=1.52
      par_a(3)=1.56
      par_a(4)=1.61
      par_a(5)=1.87
      par_a(6)=2.13
      par_a(7)=2.30
      par_a(8)=2.53
      par_a(9)=2.66
      par_b(1)=2.57
      par_b(2)=2.25
      par_b(3)=2.05
      par_b(4)=1.87
      par_b(5)=1.59
      par_b(6)=1.51
      par_b(7)=1.46
      par_b(8)=1.44
      par_b(9)=1.41
      par_c(1)=1.19
      par_c(2)=1.27
      par_c(3)=1.34
      par_c(4)=1.45
      par_c(5)=1.58
      par_c(6)=1.80
      par_c(7)=1.97
      par_c(8)=2.24
      par_c(9)=2.44
      der_aa(1)=0.00
      der_aa(2)=0.50
      der_aa(3)=-1.56
      der_aa(4)=3.05
      der_aa(5)=-2.95
      der_aa(6)=1.07
      der_aa(7)=-0.71
      der_aa(8)=0.21
      der_aa(9)=0.00
      der_a(1)=0.00
      der_a(2)=1.19
      der_a(3)=-6.34
      der_a(4)=21.36
      der_a(5)=-10.95
      der_a(6)=2.59
      der_a(7)=-0.85
      der_a(8)=-2.07
      der_a(9)=0.00
      der_b(1)=0.00
      der_b(2)=-1.08
      der_b(3)=12.61
      der_b(4)=-20.96
      der_b(5)=24.08
      der_b(6)=-6.64
      der_b(7)=3.84
      der_b(8)=-2.09
      der_b(9)=0.00
      der_c(1)=0.00
      der_c(2)=0.94
      der_c(3)=-0.43
      der_c(4)=4.61
      der_c(5)=0.01
      der_c(6)=1.21
      der_c(7)=1.43
      der_c(8)=0.33
      der_c(9)=0.00

      do i=1,9
         del(i)=log10(del(i))
      enddo

      ! radius of shell of mass M with density rhom0
      R = (0.75_dl*M/Pi/rhom0)**(1._dl/3._dl) ! R in units of h^-1 Mpc
      ! check if powerspectrum for sigma(R) is actually calculated to large
      ! enough k
      ! maybe put this later outside dndlnM, for lowest mass limit and test ahead

      ! dM/dR
      dMdR = 3*M/R

      sR = cosmopar%sigmaR%Value(R)
      !    write(*,*) cosmopar%omegam,sR
      ! g = delta(z)
      ! Parameters from Jenkins et al. for LCDM (Appendix) and SO(324)
      ! Amf = 0.316 ; Bmf = 0.67 ; epsmf = 3.82
      ! Note SO(324) is mass in mass density units or 0.3*324=97.2

      dsoz = massfnpar%dso/Omegam(z)
      if (massfnpar%psind==1) then  !Jenkins et al.
         fJen = massfnpar%Amf*exp(-abs(-dlog(g*sR)+massfnpar%Bmf)**massfnpar%epsmf)
         dndlnMy = -rhom0*fJen*cosmopar%sigmaR%Derivative(R)/dMdR/sR
      elseif (massfnpar%psind==2) then  !Tinker et al.
         call SPLINTNR(del,par_aa,der_aa,total,log10(dsoz),par1)
         call SPLINTNR(del,par_a,der_a,total,log10(dsoz),par2)
         call SPLINTNR(del,par_b,der_b,total,log10(dsoz),par3)
         call SPLINTNR(del,par_c,der_c,total,log10(dsoz),par4)

         alpha=10**(-((0.75_dl/dlog10(dsoz/75.0_dl))**1.2_dl))
         massfnpar%Amf = par1*((1.0_dl+z)**(-0.14_dl))
         massfnpar%epsmf = par2*((1.0_dl+z)**(-0.06_dl))
         massfnpar%Bmf = par3*((1.0_dl+z)**(-alpha))
         massfnpar%pmf = par4

         fTink = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)

         dndlnMy = -rhom0*fTink*cosmopar%sigmaR%Derivative(R)/dMdR/sR

      elseif (massfnpar%psind==3) then !Watson et al. 2013

         !FOF

         par1=0.282_dl !A
         par2=2.163_dl !alpha
         par3=1.406_dl !beta
         par4=1.210_dl !gamma

         !CPMSO
         !           par1=0.287_dl !A
         !           par2=2.234_dl !alpha
         !           par3=1.478_dl !beta
         !           par4=1.318_dl !gamma

         !!$           !redshift dependence:
         !par1=Omegam(z)*(0.990_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !A
         !           par1=Omegam(z)*(1.097_dl*(1.0_dl+z)**(-3.216_dl)+0.074_dl) !gaz
         !           par2=Omegam(z)*(5.907_dl*(1.0_dl+z)**(-3.058_dl)+2.349_dl) !bz
         !           par3=Omegam(z)*(3.136_dl*(1.0_dl+z)**(-3.599_dl)+2.344_dl) !az

         massfnpar%Amf = par1
         massfnpar%epsmf = par2
         massfnpar%Bmf = par3
         massfnpar%pmf = par4

         fWat = massfnpar%Amf*((g*sR/massfnpar%Bmf)**(-massfnpar%epsmf)+1.0_dl)*exp(-massfnpar%pmf/sR/sR/g/g)

         !Delta dependence:
         dsoz=massfnpar%dso/Omegam(z)
         ddz=-0.456_dl*Omegam(z)-0.139
         CDelta=exp(0.023_dl*(dsoz/178._dl-1.0_dl))*0.947
         p=0.072_dl
         q=2.130_dl

         Gamma = CDelta*(dsoz/178._dl)**ddz*exp(p*(1.0_dl-dsoz/178._dl)/(sR*g)**q)

         fwat = fwat*Gamma

         dndlnMy = -rhom0*fWat*cosmopar%sigmaR%Derivative(R)/dMdR/sR
      else
         write(*,*) 'Invalid mass function indicator: ',massfnpar%psind
         stop
      end if

      RETURN
   end function dndlnMy

   !---------------------------------------------------------------------------!

   ! :D ! integration over ACT mass with y0lim_const
   function intmassfn_act(z, lnMmin, lnMmax)
      real(dl), intent(in) :: z, lnMmin, lnMmax
      real(dl) :: g, intmassfn_act, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM) ! M=[1.e13, 1.e16] :lnM=[29.9, 36.8] :Nm = 138
      intmassfn_act = 0.
      mass = lnMmin

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)

         if (scatter_mlim_switch == 1) then       ! without scatter (+ ylim map without scatter)
            sum = 0.5_dl*(dndlnM(z, mass0, g) + dndlnM(z, mass1, g))
         elseif (scatter_mlim_switch == 2) then   ! with scatter
            sum = 0.5_dl*(dndlnM(z, mass0, g)*mlim_scatter_act(mass0, z) + dndlnM(z, mass1, g)*mlim_scatter_act(mass1, z))
         endif

         intmassfn_act = intmassfn_act + sum*dlnM
      end do

      return
   end function intmassfn_act

   !---------------------------------------------------------------------------!

   ! :D ! integration over ACT mass with limiting y0 map with scatter
   ! this is main
   !function intmassfn_act_ymap_scatter(z, lnMmin, lnMmax, y0lim, y0lim_arr)
   function intmassfn_act_ymap_scatter(z, lnMmin, lnMmax, noise_act, noise_act_arr, surveydeg2_act)
      !real(dl), intent(in) :: z, lnMmin, lnMmax, y0lim, y0lim_arr(:)
      real(dl), intent(in) :: z, lnMmin, lnMmax, noise_act, noise_act_arr(:), surveydeg2_act(:)
      real(dl) :: g, intmassfn_act_ymap_scatter, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM) ! M=[1.e13, 1.e16] :lnM=[29.9, 36.8] :Nm = 138
      intmassfn_act_ymap_scatter = 0.
      mass = lnMmin

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)

!         if (completeness_act_switch == 1) then
!            sum = 0.5_dl*(dndlnM(z, mass0, g)*completeness_act(mass0, z, y0lim, y0lim_arr) + dndlnM(z, mass1, g)*completeness_act(mass1, z, y0lim, y0lim_arr))
!         else
!            sum = 0.5_dl*(dndlnM(z, mass0, g)*mlim_scatter_act_ymap(mass0, z, y0lim, y0lim_arr) + dndlnM(z, mass1, g)*mlim_scatter_act_ymap(mass1, z, y0lim, y0lim_arr))
!         end if

         sum = 0.5_dl*(dndlnM(z, mass0, g)*completeness_act(mass0, z, noise_act, noise_act_arr, surveydeg2_act) + dndlnM(z, mass1, g)*completeness_act(mass1, z, noise_act, noise_act_arr, surveydeg2_act))

         intmassfn_act_ymap_scatter = intmassfn_act_ymap_scatter + sum*dlnM
         !if (z == 0.5_dl) then
         !print*, mass0, mlim_scatter_act_ymap(mass0, z, y0lim, y0lim_arr)
         !print*, mass0, completeness_act(mass0, z, y0lim, y0lim_arr)
         !end if
      end do

      return
   end function intmassfn_act_ymap_scatter

   !---------------------------------------------------------------------------!

   function intmassfn_act2D(z, y0, lnMmin, lnMmax)
      ! :D ! integration over ACT mass in 2D case and here passing y0
      ! not using it now
      real(dl), intent(in) :: z, y0, lnMmin, lnMmax
      real(dl) :: g, intmassfn_act2D, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM)
      intmassfn_act2D = 0.
      mass = lnMmin

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)
         sum = 0.5_dl*(dndlnMy(z, y0, mass0, g) + dndlnMy(z, y0, mass1, g))
         intmassfn_act2D = intmassfn_act2D + sum*dlnM
      end do

      return
   end function intmassfn_act2D

   !---------------------------------------------------------------------------!

   function intmassfn_act_mass(z, lnMmin, lnMmax, SZcat)
      real(dl), intent(in) :: z, lnMmin, lnMmax, SZcat(:,:)
      real(dl) :: g, intmassfn_act_mass, mass, mass0, mass1, normal
      integer :: i, Nm, j
      real, parameter :: dlnM = 0.05_dl ! stepsize for m
      real(dl) :: z_cat, y0_cat, y0err_cat, sum0, sum1, mass_cat, mass_unc

      do j = 1, 182

         z_cat = SZcat(j,1)
         y0_cat = SZcat(j,4)
         y0err_cat = SZcat(j,6)
         mass_unc = SZcat(j,7)
         mass_cat = SZcat(j,8)
         !print*, mass_cat

         y0_cat = y0_cat*1.e-4
         y0err_cat = y0err_cat*1.e-4

         g  = delta(z_cat) ! growth factor
         Nm = int((lnMmax - lnMmin)/dlnM)
         intmassfn_act_mass = 0.
         normal = 0.
         mass = lnMmin


         do i = 1, Nm-1
            sum0 = 0.
            sum1 = 0.
            mass0 = exp(mass + dlnM*(i-1))
            mass1 = exp(mass + dlnM*i)

            if (scatter_mlim_switch == 1) then ! no mlim scatter
               sum0 = 0.5_dl*(dndlnM(z_cat, mass0, g) + dndlnM(z_cat, mass1, g)) ! don't need it actually
            elseif (scatter_mlim_switch == 2) then

               !sum0 = 0.5_dl*(mass0*dndlnM(z_cat, mass0, g)*scatter_y0(mass0, z_cat, SZcat, y0_cat, y0err_cat) &
               !	& + mass1*dndlnM(z_cat, mass1, g)*scatter_y0(mass1, z_cat, SZcat, y0_cat, y0err_cat))

               !sum1 = 0.5_dl*(dndlnM(z_cat, mass0, g)*scatter_y0(mass0, z_cat, SZcat, y0_cat, y0err_cat) &
               !  & + dndlnM(z_cat, mass1, g)*scatter_y0(mass1, z_cat, SZcat, y0_cat, y0err_cat))

               sum0 = 0.5_dl*(dndlnM(z_cat, mass0, g)*scatter_y0(mass0, z_cat, SZcat, y0_cat, y0err_cat) &
               & + dndlnM(z_cat, mass1, g)*scatter_y0(mass1, z_cat, SZcat, y0_cat, y0err_cat))

               sum1 = 0.5_dl*(dndlnM(z_cat, mass0, g)*scatter_y0(mass0, z_cat, SZcat, y0_cat, y0err_cat)/mass0 &
               & + dndlnM(z_cat, mass1, g)*scatter_y0(mass1, z_cat, SZcat, y0_cat, y0err_cat)/mass1)

            endif

            intmassfn_act_mass = intmassfn_act_mass + sum0*dlnM
            normal = normal + sum1*dlnM

            !print*, log(mass0), theta_act(mass0*msun, z_cat), Qfunc(mass0*msun, z_cat)

         end do

         !print*, dlog(intmassfn_act_mass/normal), dlog(mass_cat*1.e14) ! UPP mass estimates
         !print*, dlog(Mlim_act_ymap(z_cat, y0_cat*1.e4)/msun), dlog(mass_unc*1.e14) ! recovering mass estimate from a simple inversion

      end do
      !stop

      return
   end function intmassfn_act_mass

   !---------------------------------------------------------------------------!

   ! :D ! integration over SPT mass with constant xi
   function intmassfn_spt(z, lnMmin, lnMmax)
      real(dl), intent(in) :: z, lnMmin, lnMmax
      real(dl) :: g, intmassfn_spt, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM) ! M=[1.e13, 1.e16] :lnM=[29.9, 36.8] :Nm = 138
      intmassfn_spt = 0.
      mass = lnMmin

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)

         if (scatter_mlim_switch == 1) then
            sum = 0.5_dl*(dndlnM(z, mass0, g) + dndlnM(z, mass1, g))
         elseif (scatter_mlim_switch == 2) then
            sum = 0.5_dl*(dndlnM(z, mass0, g)*completeness_spt(mass0, z) + dndlnM(z, mass1, g)*completeness_spt(mass1, z))
         endif
         !print*, z, dlog(mass0), completeness_spt(mass0,z)

         intmassfn_spt = intmassfn_spt + sum*dlnM

      end do
      !print*, intmassfn_spt

      return
   end function intmassfn_spt

   !---------------------------------------------------------------------------!

   function intmassfn_spt_mass(z, lnMmin, lnMmax, SZcat)
      real(dl), intent(in) :: z, lnMmin, lnMmax, SZcat(:,:)
      real(dl) :: g, intmassfn_spt_mass, sum0, sum1, mass, mass0, mass1
      integer :: i, Nm, j
      real, parameter :: dlnM = 0.05_dl ! stepsize for m
      real(dl) :: z_cat, xi_cat, normal, mass_cat, gamma_cat

      do j = 1, 343

         z_cat = SZcat(j,1)
         xi_cat = SZcat(j,3)
         gamma_cat = SZcat(j,4)
         mass_cat = SZcat(j,5)
         !print*, mass_cat, dlog(mass_cat)

         g  = delta(z_cat) ! growth factor
         Nm = int((lnMmax - lnMmin)/dlnM) ! usually lnMmin = 33 and lnMmax = 36
         intmassfn_spt_mass = 0.
         normal = 0.
         mass = lnMmin

         do i = 1, Nm-1
            sum0 = 0.
            sum1 = 0.
            mass0 = exp(mass + dlnM*(i-1))
            mass1 = exp(mass + dlnM*i)

            if (scatter_mlim_switch == 1) then
               sum0 = 0.5_dl*(dndlnM(z_cat, mass0, g) + dndlnM(z_cat, mass1, g))
            elseif (scatter_mlim_switch == 2) then

               sum0 = 0.5_dl*(dndlnM(z_cat, mass0, g)*scatter_zeta(mass0, z_cat, SZcat, xi_cat, gamma_cat) &
               & + dndlnM(z_cat, mass1, g)*scatter_zeta(mass1, z_cat, SZcat, xi_cat, gamma_cat))

               sum1 = 0.5_dl*(dndlnM(z_cat, mass0, g)*scatter_zeta(mass0, z_cat, SZcat, xi_cat, gamma_cat)/mass0 &
               & + dndlnM(z_cat, mass1, g)*scatter_zeta(mass1, z_cat, SZcat, xi_cat, gamma_cat)/mass1)

            endif

            intmassfn_spt_mass = intmassfn_spt_mass + sum0*dlnM
            normal = normal + sum1*dlnM

         end do
         print*, intmassfn_spt_mass/normal/1.e14, mass_cat

      end do
      stop

      return
   end function intmassfn_spt_mass

   !---------------------------------------------------------------------------!

   ! SO trial

   function intmassfn_SO(z, lnMmin, lnMmax, noise_SO, noise_SOtest, surveydeg2_SO)
      real(dl), intent(in) :: z, lnMmin, lnMmax, noise_SO, noise_SOtest(:), surveydeg2_SO(:)
      real(dl) :: g, intmassfn_SO, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM) ! M=[1.e13, 1.e16] :lnM=[29.9, 36.8] :Nm = 138
      intmassfn_SO = 0.
      mass = lnMmin

      !print*, surveydeg2_SO
      !stop

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)

         sum = 0.5_dl*(dndlnM(z, mass0, g)*completeness_SO(mass0, z, noise_SO, noise_SOtest, surveydeg2_SO) + dndlnM(z, mass1, g)*completeness_SO(mass1, z, noise_SO, noise_SOtest, surveydeg2_SO))

         intmassfn_SO = intmassfn_SO + sum*dlnM

         !if (z == 0.5_dl) then
         !print*, mass0, mlim_scatter_act_ymap(mass0, z, y0lim, y0lim_arr)
         !print*, mass0, completeness_SO(mass0, z, noise_SO, noise_SOtest, surveydeg2_SO)
         !end if
      end do

      return
   end function intmassfn_SO

   !---------------------------------------------------------------------------!

   ! SO tile trial

   function intmassfn_SO_tiles(z, lnMmin, lnMmax, noise_SO_tiles, noise_SOfull, surveydeg2_SO_tiles, tile)
      real(dl), intent(in) :: z, lnMmin, lnMmax, noise_SO_tiles, noise_SOfull(:), surveydeg2_SO_tiles(:)
      integer, intent(in) :: tile
      real(dl) :: g, intmassfn_SO_tiles, sum, mass, mass0, mass1
      integer :: i, Nm
      real, parameter :: dlnM = 0.05_dl ! stepsize for m

      g  = delta(z) ! growth factor
      Nm = int((lnMmax - lnMmin)/dlnM) ! M=[1.e13, 1.e16] :lnM=[29.9, 36.8] :Nm = 138
      intmassfn_SO_tiles = 0.
      mass = lnMmin

      do i = 1, Nm-1
         sum = 0.
         mass0 = exp(mass + dlnM*(i-1))
         mass1 = exp(mass + dlnM*i)

         sum = 0.5_dl*(dndlnM(z, mass0, g)*completeness_SO_tiles(mass0, z, noise_SO_tiles, noise_SOfull, surveydeg2_SO_tiles, tile) + dndlnM(z, mass1, g)*completeness_SO_tiles(mass1, z, noise_SO_tiles, noise_SOfull, surveydeg2_SO_tiles, tile))

         intmassfn_SO_tiles = intmassfn_SO_tiles + sum*dlnM

         !if (z == 0.5_dl) then
         !print*, mass0, mlim_scatter_act_ymap(mass0, z, y0lim, y0lim_arr)
         !print*, mass0, completeness_SO(mass0, z, noise_SO, noise_SOtest)
         !end if
      end do

      return
   end function intmassfn_SO_tiles

   !---------------------------------------------------------------------------!

   SUBROUTINE SPLINTNR(XA,YA,Y2A,N,X,Y)
      INTEGER :: N
      REAL(DL) :: XA(N),YA(N),Y2A(N),X,Y,H,A,B
      INTEGER :: KLO,KHI,K
      KLO=1
      KHI=N
      1   IF (KHI-KLO.GT.1) THEN
      K=(KHI+KLO)/2
      IF(XA(K).GT.X)THEN
      KHI=K
      ELSE
      KLO=K
      ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) stop 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
   END SUBROUTINE SPLINTNR

   !---------------------------------------------------------------------------!

end module massfunction

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module numbercounts
USE PRECISION
USE CONSTANTS_sz
USE POWER_SZ
USE MASSOBSERVABLE
USE MASSFUNCTION

implicit none
integer :: mlim_switch, survey_switch, y0bin ! :D !
contains

   !---------------------------------------------------------------------------!

   function next_z(zi,binz)
      real(dl):: zi,binz,next_z,dzi,hr
      !bins in redshifts are defined with higher resolution for z<0.2
      hr=0.2
      if (zi <hr) then
      dzi=1.e-3
      else if ((zi >=hr) .and. (zi <=1.)) then
      dzi=1.e-2
      else
      dzi=binz
      endif
      next_z=zi+dzi
      return
   end function next_z

   !---------------------------------------------------------------------------!

   ! :D ! dNdz for ACT with y0lim_const and for SPT with constant xi

   function dNdz(z, m11, m12)
      real(dl), intent(in) :: z, m11, m12
      real(dl) :: dNdz

      if ((survey_switch == 2) .or. (survey_switch == 3)) then     ! 2 = ACT, 3 = Planck + ACT
         dNdz = surveypar%deg2_a * dVdzdO(z) * intmassfn_act(z, m11, m12)
      elseif ((survey_switch == 4) .or. (survey_switch == 5)) then ! 4 = SPT, 5 = Planck + SPT
         dNdz = surveypar%deg2_s * dVdzdO(z) * intmassfn_spt(z, m11, m12)
      endif

      return
   end function dNdz

   !---------------------------------------------------------------------------!

   ! :D ! dNdz for ACT with ylim map without scatter
   ! should be integrated over skypatch
   ! surveypar%deg2_a should be surveypar%deg2_a/(number of skypatch)

   !function dNdz_act_ymap(z, m11, m12, y0lim_arr)
   function dNdz_act_ymap(z, m11, m12, noise_act_arr, surveydeg2_act)
      !real(dl), intent(in) :: z, m11, m12, y0lim_arr(:)
      real(dl), intent(in) :: z, m11, m12, noise_act_arr(:), surveydeg2_act(:)
      real(dl) :: mnew, sum, dNdz_act_ymap
      integer :: i, znew
!!!      integer, parameter :: Npatch = 48
      integer :: Npatch

      Npatch = size(surveydeg2_act)

      dNdz_act_ymap = 0.

      do i = 1, Npatch
         sum = 0.
!         if (scatter_mlim_switch == 1) then        ! without scattering
!            if (allocated(Mlim_arr_a)) then
!               znew = int(10.*z) + 1
!               mnew = log(Mlim_arr_a(znew, i)/msun)
!               sum  = surveypar%deg2_a * dVdzdO(z) * intmassfn_act(z, mnew, m12)
!            !print*, z, mnew
!            else
!               print*, 'error in scatter_act'
!               call Mpistop()
!            end if

!         elseif (scatter_mlim_switch == 2) then    ! with scattering
!             sum = surveypar%deg2_a * dVdzdO(z) * intmassfn_act_ymap_scatter(z, m11, m12, y0lim_arr(i), y0lim_arr)
!         end if

         sum = surveydeg2_act(i) * dVdzdO(z) * intmassfn_act_ymap_scatter(z, m11, m12, noise_act_arr(i), noise_act_arr, surveydeg2_act)
         dNdz_act_ymap = dNdz_act_ymap + sum

      end do

      return
   end function dNdz_act_ymap

   !---------------------------------------------------------------------------!

   function dNdz_act2D(z, y0, m11, m12)
      ! :D ! dNdz for 2D ACT case : dN(z, y0) - not using for now
      real(dl), intent(in) :: z, y0, m11, m12
      real(dl) :: dNdz_act2D
      dNdz_act2D = surveypar%deg2_a * dVdzdO(z) * intmassfn_act2D(z, y0, m11, m12)
      return
   end function dNdz_act2D

   !---------------------------------------------------------------------------!

   function dNdz_SO(z, m11, m12, noise_SOtest, surveydeg2_SO)
      real(dl), intent(in) :: z, m11, m12, noise_SOtest(:), surveydeg2_SO(:)
      real(dl) :: sum, dNdz_SO
      integer :: i
      integer :: Npatch

      Npatch = size(surveydeg2_SO)

      dNdz_SO = 0.
      do i = 1, Npatch
         sum = 0.
         sum = surveydeg2_SO(i) * dVdzdO(z) * intmassfn_SO(z, m11, m12, noise_SOtest(i), noise_SOtest, surveydeg2_SO)
         dNdz_SO = dNdz_SO + sum
         !print*, 'sky patch number = ', i, 'dNdz = ', dNdz_SO
      end do

      return
   end function dNdz_SO

   !---------------------------------------------------------------------------!

   function dNdz_SO_tiles(z, m11, m12, noise_SOfull, surveydeg2_SO_tiles, tilenames)
      real(dl), intent(in) :: z, m11, m12, noise_SOfull(:), surveydeg2_SO_tiles(:)
      integer, intent(in) :: tilenames(:)
      real(dl) :: sum, dNdz_SO_tiles
      integer :: i
      integer :: Npatch

      Npatch = size(surveydeg2_SO_tiles)

      dNdz_SO_tiles = 0.
      do i = 1, Npatch
         sum = 0.
         sum = surveydeg2_SO_tiles(i) * dVdzdO(z) * intmassfn_SO_tiles(z, m11, m12, noise_SOfull(i), noise_SOfull, surveydeg2_SO_tiles, tilenames(i))
         dNdz_SO_tiles = dNdz_SO_tiles + sum
      end do

      return
   end function dNdz_SO_tiles

   !---------------------------------------------------------------------------!

   ! :D ! computing number counts in redshift

   subroutine deltaN(Nz, z, delN)
      real(dl), intent(in) :: z(:)
      integer, intent(in) :: Nz
      real(dl) :: m11, m12
      real(dl), allocatable :: delN(:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i

      if (allocated(delN)) deallocate(delN)
      allocate(delN(Nz))

      m12 = log(1.e16) ! mass upper limit

      ! integral over redshift
      do i = 1, Nz

         ! set the mass lower limit

         ! without scattering in Mlim - so integral just starts from Mlim
         if (scatter_mlim_switch == 1) then

            ! constant Mlim for all z
            if (mlim_switch == 1) then
               m11 = log(2.35e14)  ! Battye & Weller (2003) lnMmin

            ! Mlim as a function of redshift with constant y0lim or xi
            elseif (mlim_switch == 2) then
               if ((survey_switch == 2) .or. (survey_switch == 3)) then        ! 2 = ACT, 3 = Planck + ACT
                  m11 = log(Mlim_act(z(i))/msun)
               elseif ((survey_switch == 4) .or. (survey_switch == 5)) then    ! 4 = SPT, 5 = Planck + SPT
                  m11 = log(Mlim_spt(z(i))/msun)
               endif
            endif

         ! including scattering in Mlim
         ! becomes gaussian distribution with a central value of Mlim
         ! integral starts from very low value of M
         elseif (scatter_mlim_switch == 2) then
            m11 = log(1.e13)
         endif

         delN(i) = 0.5_dl*(dNdz(z(i)-0.5_dl*dz, m11, m12) + dNdz(z(i)+0.5_dl*dz, m11, m12))*dz

      end do

      return
   end subroutine deltaN

   !---------------------------------------------------------------------------!

   ! :D ! computing number counts in redshift for ACT with ylim map

   !subroutine deltaN_act_ymap(Nz, z, delN, y0lim_arr)
   subroutine deltaN_act_ymap(Nz, z, delN, noise_act_arr, surveydeg2_act)
      !real(dl), intent(in) :: z(:), y0lim_arr(:)
      real(dl), intent(in) :: z(:), noise_act_arr(:), surveydeg2_act(:)
      integer, intent(in) :: Nz
      real(dl) :: m11, m12
      real(dl), allocatable :: delN(:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i

      if (allocated(delN)) deallocate(delN)
      allocate(delN(Nz))

      mlim_switch = 2    ! varied Mlim = 2             does it matter?????
      survey_switch = 3  ! ACT alone = 2  Planck + ACT = 3  ! now you have to do it manually ! this needs to be sorted out !!!!!

      m11 = log(1.0e13)
      m12 = log(1.0e16)

      ! integral over redshift
      do i = 1, Nz

         !delN(i) = 0.5_dl*(dNdz_act_ymap(z(i)-0.5_dl*dz, m11, m12, y0lim_arr) + dNdz_act_ymap(z(i)+0.5_dl*dz, m11, m12, y0lim_arr))*dz
         delN(i) = 0.5_dl*(dNdz_act_ymap(z(i)-0.5_dl*dz, m11, m12, noise_act_arr, surveydeg2_act) + dNdz_act_ymap(z(i)+0.5_dl*dz, m11, m12, noise_act_arr, surveydeg2_act))*dz

      end do

      return
   end subroutine deltaN_act_ymap

   !---------------------------------------------------------------------------!

   subroutine deltaN2D(Nz, z, Ny0, y0, delN2D)
      ! computing number counts in redshift and y0 for ACT2D  - not using it now
      real(dl), intent(in) :: z(:), y0(:)
      integer, intent(in) :: Nz, Ny0
      real(dl) :: m11, m12
      real(dl), allocatable :: delN2D(:,:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i, j

      if (allocated(delN2D)) deallocate(delN2D)
      allocate(delN2D(Nz,Ny0))

      if (mlim_switch == 1) then  ! constant Mlim case
         m11 = log(2.35e14) ! Battye & Weller (2003) lnMmin
         m12 = log(1.0e16)  ! lnMmax
      elseif (mlim_switch == 2) then ! Mlim as a function of y0 and z
         do j = 1, Ny0 ! number of y0 bins
            do i = 1, Nz
               m11 = log(Mlim_act_ymap(z(i), y0(j))/msun)
               if (j < Ny0) then
                  m12 = log(Mlim_act_ymap(z(i), y0(j+1))/msun)
               elseif (j == Ny0) then ! last bin
                  m12 = log(1.0e16)
               endif
               delN2D(i, j) = 0.5_dl*(dNdz_act2D(z(i)-0.5_dl*dz, y0(j), m11, m12) + dNdz_act2D(z(i)+0.5_dl*dz, y0(j), m11, m12))*dz
            end do
         end do
      end if

      return
   end subroutine deltaN2D

   !---------------------------------------------------------------------------!

   function dNdz_mass(z, m11, m12, SZcat)
      ! for calculating mass estimates
      real(dl), intent(in) :: z, m11, m12, SZcat(:,:)
      real(dl) :: dNdz_mass

      if ((survey_switch == 2) .or. (survey_switch == 3)) then     ! 2 = ACT, 3 = Planck + ACT
         dNdz_mass = surveypar%deg2_a * dVdzdO(z) * intmassfn_act_mass(z, m11, m12, SZcat)
      elseif ((survey_switch == 4) .or. (survey_switch == 5)) then ! 4 = SPT, 5 = Planck + SPT
         dNdz_mass = surveypar%deg2_s * dVdzdO(z) * intmassfn_spt_mass(z, m11, m12, SZcat)
      endif

      return
   end function dNdz_mass

   !---------------------------------------------------------------------------!

   subroutine deltaN_mass(Nz, z, delN, SZcat)
      ! for calculating mass estimates
      real(dl), intent(in) :: z(:), SZcat(:,:)
      integer, intent(in) :: Nz
      real(dl) :: m11, m12
      real(dl), allocatable :: delN(:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i

      if (allocated(delN)) deallocate(delN)
      allocate(delN(Nz))

      m12 = log(1.e16) ! mass upper limit

      ! integral over redshift
      do i = 1, Nz

         ! set the mass lower limit

         ! without scattering in Mlim - so integral just starts from Mlim
         if (scatter_mlim_switch == 1) then

            ! constant Mlim for all z
            if (mlim_switch == 1) then
               m11 = log(2.35e14)  ! Battye & Weller (2003) lnMmin

            ! Mlim as a function of redshift with constant ylim or xi
            elseif (mlim_switch == 2) then
               if ((survey_switch == 2) .or. (survey_switch == 3)) then        ! 2 = ACT, 3 = Planck + ACT
                  m11 = log(Mlim_act(z(i))/msun)
               elseif ((survey_switch == 4) .or. (survey_switch == 5)) then    ! 4 = SPT, 5 = Planck + SPT
                   m11 = log(Mlim_spt(z(i))/msun)
               endif
            endif

         ! including scattering in Mlim
         ! becomes gaussian distribution with a central value of Mlim
         ! integral starts from very low value of M
         elseif (scatter_mlim_switch == 2) then
            m11 = log(1.e13)
         endif

         delN(i) = 0.5_dl*(dNdz_mass(z(i)-0.5_dl*dz, m11, m12, SZcat) + dNdz_mass(z(i)+0.5_dl*dz, m11, m12, SZcat))*dz

      end do

      return
   end subroutine deltaN_mass

   !---------------------------------------------------------------------------!

   subroutine deltaN_SO(Nz, z, delN, noise_SOtest, surveydeg2_SO)
      real(dl), intent(in) :: z(:), noise_SOtest(:), surveydeg2_SO(:)
      integer, intent(in) :: Nz
      real(dl) :: m11, m12
      real(dl), allocatable :: delN(:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i

      if (allocated(delN)) deallocate(delN)
      allocate(delN(Nz))

      survey_switch = 6

      m11 = log(1.0e13)
      m12 = log(1.0e16)

      ! integral over redshift
      do i = 1, Nz
         delN(i) = 0.5_dl*(dNdz_SO(z(i)-0.5_dl*dz, m11, m12, noise_SOtest, surveydeg2_SO) + dNdz_SO(z(i)+0.5_dl*dz, m11, m12, noise_SOtest, surveydeg2_SO))*dz
         !print*, 'redshift = ', i/10., 'delN = ', delN(i)
         !print*, i/10., delN(i)
      end do

      return
   end subroutine deltaN_SO

   !---------------------------------------------------------------------------!

   subroutine deltaN_SO_tiles(Nz, z, delN, noise_SOfull, surveydeg2_SO_tiles, tilenames)
      real(dl), intent(in) :: z(:), noise_SOfull(:), surveydeg2_SO_tiles(:)
      integer, intent(in) :: Nz, tilenames(:)
      real(dl) :: m11, m12
      real(dl), allocatable :: delN(:)
      real(dl), parameter :: dz = 0.1_dl
      integer :: i

      if (allocated(delN)) deallocate(delN)
      allocate(delN(Nz))

      survey_switch = 7

      m11 = log(1.0e13)
      m12 = log(1.0e16)

      ! integral over redshift
      do i = 1, Nz
         delN(i) = 0.5_dl*(dNdz_SO_tiles(z(i)-0.5_dl*dz, m11, m12, noise_SOfull, surveydeg2_SO_tiles, tilenames) + dNdz_SO_tiles(z(i)+0.5_dl*dz, m11, m12, noise_SOfull, surveydeg2_SO_tiles, tilenames))*dz
      end do

      return
   end subroutine deltaN_SO_tiles

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !--- using grid ------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine deltaN_SOtest(z, Nz, delN, skyfracs_SOtest, noise_SOtest, theta, Q0)
      !subroutine computing numbercounts in z
      real(dl), intent(in) :: z(:), skyfracs_SOtest(:), noise_SOtest(:), theta(:), Q0(:)
      integer, intent(in) :: Nz
      real(dl), allocatable :: grid(:,:), steps_z(:), steps_m(:)
      real(sp), allocatable :: completeness(:,:)
      real(dl) :: delN(:)
      real(dl) :: lnmmin, lnmmax, lnm, binz, min_z, max_z, zi, z1, z2
      integer :: i, j, nsteps_m, nsteps_z, iostat, npatches
      real(dl), parameter :: dlnm = 0.05_dl

      npatches = size(skyfracs_SOtest)

      lnmmin = log(1.e13) ! = 29.9 this is natural log!!!
      lnmmax = log(1.e16) ! = 36.8

      binz = z(3) - z(2)

      !size of grid for mass
      nsteps_m = dint((lnmmax - lnmmin)/dlnm)

      !size of grid for z
      min_z = z(1) - 0.5_dl*binz
      max_z = z(Nz) + 0.5_dl*binz

      zi = min_z
      if (zi < 0) zi = 0.

      nsteps_z = 0
      do while (zi < max_z)
        !zi = next_z(zi, binz)
        zi = zi + binz
        nsteps_z = nsteps_z + 1
      enddo

      allocate(steps_m(nsteps_m), steps_z(nsteps_z), stat=iostat)
      if (iostat/=0) then
        print*, 'allocation error'
      endif

      !grid for mass
      lnm = lnmmin
      do i = 1, nsteps_m
        steps_m(i) = lnm + dlnm/2.
        lnm = lnm + dlnm
      enddo

      !grid for z
      zi = min_z
      if (zi < 0) zi = 0.

      do i = 1, nsteps_z
        steps_z(i) = zi
        !zi = next_z(zi, binz)
        zi = zi + binz
      enddo

      if (steps_z(1) == 0.) steps_z(1) = 1.e-5

      allocate(grid(nsteps_m, nsteps_z), completeness(nsteps_m, nsteps_z), stat=iostat)
      if (iostat/=0) then
          print*, 'allocation error'
          stop
      endif
      completeness(:,:) = 0.

      call grid_C_SOtest(steps_z, steps_m, skyfracs_SOtest, noise_SOtest, completeness, theta, Q0)
      !tabulate completeness for m and z

      grid(:,:) = 0.
      call get_grid_SO_test(nsteps_z, nsteps_m, steps_z, steps_m, grid)
      !tabulate number count for m and z

      delN(:) = 0.
      do i = 1, Nz
          z1 = z(i) - 0.5_dl*binz
          z2 = z(i) + 0.5_dl*binz
          call integrate(grid, completeness, steps_z, nsteps_z, nsteps_m, z1, z2, dlnm, delN(i))
          !integrate on mass and z within bins accounting for completeness
      enddo

      deallocate(grid, steps_m, steps_z, completeness, stat=iostat)
      if (iostat/=0) then
          print*, 'allocation error'
          stop
      endif

      return
   end subroutine deltaN_SOtest


   subroutine grid_C_SOtest(steps_z, steps_m, skyfracs_SOtest, noise_SOtest, completeness, theta, Q0)
      !tabulate completeness for SO 1D likelihood
      real(dl), intent(in) :: steps_z(:), steps_m(:), skyfracs_SOtest(:), noise_SOtest(:), theta(:), Q0(:)
      integer :: nsteps_z, nsteps_m, npatches, Ny0
      real(sp) :: completeness(:,:)
      real(dl) :: D_a, totalArea, mp, zp, y0p, snrcut, noise, c
      real(dl) :: sum0, fac, lny, lnymin, lnymax, dlny, yy, mu
      real(dl) :: intsc, y0, y1, dy, arg0, arg1, py
      integer :: i, ii, jj, k, iostat
      real(dl), allocatable :: tempf(:)

      D_a = cosmopar%d_act

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)
      npatches = size(noise_SOtest)
      totalArea = sum(skyfracs_SOtest)

      if (D_a == 0.) then
        ! no intrinsic scatter
        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = dexp(steps_m(ii))
                zp = steps_z(jj)
                y0p = y0_SO(mp, zp, theta, Q0)

                completeness(ii,jj) = 0.

                do i = 1, npatches
                    snrcut = 5.
                    noise = noise_SOtest(i)
                    c = erf_compl(y0p, noise, snrcut)
                    completeness(ii, jj) = completeness(ii, jj) + c*skyfracs_SOtest(i)/totalArea
                enddo
            enddo
        enddo

      else
        ! intrinsic scatter in y0-m relation
        fac = 1./sqrt(2.*pi*D_a**2)
        lnymin = -25.    !ln(1.e-10) = -23
        lnymax = 0.      !ln(1.e-2) = -4.6
        dlny = 0.1_dl

        Ny0 = int((lnymax - lnymin)/dlny) + 1

        if (allocated(tempf)) deallocate(tempf)
        allocate(tempf(Ny0), stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
        endif

        tempf(:) = 0.

        lny = lnymin
        do k = 1, Ny0
            yy = dexp(lny)
            lny = lny + dlny
            sum0 = 0.
            do i = 1, npatches
                snrcut = 5.
                noise = noise_SOtest(i)
                c = erf_compl(yy, noise, snrcut)
                sum0 = sum0 + c*skyfracs_SOtest(i)/totalArea
            enddo
            tempf(k) = sum0
        enddo

        do jj=1,nsteps_z
            do ii=1,nsteps_m
                mp = dexp(steps_m(ii))
                zp = steps_z(jj)
                mu = dlog(y0_SO(mp, zp, theta, Q0))

                intsc = 0.
                lny = lnymin
                do k = 1, Ny0-1
                    y0 = dexp(lny)
                    y1 = dexp(lny+dlny)
                    dy = y1 - y0
                    arg0 = ((lny - mu)/(sqrt(2.)*D_a))
                    lny = lny + dlny
                    arg1 = ((lny - mu)/(sqrt(2.)*D_a))
                    py = (tempf(k)*fac/y0*dexp(-arg0**2) + tempf(k+1)*fac/y1*dexp(-arg1**2))*0.5
                    intsc = intsc + py*dy
                enddo

                if (intsc > 1.) intsc = 1.
                if (intsc < 0.) intsc = 0.
                completeness(ii,jj) = intsc

            enddo
        enddo

        deallocate(tempf)

      endif

   end subroutine grid_C_SOtest

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine deltaN_SO_fullmap(z, Nz, delN, skyfracs_SOfull, noise_SOfull, tilenames, theta, Q00)
      !subroutine computing numbercounts in z and xi
      real(dl), intent(in) :: z(:), skyfracs_SOfull(:), noise_SOfull(:), theta(:), Q00(:,:)
      integer, intent(in) :: Nz, tilenames(:)
      real(dl), allocatable :: grid(:,:), steps_z(:), steps_m(:)
      real(sp), allocatable :: completeness(:,:)
      real(dl) :: delN(:)
      real(dl) :: lnmmin, lnmmax, lnm, binz, min_z, max_z, zi, z1, z2
      integer :: i, j, nsteps_m, nsteps_z, iostat, npatches
      real(dl), parameter :: dlnm = 0.05_dl

      npatches = size(skyfracs_SOfull)

      lnmmin = log(1.e13) ! = 29.9 this is natural log!!!
      lnmmax = log(1.e16) ! = 36.8

      binz = z(3) - z(2)

      !size of grid for mass
      nsteps_m = dint((lnmmax - lnmmin)/dlnm)

      !size of grid for z
      min_z = z(1) - 0.5_dl*binz
      max_z = z(Nz) + 0.5_dl*binz

      zi = min_z
      if (zi < 0) zi = 0.

      nsteps_z = 0
      do while (zi < max_z)
        !zi = next_z(zi, binz)
        zi = zi + binz
        nsteps_z = nsteps_z + 1
      enddo

      allocate(steps_m(nsteps_m), steps_z(nsteps_z), stat=iostat)
      if (iostat/=0) then
        print*, 'allocation error'
      endif

      !grid for mass
      lnm = lnmmin
      do i = 1, nsteps_m
        steps_m(i) = lnm + dlnm/2.
        lnm = lnm + dlnm
      enddo

      !grid for z
      zi = min_z
      if (zi < 0) zi = 0.

      do i = 1, nsteps_z
        steps_z(i) = zi
        !zi = next_z(zi, binz)
        zi = zi + binz
      enddo

      if (steps_z(1) == 0.) steps_z(1) = 1.e-5

      allocate(grid(nsteps_m, nsteps_z), completeness(nsteps_m, nsteps_z), stat=iostat)
      if (iostat/=0) then
          print*, 'allocation error'
          stop
      endif

      call precalc_y0_allocate(steps_z, steps_m, tilenames, precalc_y0, theta, Q00)

      completeness(:,:) = 0.
      call grid_C_SO(steps_z, steps_m, skyfracs_SOfull, noise_SOfull, completeness, tilenames)
      !tabulate completeness for m and z

      grid(:,:) = 0.
      call get_grid_SO(nsteps_z, nsteps_m, steps_z, steps_m, grid)
      !tabulate number count for m and z

      delN(:) = 0.
      do i = 1, Nz
          z1 = z(i) - 0.5_dl*binz
          z2 = z(i) + 0.5_dl*binz
          call integrate(grid, completeness, steps_z, nsteps_z, nsteps_m, z1, z2, dlnm, delN(i))
          !integrate on mass and z within bins accounting for completeness
      enddo

      call precalc_y0_deallocate

      deallocate(grid, steps_m, steps_z, completeness, stat=iostat)
      if (iostat/=0) then
          print*, 'allocation error'
          stop
      endif

      return
   end subroutine deltaN_SO_fullmap


   subroutine precalc_y0_allocate(steps_z, steps_m, tilenames, precalc_y0, theta, Q00)
      real(dl), intent(in) :: steps_z(:), steps_m(:), theta(:), Q00(:,:)
      real(dl), allocatable :: precalc_y0(:,:,:)
      integer, intent(in) :: tilenames(:)
      integer :: nsteps_z, nsteps_m, nsteps_t, iostat, i, j, k, kp
      real(dl) :: mp, zp

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)
      nsteps_t = size(tilenames)

      if (allocated(precalc_y0)) deallocate(precalc_y0)
      allocate(precalc_y0(nsteps_m, nsteps_z, nsteps_t), stat=iostat)
      if (iostat/=0) then
        print*, 'allocation error'
      endif

      do i = 1, nsteps_m
        do j = 1, nsteps_z
          do k = 1, nsteps_t
             mp = dexp(steps_m(i))
             zp = steps_z(j)
             kp = tilenames(k)

             precalc_y0(i, j, k) = y0_SO_tiles(mp, zp, kp, theta, Q00)
          enddo
        enddo
      enddo

      return
   end subroutine precalc_y0_allocate


   subroutine precalc_y0_deallocate
    if (allocated(precalc_y0)) deallocate(precalc_y0)
   end subroutine precalc_y0_deallocate


   subroutine grid_C_SO(steps_z, steps_m, skyfracs_SOfull, noise_SOfull, completeness, tile)
      !tabulate completeness for SO 1D likelihood
      real(dl), intent(in) :: steps_z(:), steps_m(:), skyfracs_SOfull(:), noise_SOfull(:)
      integer, intent(in) :: tile(:)
      integer :: nsteps_z, nsteps_m, npatches, Ny0
      real(sp) :: completeness(:,:)
      real(dl) :: D_a, totalArea, y0p, snrcut, noise, c, sum0, fac, lny, lnymin, lnymax, dlny, yy, mu
      real(dl) :: intsc, y0, y1, dy, arg0, arg1, py, mp, zp
      integer :: i, ii, jj, k, iostat
      real(dl), allocatable :: tempf(:)

      D_a = cosmopar%d_act

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)
      npatches = size(noise_SOfull)
      totalArea = sum(skyfracs_SOfull)

      if (D_a == 0.) then
        ! no intrinsic scatter
        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                do i = 1, npatches
                    y0p = precalc_y0(ii, jj, tile(i))
                    snrcut = 5.
                    noise = noise_SOfull(i)
                    c = erf_compl(y0p, noise, snrcut)
                    completeness(ii, jj) = completeness(ii, jj) + c*skyfracs_SOfull(i)/totalArea
                enddo
            enddo
        enddo

      else
        ! intrinsic scatter in y0-m relation
        fac = 1./sqrt(2.*pi*D_a**2)
        lnymin = -25.    !ln(1.e-10) = -23
        lnymax = 0.      !ln(1.e-2) = -4.6
        dlny = 0.1_dl

        Ny0 = int((lnymax - lnymin)/dlny) + 1

        if (allocated(tempf)) deallocate(tempf)
        allocate(tempf(Ny0), stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
        endif

        tempf(:) = 0.
        lny = lnymin
        do k = 1, Ny0
            yy = dexp(lny)
            lny = lny + dlny
            sum0 = 0.
            do i = 1, npatches
                snrcut = 5.
                noise = noise_SOfull(i)
                c = erf_compl(yy, noise, snrcut)
                sum0 = sum0 + c*skyfracs_SOfull(i)/totalArea
            enddo
            tempf(k) = sum0
        enddo

        do jj=1,nsteps_z
            do ii=1,nsteps_m
                do i = 1, npatches
                    mu = dlog(precalc_y0(ii, jj, tile(i)))
                    intsc = 0.
                    lny = lnymin
                    do k = 1, Ny0-1
                        y0 = dexp(lny)
                        y1 = dexp(lny+dlny)
                        dy = y1 - y0
                        arg0 = ((lny - mu)/(sqrt(2.)*D_a))
                        lny = lny + dlny
                        arg1 = ((lny - mu)/(sqrt(2.)*D_a))
                        py = (tempf(k)*fac/y0*dexp(-arg0**2) + tempf(k+1)*fac/y1*dexp(-arg1**2))*0.5
                        intsc = intsc + py*dy
                    enddo

                 enddo

                if (intsc > 1.) intsc = 1.
                if (intsc < 0.) intsc = 0.
                completeness(ii,jj) = intsc

            enddo
        enddo

        deallocate(tempf)

      endif

   end subroutine grid_C_SO

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine deltaN_act(z, Nz, delN, y0lims)
      !subroutine computing numbercounts in z and xi
      real(dl), intent(in) :: z(:), y0lims(:)
      integer, intent(in) :: Nz
      real(dl), allocatable :: grid(:,:), steps_z(:), steps_m(:)
      real(sp), allocatable :: completeness(:,:)
      real(dl) :: delN(:)
      real(dl) :: lnmmin, lnmmax, lnm, binz, min_z, max_z, zi, z1, z2
      integer :: i, j, nsteps_m, nsteps_z, iostat, npatches
      real(dl), parameter :: dlnm = 0.05_dl

      npatches = size(y0lims)

      lnmmin = log(1.e13) ! = 29.9 this is natural log!!!
      lnmmax = log(1.e16) ! = 36.8

      binz = z(3) - z(2)

      !size of grid for mass
      nsteps_m = dint((lnmmax - lnmmin)/dlnm)

      !size of grid for z
      min_z = z(1) - 0.5_dl*binz
      max_z = z(Nz) + 0.5_dl*binz

      zi = min_z
      if (zi < 0) zi = 0.

      nsteps_z = 0
      do while (zi < max_z)
        !zi = next_z(zi, binz)
        zi = zi + binz
        nsteps_z = nsteps_z + 1
      enddo

      allocate(steps_m(nsteps_m), steps_z(nsteps_z), stat=iostat)
      if (iostat/=0) then
        print*, 'allocation error'
      endif

      !grid for mass
      lnm = lnmmin
      do i = 1, nsteps_m
        steps_m(i) = lnm + dlnm/2.
        lnm = lnm + dlnm
      enddo

      !grid for z
      zi = min_z
      if (zi < 0) zi = 0.

      do i = 1, nsteps_z
        steps_z(i) = zi
        !zi = next_z(zi, binz)
        zi = zi + binz
      enddo

      if (steps_z(1) == 0.) steps_z(1) = 1.e-5

      !select case(switch)
      !case(1) !1D likelihood

        allocate(grid(nsteps_m, nsteps_z), completeness(nsteps_m, nsteps_z), stat=iostat)
        if (iostat/=0) then
            print*, 'allocation error'
            stop
        endif

        completeness(:,:) = 0.
        call grid_C_act(steps_z, steps_m, y0lims, completeness) !surveydeg2_act,
        !tabulate completeness for m and z

        grid(:,:) = 0.
        call get_grid_act(nsteps_z, nsteps_m, steps_z, steps_m, grid)
        !tabulate number count for m and z

        delN(:) = 0.
        do i = 1, Nz
            z1 = z(i) - 0.5_dl*binz
            z2 = z(i) + 0.5_dl*binz
            call integrate(grid, completeness, steps_z, nsteps_z, nsteps_m, z1, z2, dlnm, delN(i))
            !integrate on mass and z within bins accounting for completeness
        enddo

        deallocate(grid, steps_m, steps_z, completeness, stat=iostat)
        if (iostat/=0) then
            print*, 'allocation error'
            stop
        endif

      return
   end subroutine deltaN_act


   subroutine grid_C_act(steps_z, steps_m, y0lims, completeness)
      !tabulate completeness for ACT 1D likelihood
      real(dl), intent(in) :: steps_z(:), steps_m(:), y0lims(:)
      integer :: nsteps_z, nsteps_m, npatches, Ny0
      real(sp) :: completeness(:,:)
      real(dl) :: D_a, totalArea, y0p, snrcut, noise, c, sum0, fac, lny, lnymin, lnymax, dlny, yy, mu
      real(dl) :: intsc, y0, y1, dy, arg0, arg1, py, mp, zp
      integer :: i, ii, jj, k, iostat
      real(dl), allocatable :: tempf(:)

      D_a = cosmopar%d_act

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)
      npatches = size(y0lims)

      if (D_a == 0.) then
        ! no intrinsic scatter
        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = dexp(steps_m(ii))
                zp = steps_z(jj)
                y0p = y0_act(mp, zp)

                do i = 1, npatches
                    snrcut = 5.
                    noise = y0lims(i)*1.e-4/snrcut
                    c = erf_compl(y0p, noise, snrcut)
                    !completeness(ii, jj) = completeness(ii, jj) + c*surveydeg2_act(i)/totalArea
                    completeness(ii, jj) = completeness(ii, jj) + c/npatches
                enddo
            enddo
        enddo

      else
        ! intrinsic scatter in y0-m relation
        fac = 1./sqrt(2.*pi*D_a**2)
        lnymin = -25.    !ln(1.e-10) = -23
        lnymax = 0.      !ln(1.e-2) = -4.6
        dlny = 0.1_dl

        Ny0 = int((lnymax - lnymin)/dlny) + 1

        if (allocated(tempf)) deallocate(tempf)
        allocate(tempf(Ny0), stat=iostat)
        if (iostat/=0) then
            print*,'allocation error'
        endif

        tempf(:) = 0.

        lny = lnymin
            do k = 1, Ny0
                yy = dexp(lny)
                lny = lny + dlny
                sum0 = 0.
                do i = 1, npatches
                    snrcut = 5.
                    noise = y0lims(i)*1.e-4/snrcut
                    c = erf_compl(yy, noise, snrcut)
                    !sum0 = sum0 + c*surveydeg2_act(i)/totalArea
                    sum0 = sum0 + c/npatches
                enddo
                tempf(k) = sum0
            enddo

        do jj=1,nsteps_z
            do ii=1,nsteps_m
                mp = dexp(steps_m(ii))
                zp = steps_z(jj)
                mu = dlog(y0_act(mp, zp))

                intsc = 0.
                lny = lnymin
                do k = 1, Ny0-1
                    y0 = dexp(lny)
                    y1 = dexp(lny+dlny)
                    dy = y1 - y0
                    arg0 = ((lny - mu)/(sqrt(2.)*D_a))
                    lny = lny + dlny
                    arg1 = ((lny - mu)/(sqrt(2.)*D_a))
                    py = (tempf(k)*fac/y0*dexp(-arg0**2) + tempf(k+1)*fac/y1*dexp(-arg1**2))*0.5
                    intsc = intsc + py*dy
                enddo

                if (intsc > 1.) intsc = 1.
                if (intsc < 0.) intsc = 0.
                completeness(ii,jj) = intsc

            enddo
        enddo

        deallocate(tempf)

      endif

   end subroutine grid_C_act

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   subroutine deltaN_spt(z, Nz, logxi, Nxi, delN, delN2D, switch)
      !subroutine computing numbercounts in z and xi
      real(dl), intent(in) :: z(:), logxi(:)
      integer, intent(in) :: Nz, Nxi, switch
      real(dl), allocatable :: grid(:,:), steps_z(:), steps_m(:)
      real(sp), allocatable :: completeness(:,:), completeness_2d(:,:,:)
      real(dl) :: delN(:), delN2D(:,:)
      real(dl) :: lnmmin, lnmmax, lnm, binz, binx, min_z, max_z, zi, z1, z2
      integer :: i, j, nsteps_m, nsteps_z, iostat
      real(dl), parameter :: dlnm = 0.05_dl


      lnmmin = log(1.e13) ! = 29.9 this is natural log!!!
      lnmmax = log(1.e16) ! = 36.8

      binz = z(3) - z(2)
      binx = logxi(2) - logxi(1)

      !size of grid for mass
      nsteps_m = dint((lnmmax - lnmmin)/dlnm)

      !size of grid for z
      min_z = z(1) - 0.5_dl*binz
      max_z = z(Nz) + 0.5_dl*binz

      zi = min_z
      if (zi < 0) zi = 0.

      nsteps_z = 0
      do while (zi < max_z)
        !zi = next_z(zi, binz)
        zi = zi + binz
        nsteps_z = nsteps_z + 1
      enddo

      allocate(steps_m(nsteps_m), steps_z(nsteps_z), stat=iostat)
      if (iostat/=0) then
        print*, 'allocation error'
      endif

      !grid for mass
      lnm = lnmmin
      do i = 1, nsteps_m
        steps_m(i) = lnm + dlnm/2.
        lnm = lnm + dlnm
      enddo

      !grid for z
      zi = min_z
      if (zi < 0) zi = 0.

      do i = 1, nsteps_z
        steps_z(i) = zi
        !zi = next_z(zi, binz)
        zi = zi + binz
      enddo

      if (steps_z(1) == 0.) steps_z(1) = 1.e-5

      select case(switch)
      case(1) !1D likelihood

          allocate(grid(nsteps_m, nsteps_z), completeness(nsteps_m, nsteps_z), stat=iostat)
          if (iostat/=0) then
            print*, 'allocation error'
            stop
          endif

          completeness(:,:) = 0.
          call grid_C_spt(steps_z, steps_m, completeness)
          !tabulate completeness for m and z

          grid(:,:) = 0.
          call get_grid_spt(nsteps_z, nsteps_m, steps_z, steps_m, grid)
          !tabulate number count for m and z
          delN(:) = 0.
          do i = 1, Nz
            z1 = z(i) - 0.5_dl*binz
            z2 = z(i) + 0.5_dl*binz

            call integrate(grid, completeness, steps_z, nsteps_z, nsteps_m, z1, z2, dlnm, delN(i))
            !integrate on mass and z within bins accounting for completeness

          enddo

          deallocate(grid, steps_m, steps_z, completeness, stat=iostat)
          if (iostat/=0) then
            print*, 'allocation error'
            stop
          endif

      case(2) !2D likelihood

          allocate(grid(nsteps_m, nsteps_z), completeness_2d(nsteps_m, nsteps_z, Nxi+1), stat=iostat)
          if (iostat/=0) then
            print*, 'allocation error'
            stop
          endif

          completeness_2d(:,:,:) = 0.
          do j = 1, Nxi+1
            call grid_C_2d_spt(steps_z, steps_m, completeness_2d, logxi, j)
            !tabulate completeness for m and z
          enddo

          grid(:,:) = 0.
          call get_grid_spt(nsteps_z, nsteps_m, steps_z, steps_m, grid)
          !tabulate number count for m and z

          do j = 1, Nxi+1
            do i = 1, Nz+1
              z1 = z(i) - 0.5_dl*binz
              z2 = z(i) + 0.5_dl*binz

              call integrate_2d_spt(grid, completeness_2d, steps_z, nsteps_z, nsteps_m, z1, z2, j, dlnm, delN2D(i,j))
              !integrate on mass and z within bins accounting for completeness

            enddo
          enddo

          deallocate(grid, steps_m, steps_z, completeness_2d, stat=iostat)
          if (iostat/=0) then
            print*, 'allocation error'
            stop
          endif

      end select

      return
   end subroutine deltaN_spt


   subroutine integrate(grid, compl, steps_z, nsteps_z, nsteps_m, z1, z2, dlnm, sum0)
      !integration for SPT 1D likelihood
      real(dl), intent(in) :: grid(:,:), steps_z(:), z1, z2, dlnm
      real(sp), intent(in) :: compl(:,:)
      real(dl) :: sum0, psum, c1, c2, f1, f2, x1, x2, test(nsteps_z)
      integer :: nsteps_z, nsteps_m, ii, jj, p(1), j1, j2

      !find range i1->i2 for integration in z
      test = abs(steps_z - z1)
      p = minloc(test)
      j1 = p(1)
      test = abs(steps_z - z2)
      p = minloc(test)
      j2 = p(1)

      sum0 = 0._dl
      psum = 0._dl
      do jj = j1, j2-1
        do ii = 1, nsteps_m
            x1 = steps_z(jj)
            x2 = steps_z(jj+1)
            f1 = grid(ii, jj)
            f2 = grid(ii, jj+1)
            c1 = compl(ii, jj)
            c2 = compl(ii, jj+1)
            psum = psum + 0.5*(f1*c1 + f2*c2)*(x2-x1)*dlnm
        enddo

      enddo
      sum0 = sum0 + psum

   end subroutine integrate


   subroutine integrate_2d_spt(grid, compl, steps_z, nsteps_z, nsteps_m, z1, z2, binx, dlnm, sum0)
      !integration for SPT 2D likelihood
      real(dl), intent(in) :: grid(:,:), steps_z(:), z1, z2, dlnm
      real(sp), intent(in) :: compl(:,:,:)
      real(dl) :: sum0, psum, c1, c2, f1, f2, x1, x2, test(nsteps_z)
      integer :: nsteps_z, nsteps_m, binx, ii, jj, p(1), j1, j2

      !find range i1->i2 for integration in z
      test = abs(steps_z - z1)
      p = minloc(test)
      j1 = p(1)
      test = abs(steps_z - z2)
      p = minloc(test)
      j2 = p(1)

      sum0 = 0._dl
      psum = 0._dl
      do jj = j1, j2-1
        do ii = 1, nsteps_m
            x1 = steps_z(jj)
            x2 = steps_z(jj+1)
            f1 = grid(ii, jj)
            f2 = grid(ii, jj+1)
            c1 = compl(ii, jj, binx)
            c2 = compl(ii, jj+1, binx)
            psum = psum + 0.5*(f1*c1 + f2*c2)*(x2-x1)*dlnm
        enddo
      enddo
      sum0 = sum0 + psum

   end subroutine integrate_2d_spt


   subroutine grid_C_spt(steps_z, steps_m, completeness)
      !tabulate completeness for SPT 1D likelihood
      real(dl), intent(in) :: steps_z(:), steps_m(:)
      integer :: nsteps_z, nsteps_m, n
      real(sp) :: completeness(:,:)
      integer :: i, ii, jj, k, kk, iostat, Nzeta
      real(dl) :: dlogxi, mp, zp, xp, c, s, D_s, xi_cut, zt_cut, sum0
      real(dl) :: fac, lnzmin, lnzmax, dlnz, zz, intsc, lnzeta, mu, zeta0, zeta1, dzeta, pzeta, arg0, arg1
      real(dl), allocatable :: tempf(:)

      D_s = cosmopar%d_spt

      xi_cut = 5.
      zt_cut = sqrt(xi_cut**2. - 3.) ! zeta_cut = 4.69
      s = 1.                         ! assuming a gaussian noise

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)

      if (D_s == 0.) then
        !no intrinsic scatter
        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = exp(steps_m(ii))
                zp = steps_z(jj)
                xp = zeta(mp,zp)

                c = erf_compl(xp,s,zt_cut)
                completeness(ii,jj) = completeness(ii,jj) + c
            enddo
        enddo

      else
        !with intrinsic scatter
        fac = 1./sqrt(2.*pi*D_s**2.)
        lnzmin = -10.   ! zeta_min = 4.53e-5
        lnzmax = 10.    ! zeta_max = 2.20e+4
        dlnz = 0.1_dl

        Nzeta = int((lnzmax - lnzmin)/dlnz)

        if (allocated(tempf)) deallocate(tempf)
        allocate(tempf(Nzeta), stat=iostat)
        if (iostat/=0) then
            print*, 'allocation error'
        endif

        tempf(:) = 0.

        lnzeta = lnzmin
        do i = 1, Nzeta
            sum0 = 0.
            zz = exp(lnzeta)
            lnzeta = lnzeta + dlnz
            sum0 =  sum0 + erf_compl(zz,s,zt_cut)
            tempf(i) = sum0
        enddo

        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = exp(steps_m(ii))
                zp = steps_z(jj)
                mu = log(zeta(mp,zp))

                intsc = 0.
                lnzeta = lnzmin
                do i = 1, Nzeta-1
                    zeta0 = exp(lnzeta)
                    zeta1 = exp(lnzeta + dlnz)
                    dzeta = zeta1 - zeta0
                    arg0 = (lnzeta - mu)/(sqrt(2.)*D_s)
                    lnzeta = lnzeta + dlnz
                    arg1 = (lnzeta - mu)/(sqrt(2.)*D_s)
                    pzeta = (tempf(i)*fac/zeta0*exp(-arg0**2.) + tempf(i+1)*fac/zeta1*exp(-arg1**2.))*0.5_dl
                    intsc = intsc + pzeta*dzeta
                enddo
                if (intsc > 1.) intsc = 1.
                if (intsc < 0.) intsc = 0.
                completeness(ii,jj) = intsc
            enddo
        enddo

        deallocate(tempf)

      endif

   end subroutine grid_C_spt


   subroutine grid_C_2d_spt(steps_z, steps_m, completeness, logxi, index_xi)
      !tabulate completeness for SPT 2D likelihood
      real(dl), intent(in) :: steps_z(:), steps_m(:), logxi(:)
      integer :: nsteps_z, nsteps_m, n
      real(sp) :: completeness(:,:,:)
      integer :: i, ii, jj, k, kk, iostat, index_xi, Nzeta
      real(dl) :: dlogxi, mp, zp, xp, xi1, xi2, c, s, D_s, xi_cut, zt_cut
      real(dl) :: fac, lnzmin, lnzmax, dlnz, zz, intsc, lnzeta, mu, zeta0, zeta1, dzeta, pzeta, arg0, arg1
      real(dl), allocatable :: tempf2d(:,:)

      D_s = cosmopar%d_spt

      xi_cut = 5.
      zt_cut = sqrt(xi_cut**2. - 3.) ! zeta_cut = 4.69
      s = 1.                         ! assuming a gaussian noise

      nsteps_z = size(steps_z)
      nsteps_m = size(steps_m)
      dlogxi = logxi(2)-logxi(1)
      n = size(logxi)

      completeness(:,:,index_xi) = 0.

      if (D_s == 0.) then
        !no intrinsic scatter
        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = exp(steps_m(ii))
                zp = steps_z(jj)
                xp = zeta(mp,zp)

                k = index_xi
                xi1 = logxi(k) - dlogxi/2.
                xi2 = logxi(k) + dlogxi/2.
                xi1 = 10.**xi1
                xi2 = 10.**xi2

                c = erf_compl(xp,s,zt_cut)*erf_compl(xp,s,xi1)*(1. - erf_compl(xp,s,xi2))
                if (k == 1) c = erf_compl(xp,s,zt_cut)*(1. - erf_compl(xp,s,xi2))
                if (k == n) c = erf_compl(xp,s,zt_cut)*erf_compl(xp,s,xi1)
                completeness(ii,jj,k) = completeness(ii,jj,k) + c

            enddo
        enddo

      else
        !with intrinsic scatter
        fac = 1./sqrt(2.*pi*D_s**2.)
        lnzmin = -10.   ! zeta_min = 4.53e-5
        lnzmax = 10.    ! zeta_max = 2.20e+4
        dlnz = 0.1_dl

        Nzeta = int((lnzmax - lnzmin)/dlnz)

        if (allocated(tempf2d)) deallocate(tempf2d)
        allocate(tempf2d(Nzeta,n), stat=iostat)
        if (iostat/=0) then
            print*, 'allocation error'
        endif

        tempf2d(:,:) = 0.

        lnzeta = lnzmin
        do i = 1, Nzeta
            zz = exp(lnzeta)
            lnzeta = lnzeta + dlnz

            k = index_xi
            xi1 = logxi(k) - dlogxi/2.
            xi2 = logxi(k) + dlogxi/2.
            xi1 = 10.**xi1
            xi2 = 10.**xi2

            c = erf_compl(zz,s,zt_cut)*erf_compl(zz,s,xi1)*(1. - erf_compl(zz,s,xi2))
            if (k == 1) c = erf_compl(zz,s,zt_cut)*(1. - erf_compl(zz,s,xi2))
            if (k == n) c = erf_compl(zz,s,zt_cut)*erf_compl(zz,s,xi1)
            tempf2d(i,k) = tempf2d(i,k) + c
        enddo

        do jj = 1, nsteps_z
            do ii = 1, nsteps_m
                mp = exp(steps_m(ii))
                zp = steps_z(jj)
                mu = log(zeta(mp,zp))

                kk = index_xi
                intsc = 0.
                lnzeta = lnzmin
                do i = 1, Nzeta-1
                    zeta0 = exp(lnzeta)
                    zeta1 = exp(lnzeta + dlnz)
                    dzeta = zeta1 - zeta0
                    arg0 = (lnzeta - mu)/(sqrt(2.)*D_s)
                    lnzeta = lnzeta + dlnz
                    arg1 = (lnzeta - mu)/(sqrt(2.)*D_s)
                    pzeta = (tempf2d(i,kk)*fac/zeta0*exp(-arg0**2.) + tempf2d(i+1,kk)*fac/zeta1*exp(-arg1**2.))*0.5_dl
                    intsc = intsc + pzeta*dzeta
                enddo
                if (intsc > 1.) intsc = 1.
                if (intsc < 0.) intsc = 0.
                completeness(ii,jj,kk) = intsc
            enddo
        enddo

        deallocate(tempf2d)

      endif

   end subroutine grid_C_2d_spt


   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!


   SUBROUTINE deltaN_yz(Z,Nz,LOGY,Ny,DN,DN2D,skyfracs,thetas,ylims,switch,qa_ytot,erf_list)
      !subroutine computing counts in y and z
      REAL(dl),INTENT(IN) :: z(:),logy(:),thetas(:),skyfracs(:),ylims(:,:),qa_ytot(:,:),erf_list(:,:)
      INTEGER, INTENT(IN) :: Ny,Nz,switch
      REAl(dl),allocatable :: grid(:,:),steps_z(:),steps_z2(:),steps_m(:)
      REAl(dl),allocatable :: ytheta5_data(:,:),ytheta5_model(:),Mz_data(:,:)
      REAl(sp),allocatable :: completeness(:,:),completeness_2d(:,:,:)
      REAl(dl),allocatable :: dif1(:),dif2(:,:),mlims(:,:)
      REAL(dl) :: DN(:),DN2D(:,:)
      REAL(dl) :: rombint,DY,sigmaM,bias,dlogy,dz,z1,z2,z_ii,sum,yi
      REAL(dl) :: zi,dzi,lnm,lnmmin,lnmmax,binlogy,binz,window,zj,dzj,biny
      REAL(dl) :: ylim,lnylim,xi,xi1,sq_sigmaM,lny1,lny2,y1,y2,sum2,frac,col1,col2,zmax
      REAL(dl) :: z_min,z_max,m_zmin,m_zmax,zmin2,zmax2,m_zmin2,m_zmax2,m_max,m_min
      REAL(dl) :: test,obs,dif_t,max_z,min_z
      REAL(dl) :: zp(1)
      INTEGER :: N,iostat,nsum,nrows,reason,Nm
      INTEGER :: I,J,K,ii,jj,nsteps_z,nsteps_m,nsteps_z1
      real(dl), PARAMETER :: dlnm = 0.05_dl
      character (LEN=200) :: filename
      integer :: iunit,ylim_choice
      integer :: k_ini,k_fin,nthetas,ntab,nsteps_z2,p1(1),p2(1),ini,npatches,p(1,1),pz,pm
      INTEGER, DIMENSION(8,2) :: values_time
      REAL(SP) :: clock_time

      ntab=size(ylims)
      nthetas=size(thetas)
      npatches=size(skyfracs)

      if (npatches==0) npatches=1 !constant ylim case

      lnmmin=31.
      lnmmax=37.
      binz=z(2)-z(1)
      biny=logy(2)-logy(1)

      !size of grid for mass
      nsteps_m=(lnmmax-lnmmin)/dlnm

      !size of grid for z
      !zmax=Z(Nz)+0.5_dl*binz
      min_z=Z(1)-0.5_dl*binz
      max_z=Z(Nz)+0.5_dl*binz

      zi=min_z
      if (zi <0) zi=0.

      nsteps_z=0
      do while (zi <= max_z)
         zi=next_z(zi,binz)
         nsteps_z=nsteps_z+1
      enddo

      allocate(steps_m(nsteps_m),steps_z(nsteps_z),stat=iostat)
      if (iostat/=0) then
         print*,'allocation error'
      endif

      !grid for mass
      lnm=lnmmin
      do i=1,nsteps_m
         steps_m(i)=lnm+dlnm/2.
         lnm=lnm+dlnm
      enddo

      !grid for z
      zi=min_z
      if (zi <0) zi=0.

      do i=1,nsteps_z
         steps_z(i)=zi
         zi=next_z(zi,binz)
      enddo

      if (steps_z(1) ==0) steps_z(1)=1.e-5

      ! from first bin in z to Nz
      SELECT CASE(SWITCH)

      CASE(1)

         allocate(grid(nsteps_m,nsteps_z),completeness(nsteps_m,nsteps_z),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
            stop
         endif
         completeness(:,:)=0.

         call grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,qa_ytot,erf_list) !tabulate completeness for m and z

         grid(:,:)=0.
         call get_grid(nsteps_z,nsteps_m,steps_z,steps_m,grid) !tabulate number counts for m and z
         DN(:)=0.
         DO I=1,Nz
            !limits of bin in z
            z1=Z(I)-0.5_dl*binz
            z2=Z(I)+0.5_dl*binz

            call integrate_m_z(grid,skyfracs,completeness,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)
            ! integrate on mass accounting for completeness
            DN(I)=sum
         ENDDO

         deallocate(grid,steps_m,steps_z,completeness,stat=iostat)
         if (iostat/=0) then
            print*,'deallocation error'
            stop
         endif

      CASE(2)

         allocate(grid(nsteps_m,nsteps_z),completeness_2d(nsteps_m,nsteps_z,ny+1),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
            stop
         endif
         completeness_2d(:,:,:)=0.

         !$OMP PARALLEL
         !$OMP do
         do J=1,Ny+1
            call grid_C_2d(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness_2d,logy,qa_ytot,erf_list,J)
            !tabulate completeness for m and z and sky patch
         ENDDO
         !$OMP end PARALLEL

         grid(:,:)=0.
         call get_grid(nsteps_z,nsteps_m,steps_z,steps_m,grid) !tabulate number counts for m and z
         DN(:)=0.

         !$OMP PARALLEL
         !$OMP do
         do J=1,Ny+1
            DO I=1,Nz+1
               !limits of bin in z
               call integrate_m_zq(grid,skyfracs,completeness_2d,steps_z,steps_m,nsteps_z,&
               nsteps_m,Z(I)-0.5_dl*binz,Z(I)+0.5_dl*binz,binz,J,dlnm,DN2D(I,J))
               ! integrate on mass and z within bins accounting for completeness
            enddo
         ENDDO
         !$OMP end PARALLEL

         deallocate(grid,steps_m,steps_z,completeness_2d,stat=iostat)
         if (iostat/=0) then
            print*,'deallocation error'
            stop
         endif

      END SELECT

      RETURN
   end SUBROUTINE DeltaN_yz

   !---------------------------------------------------------------------------!

   SUBROUTINE integrate_m_z(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,dlnm,sum)
      ! integration for 1D likelihood
      REAl(dl),intent(in) :: grid(:,:),skyfracs(:)
      REAL(sp),intent(in) :: compl(:,:)
      REAl(dl),intent(in) :: steps_z(:),steps_m(:),z1,z2,binz,dlnm
      REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3
      REAl(dl):: dnnew,dnold,lny1_new,lnylim,angsize,counts
      REAL(dl):: dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
      INTEGER :: nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches
      REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
      INTEGER :: k,kk,JJ,II,j1,j2,i,index
      INTEGER :: p(1)

      sum=0._dl
      npatches=size(skyfracs)

      ! find range i1->i2 for integration in z
      !allocate(test(nsteps_z))
      test=abs(steps_z-z1)
      p=minloc(test)
      j1=p(1)
      test=abs(steps_z-z2)
      p=minloc(test)
      j2=p(1)

      sum2=0.
      do jj=j1,j2-1
         DO ii=1,nsteps_m
            x1=steps_z(jj)
            x2=steps_z(jj+1)
            f1=grid(ii,jj)
            f2=grid(ii,jj+1)
            c1=compl(ii,jj)
            c2=compl(ii,jj+1)
            sum2=sum2+0.5*(f1*c1+f2*c2)*(x2-x1)*dlnm
         enddo
      enddo
      sum=sum+sum2

   END SUBROUTINE integrate_m_z

   !---------------------------------------------------------------------------!

   SUBROUTINE integrate_m_zq(grid,skyfracs,compl,steps_z,steps_m,nsteps_z,nsteps_m,z1,z2,binz,biny,dlnm,sum)
      ! integration for 2D likelihood
      REAl(dl),intent(in) :: grid(:,:),skyfracs(:)
      REAL(sp),intent(in) :: compl(:,:,:)
      REAl(dl),intent(in) :: steps_z(:),steps_m(:),z1,z2,binz,dlnm
      REAl(dl):: sum,xi,xi1,zi,sigmaM,sq_sigmaM,lnm,window,bias,sum2,dz,sum3
      REAl(dl):: dnnew,dnold,lny1_new,lnylim,angsize,counts
      REAL(dl):: dif,dif_t,frac,ylim,frac_same,limit_m,test_new,h,C
      INTEGER:: nsteps_z,nsteps_m,nsum,nthetas,col,k_ini,k_fin,t,k_ininew,col_old,npatches,biny
      REAl(dl):: test(nsteps_z),c1,c2,f1,f2,x1,x2
      INTEGER:: k,kk,JJ,II,j1,j2,i,index
      INTEGER:: p(1)

      !compl(nteps_m,nsteps_z,npatches)
      sum=0._dl
      npatches=size(skyfracs)

      ! find range i1->i2 for integration in z
      !allocate(test(nsteps_z))
      test=abs(steps_z-z1)
      p=minloc(test)
      j1=p(1)
      test=abs(steps_z-z2)
      p=minloc(test)
      j2=p(1)

      sum2=0.
      do jj=j1,j2-1
         DO ii=1,nsteps_m
            x1=steps_z(jj)
            x2=steps_z(jj+1)
            f1=grid(ii,jj)
            f2=grid(ii,jj+1)
            c1=compl(ii,jj,biny)
            c2=compl(ii,jj+1,biny)
            sum2=sum2+0.5*(f1*c1+f2*c2)*(x2-x1)*dlnm
         enddo
      enddo
      sum=sum+sum2

   END SUBROUTINE integrate_m_zq

   !---------------------------------------------------------------------------!

   SUBROUTINE grid_C_2d(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,logy,qa_ytot,erf_list,iy)
      !tabulate completeness for 2D likelihood
      REAl(dl),intent(in) :: steps_z(:),steps_m(:),thetas(:),ylims(:,:),skyfracs(:),logy(:),qa_ytot(:,:),erf_list(:,:)
      INTEGER :: nsteps_z,nsteps_m,nthetas,ntab,npatches,iostat,switch_comp,nsteps_y,nsteps_t,nerf,iy
      real(sp):: completeness(:,:,:) !different completeness for different bin in q
      integer :: i,j,ii,jj,index1,index2,count,P(1),k1,k2,l1,l2,k,N,nthetas2,nrows,indt(1),indy(1),it,i1y,i2y,i3y,kk
      integer :: indminy(2),indmaxy(2),iminy,imaxy,nq
      real (dl):: dif_old,dif,max,min,dlm,binz,m_min,m_max,mp,yp,zp,thp,xk1,xk2,xk3,yk1,yk2,yk3,fact,qmin,qmax,dlogy
      real(dl),allocatable:: dif_y(:),dif_theta(:),difft(:),diffy(:)
      real(dl):: min_thetas,max_thetas,min_y,max_y
      real(dl):: c1,c2,th1,th2,c,y12,y1,y2,y,col1,col0,csum,ytiny
      real(dl),allocatable:: thetas2(:),ylims2(:,:),erfs(:,:,:)
      real(dl):: win0,win,arg,arg0,y0,py,lnymax,lnymin,lny,dy,fac,mu,int,dlny,fsky
      real(dl):: sigmaM

      sigmaM=cosmopar%sigmaM
      switch_comp = 1 !ERF  (2D not implemented with QA selection function)
      ntab=size(ylims)
      nsteps_z=size(steps_z)
      nsteps_m=size(steps_m)
      nthetas=size(thetas)
      dlogy=logy(2)-logy(1)
      nq=size(logy)
      allocate(dif_y(ntab),dif_theta(nthetas))
      completeness(:,:,iy)=0.

      min_thetas=minval(thetas)
      max_thetas=maxval(thetas)

      min_y=minval(qa_ytot)
      max_y=maxval(qa_ytot)

      nerf =size(qa_ytot,1)

      if (sigmaM==0) then

         do jj=1,nsteps_z
            do ii=1,nsteps_m
               mp=exp(steps_m(ii))
               zp=steps_z(jj)
               thp=theta500(mp,zp)
               yp=y500(mp,zp)

               if (thp > max_thetas) then
                  l1=nthetas
                  l2=nthetas-1
                  th1=thetas(l1)
                  th2=thetas(l2)

               else if  (thp < min_thetas) then
                  l1=1
                  l2=2
                  th1=thetas(l1)
                  th2=thetas(l2)
               else
                  dif_theta=abs(thetas-thp)
                  P=minloc(dif_theta)
                  l1=P(1)
                  th1=thetas(l1)
                  l2=l1+1
                  if (th1 > thp) l2=l1-1
                  th2=thetas(l2)

               endif

               do i=1,npatches
                  y1=ylims(i,l1)
                  y2=ylims(i,l2)
                  y=y1+(y2-y1)/(th2-th1)*(thp-th1)!sigma at the relevant scale for the patch

                  k=iy
                  qmin=logy(k)-dlogy/2.
                  qmax=logy(k)+dlogy/2.
                  qmin=10.**qmin
                  qmax=10.**qmax
                  c2=erf_compl(yp,y,q)*erf_compl(yp,y,qmin)*(1.-erf_compl(yp,y,qmax))
                  if (k==1) c2=erf_compl(yp,y,q)*(1.-erf_compl(yp,y,qmax))
                  if (k==nq) c2=erf_compl(yp,y,q)*erf_compl(yp,y,qmin)
                  completeness(ii,jj,k)=completeness(ii,jj,k)+c2*skyfracs(i)
               enddo
            enddo
         enddo

      else
         fac=1./sqrt(2.*pi*sigmaM**2)
         lnymin=-11.5
         lnymax=10.
         dlny=0.05

         N=(lnymax-lnymin)/dlny

         fsky=0
         do i=1,npatches
            fsky=fsky+skyfracs(i)
         enddo

         allocate(erfs(N,nthetas,nq),stat=iostat)!y,integrated completeness
         if (iostat/=0) then
            print*,'allocation error'
         endif

         erfs(:,:,:)=0.
         do j=1,nthetas
            lny=lnymin
            do jj=1,N
               y0=dexp(lny)
               lny=lny+dlny

               do i=1,npatches
                  y1=ylims(i,j)

                  k=iy
                  qmin=logy(k)-dlogy/2.
                  qmax=logy(k)+dlogy/2.
                  qmin=10.**qmin
                  qmax=10.**qmax
                  c2=erf_compl(y0,y1,q)*erf_compl(y0,y1,qmin)*(1.-erf_compl(y0,y1,qmax))
                  if (k==1)  c2=erf_compl(y0,y1,q)*(1.-erf_compl(y0,y1,qmax))
                  if (k==nq) c2=erf_compl(y0,y1,qmin)*erf_compl(y0,y1,q)

                  erfs(jj,j,k)=erfs(jj,j,k)+c2*skyfracs(i)

               enddo
            enddo
         enddo

         do jj=1,nsteps_z
            do ii=1,nsteps_m

               mp=exp(steps_m(ii))
               zp=steps_z(jj)
               thp=theta500(mp,zp)
               yp=y500(mp,zp)
               if (thp > max_thetas) then
                  l1=nthetas
                  l2=nthetas-1
                  th1=thetas(l1)
                  th2=thetas(l2)
               else if  (thp < min_thetas) then
                  l1=1
                  l2=2
                  th1=thetas(l1)
                  th2=thetas(l2)
               else
                  dif_theta=abs(thetas-thp)
                  P=minloc(dif_theta)
                  l1=P(1)
                  th1=thetas(l1)
                  l2=l1+1
                  if (th1 > thp) l2=l1-1
                  th2=thetas(l2)
               endif
               y=y500(mp,zp)
               mu=dlog(y500(mp,zp))

               kk=iy
               int=0.
               lny=lnymin
               do k=1,N-1
                  y0=dexp(lny)
                  y=dexp(lny+dlny)
                  dy=y-y0
                  arg0=((lny-mu)/(sqrt(2.)*sigmaM))
                  win0=erfs(k,l1,kk)+(erfs(k,l2,kk)-erfs(k,l1,kk))/(th2-th1)*(thp-th1)
                  win=erfs(k+1,l1,kk)+(erfs(k+1,l2,kk)-erfs(k+1,l1,kk))/(th2-th1)*(thp-th1)
                  lny=lny+dlny
                  arg=((lny-mu)/(sqrt(2.)*sigmaM))
                  py=(win0*fac/y0*exp(-arg0**2)+win*fac/y*exp(-arg**2))*0.5
                  int=int+py*dy
               enddo
               if (int > fsky) int=fsky
               if (int < 0.) int=0.
               completeness(ii,jj,kk)=int
            enddo
         enddo

         deallocate(erfs)

      endif

   END SUBROUTINE grid_C_2d

   !---------------------------------------------------------------------------!

   SUBROUTINE grid_C(npatches,steps_z,steps_m,thetas,ylims,skyfracs,completeness,qa_ytot,erf_list)
      !tabulate completeness for 1D likelihood
      REAl(dl),intent(in) :: steps_z(:),steps_m(:),thetas(:),ylims(:,:),skyfracs(:),qa_ytot(:,:),erf_list(:,:)
      INTEGER :: nsteps_z,nsteps_m,nthetas,ntab,npatches,iostat,switch_comp,nsteps_y,nsteps_t,nerf
      real(sp):: completeness(:,:)
      integer :: i,j,ii,jj,index1,index2,count,P(1),k1,k2,l1,l2,k,N,nthetas2,nrows,indt(1),indy(1),it,i1y,i2y,i3y
      integer :: indminy(2),indmaxy(2),iminy,imaxy
      real (dl):: dif_old,dif,max,min,dlm,binz,m_min,m_max,mp,yp,zp,thp,xk1,xk2,xk3,yk1,yk2,yk3,fact
      real(dl),allocatable:: dif_y(:),dif_theta(:),difft(:),diffy(:)
      real(dl):: min_thetas,max_thetas,min_y,max_y
      real(dl):: c1,c2,th1,th2,c,y12,y1,y2,y,col1,col0
      real(dl),allocatable   :: thetas2(:),ylims2(:,:),erfs(:,:)
      real(dl):: win0,win,arg,arg0,y0,py,lnymax,lnymin,lny,dy,fac,mu,int,dlny,fsky
      real(dl):: sigmaM

      sigmaM=cosmopar%sigmaM
      !choice to use error function completeness or QA
      switch_comp = 1 !ERF

      !switch_comp = 2 !QA
      !end choice to use error function completeness or QA

      ntab=size(ylims)
      nsteps_z=size(steps_z)
      nsteps_m=size(steps_m)
      nthetas=size(thetas)

      allocate(dif_y(ntab),dif_theta(nthetas))

      min_thetas=minval(thetas)
      max_thetas=maxval(thetas)

      min_y=minval(qa_ytot)
      max_y=maxval(qa_ytot)

      nerf = size(qa_ytot,1)

      SELECT CASE(switch_comp)

      CASE(1)
         !error function
         if (sigmaM==0) then
            ! no scatter in y-m relation
            do jj=1,nsteps_z
               do ii=1,nsteps_m
                  mp=exp(steps_m(ii))
                  zp=steps_z(jj)
                  thp=theta500(mp,zp)
                  yp=y500(mp,zp)
                  if (thp > max_thetas) then
                     l1=nthetas
                     l2=nthetas-1
                     th1=thetas(l1) ! last one = largest theta = thp
                     th2=thetas(l2) ! before that

                  else if  (thp < min_thetas) then
                     l1=1
                     l2=2
                     th1=thetas(l1) ! first one = smallest theta  = thp
                     th2=thetas(l2) ! next
                  else
                     dif_theta=abs(thetas-thp)
                     P=minloc(dif_theta)
                     l1=P(1)
                     th1=thetas(l1)  ! the closest one = thp
                     l2=l1+1
                     if (th1 > thp) l2=l1-1
                     th2=thetas(l2)
                     ! consequently, set th1 = closest to thp
                     !                   th2 = slightly larger or smaller than th1
                  endif
                  completeness(ii,jj)=0.

                  do i=1,npatches
                     y1=ylims(i,l1)
                     y2=ylims(i,l2)
                     y=y1+(y2-y1)/(th2-th1)*(thp-th1)!sigma at the relevant scale for the patch
                     c2=erf_compl(yp,y,q)
                     completeness(ii,jj)=completeness(ii,jj)+c2*skyfracs(i)
                  enddo
               enddo
            enddo

         else
            ! scatter in y-m relation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            fac=1./sqrt(2.*pi*sigmaM**2)
            lnymin=-11.5
            lnymax=10.
            dlny=0.05

            N=(lnymax-lnymin)/dlny

            fsky=0
            do i=1,npatches
               fsky=fsky+skyfracs(i)
            enddo

            allocate(erfs(N,nthetas),stat=iostat)!y,integrated completeness
            if (iostat/=0) then
               print*,'allocation error'
            endif

            do j=1,nthetas
               lny=lnymin
               do k=1,N
                  y0=dexp(lny)
                  lny=lny+dlny
                  win0=0.
                  do i=1,npatches
                     y1=ylims(i,j) ! sigma
                     win0=win0+erf_compl(y0,y1,q)*skyfracs(i)
                  enddo
                  erfs(k,j)=win0
                  !           if (erfs(k,j)>fsky) erfs(k,j)=fsky
               enddo
            enddo

            do jj=1,nsteps_z
               do ii=1,nsteps_m

                  mp=exp(steps_m(ii))
                  zp=steps_z(jj)
                  thp=theta500(mp,zp)
                  yp=y500(mp,zp)

                  if (thp > max_thetas) then
                     l1=nthetas
                     l2=nthetas-1
                     th1=thetas(l1)
                     th2=thetas(l2)
                  else if  (thp < min_thetas) then
                     l1=1
                     l2=2
                     th1=thetas(l1)
                     th2=thetas(l2)
                  else
                     dif_theta=abs(thetas-thp)
                     P=minloc(dif_theta)
                     l1=P(1)
                     th1=thetas(l1)
                     l2=l1+1
                     if (th1 > thp) l2=l1-1
                     th2=thetas(l2)
                  endif
                  y=y500(mp,zp)
                  mu=dlog(y500(mp,zp))
                  !print*, y500(mp,zp), mp, zp

                  int=0.
                  lny=lnymin
                  do k=1,N-1
                     y0=dexp(lny)
                     y=dexp(lny+dlny)
                     dy=y-y0
                     arg0=((lny-mu)/(sqrt(2.)*sigmaM))
                     win0=erfs(k,l1)+(erfs(k,l2)-erfs(k,l1))/(th2-th1)*(thp-th1)
                     win=erfs(k+1,l1)+(erfs(k+1,l2)-erfs(k+1,l1))/(th2-th1)*(thp-th1)
                     lny=lny+dlny
                     arg=((lny-mu)/(sqrt(2.)*sigmaM))
                     py=(win0*fac/y0*exp(-arg0**2)+win*fac/y*exp(-arg**2))*0.5
                     int=int+py*dy
                     !print*, int
                  enddo

                  if (int > fsky) int=fsky
                  if (int < 0.) int=0.
                  completeness(ii,jj)=int
                  !print*, int, fsky

               enddo
            enddo

            deallocate(erfs)

         endif

      CASE(2)
         !QA COMPLETENESS
         if (sigmaM==0) then
            ! no scatter in y-m relation
            do jj=1,nsteps_z
               do ii=1,nsteps_m

                  mp=exp(steps_m(ii))
                  zp=steps_z(jj)

                  thp=theta500(mp,zp)
                  yp=y500(mp,zp)

                  allocate(difft(nthetas))
                  difft = abs(thetas-thp)
                  indt = minloc(difft)
                  it = indt(1)

                  allocate(diffy(nerf))
                  diffy = abs(qa_ytot(:,it)-yp)

                  indy = minloc(diffy)

                  i1y = indy(1)-1
                  i2y = indy(1)+1

                  deallocate(difft)
                  deallocate(diffy)

                  indminy = minloc(qa_ytot)
                  indmaxy = maxloc(qa_ytot)

                  if (i1y <=0) then
                     i1y = indminy(1)
                  endif

                  if (i2y>nerf) then
                     i2y = indmaxy(1)
                  endif

                  if (yp < min_y) then
                     i1y = indminy(1)
                     i2y = indminy(1)+1
                  endif

                  if (yp > max_y) then
                     i1y = indmaxy(1)-1
                     i2y = indmaxy(1)
                  endif

                  xk1 = qa_ytot(i1y,it)
                  xk2 = qa_ytot(i2y,it)
                  yk1 = erf_list(i1y,it)
                  yk2 = erf_list(i2y,it)

                  completeness(ii,jj) = 0.

                  !INTERPOLATION POLY ORDER 1
                  completeness(ii,jj) = exp(log(yk1)+((log(yp)-log(xk1))/(log(xk2)-log(xk1)))*(log(yk2)-log(yk1)))

                  if (completeness(ii,jj)> 1.) completeness(ii,jj)=1.
                  if (completeness(ii,jj)==0.) completeness(ii,jj)=1e-20

                  fact = sum(skyfracs)
                  completeness(ii,jj) = completeness(ii,jj)*fact

               enddo
            enddo
         else
            ! scatter in y-m relation
            print*, 'QA completeness and scatter in Y-M relation'
            print*, 'CASE NOT IMPLEMENTED - STOPPING'
            stop
         endif

      END SELECT

   END SUBROUTINE grid_C

   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!
   !---------------------------------------------------------------------------!

   SUBROUTINE get_grid(nz,nm,z,lnm,grid)
      !tabulate counts from theory
      REAL(dl):: z(:),lnm(:),grid(:,:)
      INTEGER :: I,J,nz,nm
      REAL(dl):: Mnew,g,ynew,vol,theta

      DO I=1,Nz
         g = delta(z(I))
         vol = dVdzdO(z(I))
         DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2*vol
         ENDDO
      ENDDO
   END SUBROUTINE get_grid

   SUBROUTINE get_grid_act(nz,nm,z,lnm,grid)
      !tabulate counts from theory
      REAL(dl):: z(:),lnm(:),grid(:,:)
      INTEGER :: I,J,nz,nm
      REAL(dl):: Mnew,g,ynew,vol,theta

      DO I=1,Nz
         g = delta(z(I))
         vol = dVdzdO(z(I))
         DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2_a*vol
         ENDDO
      ENDDO
   END SUBROUTINE get_grid_act

   SUBROUTINE get_grid_spt(nz,nm,z,lnm,grid)
      !tabulate counts from theory
      REAL(dl):: z(:),lnm(:),grid(:,:)
      INTEGER :: I,J,nz,nm
      REAL(dl):: Mnew,g,ynew,vol,theta

      DO I=1,Nz
         g = delta(z(I))
         vol = dVdzdO(z(I))
         DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2_s*vol
         ENDDO
      ENDDO
   END SUBROUTINE get_grid_spt

   SUBROUTINE get_grid_SO(nz,nm,z,lnm,grid)
      !tabulate counts from theory
      REAL(dl):: z(:),lnm(:),grid(:,:)
      INTEGER :: I,J,nz,nm
      REAL(dl):: Mnew,g,ynew,vol,theta

      DO I=1,Nz
         g = delta(z(I))
         vol = dVdzdO(z(I))
         DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2_SOfull*vol
         ENDDO
      ENDDO
   END SUBROUTINE get_grid_SO

   SUBROUTINE get_grid_SO_test(nz,nm,z,lnm,grid)
      !tabulate counts from theory
      REAL(dl):: z(:),lnm(:),grid(:,:)
      INTEGER :: I,J,nz,nm
      REAL(dl):: Mnew,g,ynew,vol,theta

      DO I=1,Nz
         g = delta(z(I))
         vol = dVdzdO(z(I))
         DO J=1,Nm
            Mnew=exp(lnm(J))
            grid(j,i) = dndlnM(z(I),Mnew,g)*surveypar%deg2_SOtest*vol
         ENDDO
      ENDDO
   END SUBROUTINE get_grid_SO_test

   !---------------------------------------------------------------------------!

end module numbercounts

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

module szcounts
   use settings
   use CosmologyTypes
   use CosmoTheory
   use Calculator_Cosmology
   use Likelihood_Cosmology
   use IniObjects
   use massfunction
   use massobservable
   use numbercounts !:D

   implicit none
   private
   logical :: do_sz_init = .true.
   logical :: use_sz = .false.
   integer :: SZ_num
   integer :: Ny, Nz, Nxi
   integer :: Nz_a, Nz_s, Ny0 !!!!! for combined case
   integer :: dim_switch_plc, dim_switch_act, dim_switch_spt  ! :D !!!!! 1 for 1D, 2 for 2D
   integer :: ylim_switch_act !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   integer :: nmiss_switch  !0 for simple rescaling, 1 for MCMC
   integer :: errorz_switch  !0 for simple rescaling, 1 for MCMC
   integer :: print_counts !to print theory counts (for action=4)
   integer :: nmiss,nred2,ncat,nerf,nrows_qa,nrows_erf
   integer :: nmiss_a, nmiss_s, nred2_a, nred2_s, ncat_a, ncat_s !!!!! for combined case
   integer :: ylim_switch  !>1 for constant ylim, <1 for variable ylim
   real(dl), allocatable :: DNcat(:,:),DNzcat(:),DNycat(:),DN(:,:),DNz(:),DNy(:)
   real(dl), allocatable :: delN(:), delN2D(:,:), y0(:), delNcat(:), delN2Dcat(:,:) !!!!! for ACT
   real(dl), allocatable :: delNzcat(:), delNxcat(:), delNz(:), delNx(:), logxi(:) !!!!! for SPT
   real(dl), allocatable :: Z(:),LOGY(:),ylims(:,:),thetas(:),skyfracs(:),erf_list(:,:),qa_ytot(:,:)
   real(dl), allocatable :: y0lim_arr(:), noise_SOtest(:), noise_SOfull(:), noise_act_arr(:)
   integer, allocatable :: tilenames(:)
   real(dl), allocatable :: SZcat(:,:)
   real(dl), allocatable :: SZcat_a(:,:), z_a(:), SZcat_s(:,:), z_s(:) !!! for combined case
   real(dl), allocatable :: y0lims(:), skyfracs_SOtest(:), skyfracs_SOfull(:), theta(:), Q00(:,:), Q0(:)
   real(dl) :: clash,wtg,lens,cccp,pns,oh2,palpha,pystar,psigma,pbeta ! switches for priors =0 no prior, =1 prior
   real(dl) :: pA_s, pB_s, pC_s, pD_s, pH0
   real(dl) :: snr_act_cat, ylim_act_cat, snr_spt_cat, snr_so_cat !!!!! for ACT & SPT catalogue cut
   real(dl) :: y0binN !!!!! number of y0 bins for 2D ACT
   real :: sz_kmax = 4.0
   character(len=256) :: SZ_filename = ''
   character(LEN=*), parameter :: SZ_version = 'Nov_2014'
   REAL(dl), SAVE :: zmaxs,z0,dz
   real(dl), save :: zmaxs_a, zmaxs_s !!!!! for combined case
   INTEGER, SAVE :: Nred, Nscat

   type, extends(TCosmoCalcLikelihood) :: SZLikelihood
   contains
   procedure :: LogLike => SZCC_Cash !!!!!!!!! so it was loglike
   end type SZLikelihood

   PUBLIC :: SZLikelihood_Add, SZcc_Cash

   contains

   !---------------------------------------------------------------------------!

   subroutine SZLikelihood_Add(LikeList, Ini)
      ! interface of the module with cosmomcplanck code
      ! user-defined settings are interfaced here
      implicit none
      class(TLikelihoodList) :: LikeList
      class(TSettingIni) :: ini
      integer count
      Type(SZLikelihood), pointer :: this

      print_counts=0

      if (Ini%Read_Int('action')==4) print_counts=1
      if (Ini%Read_Logical('use_SZ',.false.)) then
      allocate(this)
      this%needs_background_functions = .true.
      this%needs_powerspectra = .true.
      this%kmax = sz_kmax
      this%needs_sigmaR = .true.
      this%version = SZ_version
      call this%loadParamNames(trim(DataDir)//'SZ.paramnames')
      call LikeList%Add(this)
      this%LikelihoodType = 'SZ'
      this%name='SZ'

      clash=0.  !swithes for the priors (get multiplied to priors in SZCC_Cash)
      wtg=0.    !either 0 or 1
      lens=0.
      cccp=0.
      oh2=0.
      pns=0.
      palpha=0.
      pystar=0.
      psigma=0.
      pbeta=0.
      pA_s=0.
      pB_s=0.
      pC_s=0.
      pD_s=0.
      pH0=0.

      if (Ini%Read_Logical('prior_alpha_SZ',.false.)) then
         palpha=1.
         print*,'prior on alpha_SZ'
      endif

      if (Ini%Read_Logical('prior_ystar_SZ',.false.)) then
         pystar=1.
         print*,'prior on ystar_SZ'
      endif

      if (Ini%Read_Logical('prior_scatter_SZ',.false.)) then
         psigma=1.
         print*,'prior on scatter_SZ'
      endif

      if (Ini%Read_Logical('prior_beta_SZ',.false.)) then
         pbeta=1.
         print*,'prior on beta_SZ'
      endif

      if (Ini%Read_Logical('prior_A_SZ_SPT',.false.)) then
         pA_s=1.
         print*,'prior on A_SZ for SPT'
      endif

      if (Ini%Read_Logical('prior_B_SZ_SPT',.false.)) then
         pB_s=1.
         print*,'prior on B_SZ for SPT'
      endif

      if (Ini%Read_Logical('prior_C_SZ_SPT',.false.)) then
         pC_s=1.
         print*,'prior on C_SZ for SPT'
      endif

      if (Ini%Read_Logical('prior_D_SZ_SPT',.false.)) then
         pD_s=1.
         print*,'prior on D_SZ for SPT'
      endif

      if (Ini%Read_Logical('prior_H0',.false.)) then
         pH0=1.
         print*,'prior on H0'
      endif

      !!$        if (Ini%Read_Logical('prior_clash',.false.)) then
      !!$            clash=1.
      !!$            print*,'CLASH prior on SZ mass bias'
      !!$        endif

      if (Ini%Read_Logical('prior_wtg',.false.)) then
         wtg=1.
         print*,'WtG prior on SZ mass bias'
      endif

      if (Ini%Read_Logical('prior_lens',.false.)) then
         print*,'Planck lensing prior on SZ mass bias'
         lens=1.
      endif

      if (Ini%Read_Logical('prior_cccp',.false.)) then
         print*,'CCCP prior on SZ mass bias'
         cccp=1.
      endif

      if (Ini%Read_Logical('prior_ns',.false.))then
         print*,'Prior on ns'
         pns=1.
      endif

      if (Ini%Read_Logical('prior_omegabh2',.false.)) then
         oh2=1.
         print*,'Prior on omegabh2'
      endif

      massfnpar%psind=2 !mass function = tinker mass function
      if (Ini%Read_Logical('use_watson',.false.)) then
         massfnpar%psind=3
      endif

      mlim_switch = 2
      if (Ini%Read_Logical('constant_Mlim',.false.)) then
         mlim_switch = 1
         print*,'constant Mlim'
      endif

      if (Ini%Read_Logical('Mlim_of_redshift',.false.)) then
         mlim_switch = 2
         print*,'Mlim as a function of redshift'
      endif

      ylim_switch_act = 1
      if (Ini%Read_Logical('constant_ylim_act',.false.)) then
         ylim_switch_act = 1
         print*,'constant ylim for the entire survey area'
      endif

      if (Ini%Read_Logical('ylim_act_from_map',.false.)) then
         ylim_switch_act = 2
         print*,'ylim from the map'
      endif

      survey_switch = 0
      if (Ini%Read_Logical('Planck',.false.)) then
         survey_switch = 1
         print*,'only Planck survey is chosen'
      endif

      if (Ini%Read_Logical('ACT',.false.)) then
         survey_switch = 2
         print*,'only ACT survey is chosen'
      endif

      if (Ini%Read_Logical('Planck_and_ACT',.false.)) then
         survey_switch = 3
         print*,'Both Planck and ACT surveys are chosen'
      endif

      if (Ini%Read_Logical('SPT',.false.)) then
         survey_switch = 4
         print*,'only SPT survey is chosen'
      endif

      if (Ini%Read_Logical('Planck_and_SPT',.false.)) then
         survey_switch = 5
         print*,'Both Planck and SPT surveys are chosen'
      endif

      if (Ini%Read_Logical('SO',.false.)) then
         survey_switch = 6
         print*,'only SO survey is chosen'
      endif

      if (Ini%Read_Logical('SO_tiles',.false.)) then
         survey_switch = 7
         print*,'only SO tiles survey is chosen'
      endif

      dim_switch_plc = 2 !initialization
      if ((Ini%Read_Logical('Planck',.false.) .or. Ini%Read_Logical('Planck_and_ACT',.false.) .or. Ini%Read_Logical('Planck_and_SPT',.false.)) .and. Ini%Read_Logical('Planck_1D',.false.)) then
         dim_switch_plc = 1
         print*,'Planck 1D SZ likelihood dN/dz'
      endif

      if ((Ini%Read_Logical('Planck',.false.) .or. Ini%Read_Logical('Planck_and_ACT',.false.) .or. Ini%Read_Logical('Planck_and_SPT',.false.)) .and. Ini%Read_Logical('Planck_2D',.false.)) then
         dim_switch_plc = 2
         print*,'Planck 2D SZ likelihood dN/dzdq'
      endif

      if ((Ini%Read_Logical('Planck_2D',.false.)) .and. dim_switch_plc == 1 ) then
         print*,'Error 1D and 2D likelihood both selected'
         stop
      endif

      dim_switch_act = 1 !initialization
      if ((Ini%Read_Logical('ACT',.false.) .or. Ini%Read_Logical('Planck_and_ACT',.false.)) .and. Ini%Read_Logical('ACT_1D',.false.)) then
         dim_switch_act = 1
         print*,'ACT 1D SZ likelihood dN/dz'
      endif

      if ((Ini%Read_Logical('ACT',.false.) .or. Ini%Read_Logical('Planck_and_ACT',.false.)) .and. Ini%Read_Logical('ACT_2D',.false.)) then
         dim_switch_act = 2
         y0binN = Ini%Read_Double('y0binN', 2._dl)
         y0bin = y0binN
         print*,'ACT 2D SZ likelihood dN/dy0'
      endif

      if ((Ini%Read_Logical('ACT_2D',.false.)) .and. dim_switch_act == 1 ) then
         print*,'Error 1D and 2D likelihood both selected'
         stop
      endif

      dim_switch_spt = 1 !initialization
      if ((Ini%Read_Logical('SPT',.false.) .or. Ini%Read_Logical('Planck_and_SPT',.false.)) .and. Ini%Read_Logical('SPT_1D',.false.)) then
         dim_switch_spt = 1
         print*,'SPT 1D SZ likelihood dN/dz'
      endif

      if ((Ini%Read_Logical('SPT',.false.) .or. Ini%Read_Logical('Planck_and_SPT',.false.)) .and. Ini%Read_Logical('SPT_2D',.false.)) then
         dim_switch_spt = 2
         print*,'SPT 2D SZ likelihood dN/dzdxi'
      endif

      if (Ini%Read_Logical('SPT_2D',.false.) .and. dim_switch_spt == 1) then
         print*,'Error 1D and 2D likelihood both selected'
         stop
      endif

      scatter_mlim_switch = 2
      if (Ini%Read_Logical('no_scatter_in_Mlim',.false.)) then
         scatter_mlim_switch = 1
         print*, 'no scatter in Mlim'
      endif

      if (Ini%Read_Logical('including_scatter_in_Mlim',.false.)) then
         scatter_mlim_switch = 2
         print*, 'including scatter in Mlim'
      endif

      !completeness_act_switch = 1
      if (Ini%Read_Logical('completeness_act',.false.)) then
         completeness_act_switch = 1
         print*, 'considering completeness in ACT'
      endif

      !if (dim_switch_plc == 0) then
      !    dim_switch_plc = 2 !default
      !    print*,'Planck 2D SZ likelihood dN/dzdq'
      !endif

      ! Add signal-to-noise or ylim cut from catalogue
      snr_act_cat  = Ini%Read_Double('snr_act_cat', 5._dl) ! default is 5 (SNR_{2.4} cut of ACTPol(2017) ylim map)
      snr_spt_cat  = Ini%Read_Double('snr_spt_cat', 5._dl) ! (detection significance of SPT full catalogue : 4.5)
      snr_so_cat   = Ini%Read_Double('snr_so_cat', 5._dl) ! default is 5
      ylim_act_cat = Ini%Read_Double('ylim_act_cat', 0.5_dl) ! default is 0.5 [10^-4]

      CALL SZ_init

      endif
   end subroutine SZLikelihood_Add

   !---------------------------------------------------------------------------!

   subroutine SZ_init
      !initialization of SZ module
      implicit none
      character (LEN=20):: name
      character (LEN=200):: cat_filename,skyfracs_filename,ylims_filename,thetas_filename,filename, Qfunc_filename, theta_filename
      real :: dummy
      real :: dzorg
      integer :: i,j,jj,iostat,nrows,reason,ii,iunit,nrows_old,nthetas,npatches
      integer :: nrows_a, nrows_s !!!!!  for combined case
      integer :: Nfact,error,nq
      real(dl) :: ymin,ymax,dlogy,logymax,logymin,yi,col1,col2,sum,col3,col4,col5,col6,col7,col8 !!!!!!!!!!!!!!!!!!!!!!!
      real(dl) :: y_min,y_max,z_min,z_max,factorial,qmin,qmax,dq,qbin,norm
      real(dl) :: logximin, logximax, dlogxi, xii, xi_min, xi_max
      character (Len=200) :: coltxt
      integer :: ncols

      error = 0

      select case(survey_switch)
      !1 Planck
      !2 ACT
      !3 Planck + ACT
      !4 SPT
      !5 Planck + SPT
      !6 SO single tile
      !7 SO full map

      !------------------------------------------------------------------------!

      case(1) ! survey_switch 1 = Planck

         !file names
         cat_filename='data/SZ_cat.txt'
         !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
         thetas_filename='data/SZ_thetas.txt'
         skyfracs_filename='data/SZ_skyfracs.txt'
         ylims_filename='data/SZ_ylims.txt'

         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         surveypar%deg2=41253.0  !full sky (sky fraction handled in skyfracs file)

         z0=0.0d0
         zmaxs=1.
         dz=0.1d0

         surveypar%ylimin=1.e-3
         surveypar%ymaxin=1.e-1
         logymin=0.7 !s/n=6
         logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
         dlogy=0.25

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         ylim_switch=-1  ! >0=constant ylim, <0=variable ylim
         nmiss_switch=0
         errorz_switch=0

         surveypar%deg2 = 3.046174198d-4*surveypar%deg2
         Nz = DINT((zmaxs-z0)/dz)+1
         Ny = DINT((logymax-logymin)/dlogy)+1
         print*,'Nz=',Nz
         print*,'Ny=',Ny

         !ylims file                    for the variable ylim case
         if (ylim_switch < 0) then
            ymin=-1.*ymin

            ! read catalogue
            filename=cat_filename

            nrows=0
            print*,'Reading catalogue'

            open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
            IF (Reason > 0)  THEN
               print*,'Error opening file',filename
            endif

            print*,filename
            nrows=1
            ncat=1 !number of clusters above s/n threshold

            DO
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
               !print*,Reason
               IF (Reason > 0)  THEN
                  print*,'Error in reading file'
                  print*,filename
                  stop
               ELSE IF (Reason < 0) THEN
                  exit
               ELSE
                  nrows=nrows+1
                  if (col3 >=q) ncat=ncat+1
               END IF
            END DO

            CLOSE (unit=iunit)
            print*,'done'

            Print*,'Catalogue Number of clusters=',nrows
            Print*,'Catalogue Number of clusters above the S/N threshold of',q,'=',ncat

            if (allocated(SZcat)) deallocate(SZcat)
            allocate(SZcat(ncat,3),stat=iostat)! z,y

            if (iostat/=0) then
               print*,'allocation error'
            endif

            open (newunit=iunit,file=filename,status='old',form="formatted")
            ii=1
            DO i=1,nrows
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
               if (col3 >=q) then   ! exclude clusters below s/n threshold q
                  SZcat(ii,1)=col1  ! redshift
                  SZcat(ii,2)=col2  ! error on redshift
                  SZcat(ii,3)=col3  ! detection S/N
                  ii=ii+1
               endif
            ENDDO
            CLOSE (unit=iunit)

            ! load files describing selection function

            !file with theta
            ! theta, nbins, first_index, last_index,first_index2, last_index2
            call File%LoadTxt(thetas_filename, thetas, n=nthetas)
            print*,'Number of size thetas=',nthetas

            !file with skyfracs
            ! theta, nbins, first_index, last_index,first_index2, last_index2
            call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
            print*,'Number of patches=',npatches

            filename=ylims_filename
            nrows=File%TxtFileLines(filename)
            print*,'Number of size y =',nrows
            if (nrows /= npatches*nthetas) then
               print*,'Format error for ylims.txt:'
               print*,'Expected rows:',npatches*nthetas
               print*,'Actual rows:',nrows
               stop
            endif

            if (allocated(ylims)) deallocate(ylims)
            allocate(ylims(npatches,nthetas),stat=iostat)

            if (iostat/=0) then
               print*,'allocation error'
            endif

            open (newunit=iunit,file=filename,status='unknown',form="formatted")
            i=1
            j=1
            DO ii=1,nrows
               READ(iunit,*,IOSTAT=Reason)  col1
               ylims(i,j)=col1
               i=i+1
               if (i > npatches) then
                  i=1
                  j=j+1
               endif
            ENDDO
            CLOSE (unit=iunit)
         endif

         if (allocated(z)) deallocate(z)
         if (allocated(logy)) deallocate(logy)
         allocate(z(Nz),logy(Ny+1),stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         ! logy vector
         yi=logymin+dlogy/2.
         DO I=1,Ny+1
            logy(I)=yi
            yi=yi+dlogy
         END DO

         print*,'q =',10.d0**logy

         ! z vector
         DO I=0,Nz-1
            Z(I+1)=z0+I*dz+0.5_dl*dz
         END DO
         if (z0==0._dl) Z(1)=Z(1)+1.e-8 ! for numerical problem when starting from 0.

         if (allocated(DNcat)) deallocate(DNcat)
         if (allocated(DNzcat)) deallocate(DNzcat)
         if (allocated(DNycat)) deallocate(DNycat)
         if (allocated(DN)) deallocate(DN)
         if (allocated(DNz)) deallocate(DNz)
         if (allocated(DNy)) deallocate(DNy)
         allocate(DNcat(Nz,Ny+1),DNzcat(Nz),DNycat(Ny),DN(Nz+1,Ny+1),DNz(Nz),DNy(Ny),stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         DNcat(:,:)=0.
         DNzcat(:)=0.

         SELECT CASE(dim_switch_plc)
         CASE(1) ! Planck 1D N(z)

            nrows=size(SZcat(:,1))

            !!!!! number of catalogue clusters without redshift (they are -1 in .txt file)
            nmiss=0
            DO ii=1,nrows
               if (SZcat(ii,1) <0.) nmiss=nmiss+1.
            enddo

            !!!!! number of catalogue clusters in each redshift bin
            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz
               DO ii=1,nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                     DNzcat(I)=DNzcat(I)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            !!!!! just total number of catalogue clusters
            sum=0.
            DO I=1,Nz
               sum=sum+DNzcat(I)
            ENDDO

            nred2=nrows-nmiss   !!!!! number of catalogue clusters with redshift
            print*,nrows,nmiss
            ncat=nrows          !!!!! total number of catalogue clusters

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss
            print*,'Counts:',DNzcat !!!!! array of number of catalogue clusters for 11 redshift bins

            if (nmiss==0) nmiss_switch=0 !!!!! nrows/nred2 = 1 so it's same

            SELECT CASE(nmiss_switch) !!!!! nmiss_switch=0 is already set at the beginning
            CASE(0)
               print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)   !!!!! 1.01385681293303
            CASE(1)
               print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2) !!!!! should be the case for MCMC but how?
            END SELECT

         CASE(2) ! Planck 2D N(z,q)

            ! compute catalogue counts in z and q
            ! compute P(q|qm) once and for all

            nrows=size(SZcat(:,1))

            nmiss=0
            DO ii=1,nrows
               if (SZcat(ii,1) <0.d0) nmiss=nmiss+1.
            enddo

            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz
               DO J=1,Ny
                  y_min=logY(J)-dlogy/2.
                  y_max=logY(J)+dlogy/2.
                  y_min=10.d0**y_min
                  y_max=10.d0**y_max
                  DO ii=1,nrows
                     if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                     (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                        DNcat(I,J)=DNcat(I,J)+1.
                     endif
                  ENDDO
               ENDDO
               J=Ny+1 ! the last bin contains all S/N greater than what in the previous bin
               y_min=y_max
               DO ii=1,nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,3) >= y_min)) then
                     DNcat(I,J)=DNcat(I,J)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            print*, 'Catalogue counts:' !!!!!!!!!!!
            DO j = 1, Ny+1
               sum = 0.
               DO i = 1, Nz
                  sum = sum + DNcat(I,J)
               ENDDO
               print*, j, sum ! by signal-to-noise
            END DO

            sum = 0.
            DO I = 1, Nz
               DO J = 1, Ny+1
                  sum = sum + DNcat(I,J)
               ENDDO
               !print*, i, DNcat(i,1), DNcat(i,2), DNcat(i,3), DNcat(i,4) DNcat(i,5) !!!!!!!!!!
            END DO
            print*, 'total', sum

            !missing redshifts
            DO J=1,Ny
               y_min=logY(J)-dlogy/2.
               y_max=logY(J)+dlogy/2.
               y_min=10.**y_min
               y_max=10.**y_max
               DO ii=1,nrows
                  if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                     norm=0.
                     do jj=1,Nz
                        norm=norm+DNcat(jj,J)
                     enddo
                     DNcat(:,J)=DNcat(:,J)*(norm+1.d0)/norm
                  endif
               ENDDO
            ENDDO
            J=Ny+1 ! the last bin contains all S/N greater than what in the previous bin
            y_min=y_max
            DO ii=1,nrows
               if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) >= y_min)) then
                  norm=0.
                  do jj=1,Nz
                     norm=norm+DNcat(jj,J)
                  enddo
                  DNcat(:,J)=DNcat(:,J)*(norm+1.d0)/norm
               endif
            ENDDO
            !end missing redshifts

            nred2=nrows-nmiss
            ncat=nrows

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss

         END SELECT

      !------------------------------------------------------------------------!

      case(2) ! :D survey switch 2 = ACT

         !file names
         cat_filename='data/ACT_SZ_cat.txt'                         ! :D ! (182) full catalogue with SNR > 4         z / zErr / SNR / y0 / SNR_{2.4} / y0Err / M500unc / M500upp
         ylims_filename='data/ACT_SZ_ylims.txt'                    ! :D ! (48) with SNR_{2.4} > 5 map      my map :D hahahaha
         !ylims_filename='data/actpol_RMSTab01_56556.txt'           ! :D ! (56556) with SNR_{2.4} > 5 map  ORIGINAL DATA
         !ylims_filename='data/actpol_RMSTab02_26307.txt'           ! :D ! (26307) with SNR_{2.4} > 5 map   2 x grid size
         !ylims_filename='data/actpol_RMSTab04_9836.txt'            ! :D ! (9836)  with SNR_{2.4} > 5 map   4 x grid size
         !ylims_filename='data/actpol_RMSTab08_2991.txt'            ! :D ! (2991)  with SNR_{2.4} > 5 map   8 x grid size
         !ylims_filename='data/actpol_RMSTab12_1926.txt'            ! :D ! (1926)  with SNR_{2.4} > 5 map  12 x grid size
         !ylims_filename='data/actpol_RMSTab16_943.txt'             ! :D ! (943)   with SNR_{2.4} > 5 map  16 x grid size
         !ylims_filename='data/actpol_RMSTab32_197.txt'             ! :D ! (197)   with SNR_{2.4} > 5 map  32 x grid size
         !ylims_filename='data/actpol_RMSTab48_66.txt'              ! :D ! (66)    with SNR_{2.4} > 5 map  48 x grid size
         !ylims_filename='data/actpol_RMSTab_downsample_1.txt'       ! :D ! (3305)  with SNR_{2.4} > 5 map  step size 0.0001e-4
         !ylims_filename='data/actpol_RMSTab_downsample_5.txt'      ! :D ! (789)   with SNR_{2.4} > 5 map  step size 0.0005e-4
         !ylims_filename='data/actpol_RMSTab_downsample_10.txt'     ! :D ! (418)   with SNR_{2.4} > 5 map  step size 0.001e-4
         !ylims_filename='data/actpol_RMSTab_downsample_50.txt'     ! :D ! (103)   with SNR_{2.4} > 5 map  step size 0.005e-4
         !ylims_filename='data/actpol_RMSTab_downsample_100.txt'    ! :D ! (55)    with SNR_{2.4} > 5 map  step size 0.01e-4
         !ylims_filename='data/actpol_RMSTab_downsample_500.txt'    ! :D ! (12)    with SNR_{2.4} > 5 map  step size 0.05e-4


         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         surveypar%deg2_a = 987.5 ! ACT survey area (E-D56 field) in unit of deg2

         z0 = 0.0d0
         zmaxs = 1.4    ! :D !
         dz = 0.1d0

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         nmiss_switch=0
         errorz_switch=0

         surveypar%deg2_a = 3.046174198d-4*surveypar%deg2_a
         ! from deg2 to steradian (unit of solid angle)
         ! 3.046174198d-4 = (pi/180 deg)^2

         Nz = dint((zmaxs-z0)/dz)+1
         Ny0 = y0binN

         print*, 'Nz=', Nz
         !       print*, 'Ny0=', Ny0

         ! read catalogue
         filename=cat_filename
         nrows=0
         print*,'Reading catalogue'

         open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
         IF (Reason > 0)  THEN
            print*,'Error opening file',filename
         endif

         print*,filename
         nrows=0
         ncat=0 !number of clusters above s/n or ylim threshold

         DO
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3, col4, col5, col6, col7, col8  ! :D !  z / zErr / SNR / y0 / SNR_{2.4} / y0Err / M500unc / M500upp
            !print*,Reason
            IF (Reason > 0)  THEN
               print*,'Error in reading file'
               print*,filename
               stop
            ELSE IF (Reason < 0) THEN
               exit
            ELSE
               nrows = nrows + 1
               !ncat = ncat + 1
               if (col5 >= snr_act_cat) ncat = ncat + 1  ! :D ! orginal cat (col3: snr_act_cat or col4: ylim_act_cat) / extra cat (col5: snr_act_cat)
               !if (col3 >= ylim_act_cat) ncat=ncat+1
            END IF
         END DO

         CLOSE (unit=iunit)
         print*, 'intrinsic scatter = ', cosmopar%d_act
         print*,'done'
         Print*,'Catalogue Number of clusters=', nrows                                                               ! need to check still WRONG!
         !       Print*,'Catalogue Number of clusters above the ylim threshold of', ylim_act_cat, '*10^(-4)=', ncat          ! orginal cat (ylim cut)
         !       Print*,'Catalogue Number of clusters above the SNR of', snr_act_cat, ' = ', ncat                               ! orginal cat (SNR cut)
         Print*,'Catalogue Number of clusters above the SNR2.4 of', snr_act_cat, ' = ', ncat                             ! extra cat

         if (allocated(SZcat)) deallocate(SZcat)
         ALLOCATE(SZcat(ncat,8),stat=iostat)         ! :D ! need to change when a different catalogue is used
         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit,file=filename,status='old',form="formatted")
         ii=1
         DO i=1,nrows
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3, col4, col5, col6, col7, col8  ! :D ! original cat (4 col) / extra cat (5 col)
            if (col5 >= snr_act_cat) then                   !!!!!!!! :D ! orginal cat (col3: snr_act_cat or col4: ylim_act_cat) / extra cat (col5: snr_act_cat)
               !if (col3 >= ylim_act_cat) then
               SZcat(ii,1)=col1  !redshift
               SZcat(ii,2)=col2  !error on redshift
               SZcat(ii,3)=col3  !detection S/N
               SZcat(ii,4)=col4  !y0
               SZcat(ii,5)=col5  !S/N2.4              ! :D ! from extra catalogue
               !SZcat(ii,6)=col6  !ylim from the map   ! :D ! from new catalogue (84)
               SZcat(ii,6)=col6  !y0err                ! from mass cat
               SZcat(ii,7)=col7  !M500unc              ! from mass cat
               SZcat(ii,8)=col8  !M500upp              ! from mass cat
               ii=ii+1
            endif
         ENDDO
         CLOSE (unit=iunit)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! load files describing ylim map

!         if (ylim_switch_act == 2) then

            filename = ylims_filename
            nrows = File%TxtFileLines(filename)
            print*, 'Number of sky patches =', nrows

!            surveypar%deg2_a = surveypar%deg2_a/nrows
            !          print*, surveypar%deg2_a

!            if (allocated(y0lim_arr)) deallocate(y0lim_arr)
!            allocate(y0lim_arr(nrows), stat=iostat)

            if (allocated(y0lims)) deallocate(y0lims)
            allocate(y0lims(nrows), stat=iostat)


!            if (allocated(noise_act_arr)) deallocate(noise_act_arr)
!            if (allocated(surveydeg2_act)) deallocate(surveydeg2_act)
!            allocate(noise_act_arr(nrows), stat=iostat)
!            allocate(surveydeg2_act(nrows), stat=iostat)

            if (iostat/=0) then
               print*,'allocation error'
            endif

            open (newunit=iunit, file=filename, status='unknown', form="formatted")
            DO i = 1, nrows
               READ(iunit,*,IOSTAT=Reason) col1!, col2
!               y0lim_arr(i) = col1
                y0lims(i) = col1
!               surveydeg2_act(i) = col1*3.046174198d-4
!               noise_act_arr(i)  = col2
            ENDDO
            CLOSE (unit=iunit)

!         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (allocated(z)) deallocate(z)
         ALLOCATE(Z(Nz), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         ! z vector
         DO i = 0, Nz-1
            z(i+1) = z0 + i*dz + 0.5_dl*dz
         END DO
         if (z0==0._dl) z(1) = z(1) + 1.e-8 ! for numerical problem when starting from 0.

         select case (dim_switch_act)
         case(1) ! ACT 1D N(z)

            if (allocated(delNcat)) deallocate(delNcat)
            if (allocated(delN)) deallocate(delN)
            allocate(delNcat(Nz), delN(Nz), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            delNcat(:) = 0.

            nrows=size(SZcat(:,1))

            nmiss=0
            DO ii=1,nrows
               if (SZcat(ii,1) <0.) nmiss=nmiss+1.
            enddo

            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz
               DO ii=1,nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                     delNcat(I)=delNcat(I)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            sum=0.
            DO I=1,Nz
               sum=sum+delNcat(I)
            ENDDO

            nred2=nrows-nmiss
            ncat=nrows

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss
            print*,'Counts here:',delNcat

            if (nmiss==0) nmiss_switch=0

            SELECT CASE(nmiss_switch)
            CASE(0)
               print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
            CASE(1)
               print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
            END SELECT


         case(2)  ! ACT 2D N(z, y0)

            ! load files describing selection function
            if (allocated(y0)) deallocate(y0)
            ALLOCATE(y0(Ny0), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            select case(y0bin)
            case(2)
               y0(1) = 0.8_dl
               y0(2) = 1.5_dl
            case(3)
               y0(1) = 0.8_dl
               y0(2) = 1.25_dl
               y0(3) = 1.7_dl
            end select
            print*, 'y0 = ', y0

            if (allocated(delN2Dcat)) deallocate(delN2Dcat)
            allocate(delN2Dcat(Nz, Ny0), delN2D(Nz, Ny0), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            delN2Dcat(:,:) = 0.

            nrows=size(SZcat(:,1))

            ! exclude clusters below s/n threshold q

            nmiss=0
            DO ii = 1, nrows
               if (SZcat(ii,1) <0.) nmiss=nmiss+1.
            enddo

            z_min = z0
            z_max = z_min + dz

            select case(y0bin)
            case (2)
               DO i = 1, Nz
                  DO ii = 1, nrows
                     if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                     (SZcat(ii,4) < y0(2)) .and. (SZcat(ii,4) >= y0(1))) then
                        delN2Dcat(i,1) = delN2Dcat(i,1) + 1.
                     end if
                  ENDDO
                  do ii = 1, nrows
                     if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,4) >= y0(2))) then
                        delN2Dcat(i,2) = delN2Dcat(i,2) + 1.
                     end if
                  end do
                  z_min=z_min+dz
                  z_max=z_max+dz
               END DO

               print*, 'Catalogue Counts:'
               sum=0.
               DO i = 1, Nz
                  do j = 1, Ny0
                     sum = sum + delN2Dcat(i,j)
                  end do
                  print*, i, delN2Dcat(i,1), delN2Dcat(i,2)
                  ! print out predicted number counts
                  open(unit=1, file = "ACT_catN_2D_2bins_snr24.txt")
                  write(1, 187) i/10.0, delN2Dcat(i,1), delN2Dcat(i,2)
                  187          format(F5.1, F5.1, F5.1)
               ENDDO
               close(unit=1)

            case(3)
               do i = 1, Nz
                  do ii = 1, nrows
                     do j = 1, Ny0
                        if (j < Ny0) then
                           if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                           (SZcat(ii,4) < y0(j+1)) .and. (SZcat(ii,4) >= y0(j))) then
                              delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                           end if
                        else if (j == Ny0) then
                           if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,4) >= y0(j))) then
                              delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                           end if
                        end if
                     end do
                  end do
                  z_min=z_min+dz
                  z_max=z_max+dz
               end do

               print*, 'Catalogue Counts:'
               sum=0.
               DO i = 1, Nz
                  do j = 1, Ny0
                     sum = sum + delN2Dcat(i,j)
                  end do
                  !print*, i, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3)
                  ! print out predicted number counts
                  !             open(unit=1, file = "ACT_catN_2D_3bins.txt")
                  !             write(1, 102) i/10.0, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3)
                  !102          format(F5.1, F5.1, F5.1, F5.1)
               ENDDO
               !          close(unit=1)


            end select
            print*, 'Total Catalogue Number', sum

            do i = 1, Nz
               sum = 0.
               do j = 1, Ny0
                  sum = sum + delN2Dcat(i,j)
               end do
               print*, i, sum
            end do

            print*, 'in each y0 bin:'

            do j = 1, Ny0
               sum = 0.
               do i = 1, Nz
                  sum = sum + delN2Dcat(i,j)
               end do
               print*, j, sum
            end do

            nred2=nrows-nmiss
            ncat=nrows

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss

         end select

      !------------------------------------------------------------------------!

      case(3) ! :D survey_switch 3 = both Planck and ACT

         ! Planck first

         SELECT CASE(dim_switch_plc)
         CASE(1)

            ! planck N(z)

            !file names
            cat_filename='data/SZ_cat.txt'
            !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
            thetas_filename='data/SZ_thetas.txt'
            skyfracs_filename='data/SZ_skyfracs.txt'
            ylims_filename='data/SZ_ylims.txt'

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2 = 41253.0  !full sky (sky fraction handled in skyfracs file)

            z0 = 0.0d0
            zmaxs = 1.
            dz = 0.1d0

            surveypar%ylimin=1.e-3
            surveypar%ymaxin=1.e-1
            logymin=0.7 !s/n=6
            logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
            dlogy=0.25

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            ylim_switch=-1  ! >0=constant ylim, <0=variable ylim
            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2 = 3.046174198d-4*surveypar%deg2
            Nz = DINT((zmaxs - z0)/dz)+1
            Ny = DINT((logymax-logymin)/dlogy)+1
            print*,'Nz=', Nz


            !ylims file                    for the variable ylim case
            if (ylim_switch < 0) then

               ymin=-1.*ymin

               ! read catalogue
               filename=cat_filename

               nrows = 0
               print*,'Reading catalogue'

               open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
               IF (Reason > 0)  THEN
                  print*,'Error opening file',filename
               endif

               print*,filename
               nrows = 1
               ncat = 1 !number of clusters above s/n threshold

               DO
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  !print*,Reason
                  IF (Reason > 0)  THEN
                     print*,'Error in reading file'
                     print*,filename
                     stop
                  ELSE IF (Reason < 0) THEN
                     exit
                  ELSE
                     nrows = nrows + 1
                     if (col3 >= q) ncat = ncat + 1
                  END IF
               END DO

               CLOSE (unit=iunit)
               print*,'done'

               Print*,'Catalogue Number of clusters=', nrows
               Print*,'Catalogue Number of clusters above the S/N threshold of', q, '=', ncat

               if (allocated(SZcat)) deallocate(SZcat)
               allocate(SZcat(ncat, 3), stat=iostat)! z,y
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='old',form="formatted")
               ii=1
               DO i=1,nrows
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  if (col3 >= q) then
                     SZcat(ii,1)=col1  !redshift
                     SZcat(ii,2)=col2  !error on redshift
                     SZcat(ii,3)=col3  !detection S/N
                     ii=ii+1
                  endif
               ENDDO
               CLOSE (unit=iunit)

               ! load files describing selection function

               !file with theta
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(thetas_filename, thetas, n=nthetas)
               print*,'Number of size thetas=',nthetas

               !file with skyfracs
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
               print*,'Number of patches=',npatches

               filename=ylims_filename
               nrows=File%TxtFileLines(filename)
               print*,'Number of size y =', nrows
               if (nrows /= npatches*nthetas) then
                  print*,'Format error for ylims.txt:'
                  print*,'Expected rows:',npatches*nthetas
                  print*,'Actual rows:',nrows
                  stop
               endif
               allocate(ylims(npatches,nthetas),stat=iostat)
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='unknown',form="formatted")
               i=1
               j=1
               DO ii=1, nrows
                  READ(iunit,*,IOSTAT=Reason)  col1
                  ylims(i,j)=col1
                  i=i+1
                  if (i > npatches) then
                     i=1
                     j=j+1
                  endif
               ENDDO
               CLOSE (unit=iunit)
            endif

            if (allocated(z)) deallocate(z)
            if (allocated(logy)) deallocate(logy)
            allocate(z(Nz), logy(Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! logy vector
            yi=logymin+dlogy/2.
            DO I=1,Ny+1
               logy(I)=yi
               yi=yi+dlogy
            END DO

            print*,'q =',10.d0**logy

            ! z vector
            DO I=0,Nz-1
               Z(I+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) Z(1) = Z(1) + 1.e-8 ! for numerical problem when starting from 0.

            if (allocated(DNcat)) deallocate(DNcat)
            if (allocated(DNzcat)) deallocate(DNzcat)
            if (allocated(DNycat)) deallocate(DNycat)
            if (allocated(DN)) deallocate(DN)
            if (allocated(DNz)) deallocate(DNz)
            if (allocated(DNy)) deallocate(DNy)
            allocate(DNcat(Nz, Ny+1), DNzcat(Nz), DNycat(Ny), DN(Nz+1, Ny+1), DNz(Nz), DNy(Ny), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            DNcat(:,:)=0.
            DNzcat(:)=0.

            nrows = size(SZcat(:,1))
            ! exclude clusters below s/n threshold q
            nmiss = 0
            DO ii=1,nrows
               if (SZcat(ii,1) < 0.) nmiss = nmiss + 1.
            enddo

            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz
               DO ii=1,nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                     DNzcat(I) = DNzcat(I)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            sum=0.
            DO I=1,Nz
               sum = sum + DNzcat(I)
            ENDDO

            nred2 = nrows - nmiss
            !print*, nrows, nmiss

            ncat = nrows

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss
            print*,'Counts:',DNzcat !!!!! array of catalogue cluster numbercount for each bin (10 except for 0)

            if (nmiss == 0) nmiss_switch=0

            SELECT CASE(nmiss_switch)
            CASE(0)
               print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
            CASE(1)
               print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
            END SELECT

         CASE(2) ! 2D

            ! Planck 2D

            ! compute catalogue counts in z and q
            ! compute P(q|qm) once and for all

            !file names
            cat_filename='data/SZ_cat.txt'
            !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
            thetas_filename='data/SZ_thetas.txt'
            skyfracs_filename='data/SZ_skyfracs.txt'
            ylims_filename='data/SZ_ylims.txt'

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2 = 41253.0  !full sky (sky fraction handled in skyfracs file)

            z0 = 0.0d0
            zmaxs = 1.
            dz = 0.1d0

            surveypar%ylimin=1.e-3
            surveypar%ymaxin=1.e-1
            logymin=0.7 !s/n=6
            logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
            dlogy=0.25

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            ylim_switch=-1  ! >0=constant ylim, <0=variable ylim
            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2 = 3.046174198d-4*surveypar%deg2
            Nz = DINT((zmaxs - z0)/dz)+1
            Ny = DINT((logymax-logymin)/dlogy)+1 !!!!! Nq
            print*,'Nz=',Nz !!!!!!!!!
            print*,'Ny=',Ny

            !ylims file                    for the variable ylim case
            if (ylim_switch < 0) then

               ymin=-1.*ymin

               ! read catalogue
               filename=cat_filename

               nrows = 0

               print*,'Reading catalogue'

               open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
               IF (Reason > 0)  THEN
                  print*,'Error opening file',filename
               endif

               print*,filename
               nrows = 1
               ncat  = 1 !number of clusters above s/n threshold

               DO
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  !print*,Reason
                  IF (Reason > 0)  THEN
                     print*,'Error in reading file'
                     print*,filename
                     stop
                  ELSE IF (Reason < 0) THEN
                     exit
                  ELSE
                     nrows = nrows + 1
                     if (col3 >= q) ncat = ncat + 1
                  END IF
               END DO

               CLOSE (unit=iunit)
               print*,'done'

               Print*,'Catalogue Number of clusters=', nrows
               Print*,'Catalogue Number of clusters above the S/N threshold of', q, '=', ncat

               if (allocated(SZcat)) deallocate(SZcat)
               allocate(SZcat(ncat, 3), stat=iostat)

               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='old',form="formatted")
               ii=1
               DO i=1,nrows
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  if (col3 >= q) then
                     SZcat(ii,1)=col1  !redshift
                     SZcat(ii,2)=col2  !error on redshift
                     SZcat(ii,3)=col3  !detection S/N
                     ii=ii+1
                  endif
               ENDDO
               CLOSE (unit=iunit)

               ! load files describing selection function

               !file with theta
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(thetas_filename, thetas, n=nthetas)
               print*,'Number of size thetas=',nthetas

               !file with skyfracs
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
               print*,'Number of patches=',npatches

               filename=ylims_filename
               nrows=File%TxtFileLines(filename)
               print*,'Number of size y =', nrows
               if (nrows /= npatches*nthetas) then
                  print*,'Format error for ylims.txt:'
                  print*,'Expected rows:',npatches*nthetas
                  print*,'Actual rows:',nrows
                  stop
               endif
               allocate(ylims(npatches,nthetas),stat=iostat)
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='unknown',form="formatted")
               i=1
               j=1
               DO ii=1, nrows
                  READ(iunit,*,IOSTAT=Reason)  col1
                  ylims(i,j)=col1
                  i=i+1
                  if (i > npatches) then
                     i=1
                     j=j+1
                  endif
               ENDDO
               CLOSE (unit=iunit)
            endif

            if (allocated(z)) deallocate(z)
            if (allocated(logy)) deallocate(logy)
            allocate(z(Nz), logy(Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! logy vector
            yi=logymin+dlogy/2.
            DO I=1,Ny+1
               logy(I)=yi
               yi=yi+dlogy
            END DO

            print*,'q =',10.d0**logy

            ! z vector
            DO I=0,Nz-1
               Z(I+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) Z(1) = Z(1) + 1.e-8 ! for numerical problem when starting from 0.

            if (allocated(DNcat)) deallocate(DNcat)
            if (allocated(DN)) deallocate(DN)
            allocate(DNcat(Nz,Ny+1), DN(Nz+1,Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            DNcat(:,:)=0.
            nrows = size(SZcat(:,1))

            nmiss = 0
            DO ii = 1, nrows
               if (SZcat(ii,1) < 0.d0) nmiss = nmiss + 1.
            enddo

            z_min = z0
            z_max = z_min + dz
            DO I = 1, Nz
               DO J = 1, Ny
                  y_min = logY(J) - dlogy/2.
                  y_max = logY(J) + dlogy/2.
                  y_min = 10.d0**y_min
                  y_max = 10.d0**y_max
                  DO ii = 1, nrows
                     if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                     (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                        DNcat(I,J) = DNcat(I,J)+1.
                     endif
                  ENDDO
               ENDDO
               J = Ny + 1 ! the last bin contains all S/N greater than what in the previous bin
               y_min = y_max
               DO ii = 1, nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,3) >= y_min)) then
                     DNcat(I,J) = DNcat(I,J) + 1.
                  endif
               ENDDO
               z_min = z_min + dz
               z_max = z_max + dz
            END DO

            !missing redshifts
            DO J = 1, Ny
               y_min = logY(J) - dlogy/2.
               y_max = logY(J) + dlogy/2.
               y_min = 10.**y_min
               y_max = 10.**y_max
               DO ii = 1, nrows
                  if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                     norm = 0.
                     do jj=1, Nz
                        norm = norm + DNcat(jj,J)
                     enddo
                     DNcat(:,J) = DNcat(:,J)*(norm + 1.d0)/norm
                  endif
               ENDDO
            ENDDO
            J = Ny + 1 ! the last bin contains all S/N greater than what in the previous bin
            y_min = y_max
            DO ii = 1, nrows
               if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) >= y_min)) then
                  norm = 0.
                  do jj = 1, Nz
                     norm = norm + DNcat(jj,J)
                  enddo
                  DNcat(:,J) = DNcat(:,J)*(norm + 1.d0)/norm
               endif
            ENDDO
            !end missing redshifts

            sum = 0.
            !print*, 'Catalogue counts:'
            DO I = 1, Nz
               DO J = 1, Ny+1
                  sum = sum + DNcat(I,J)
               ENDDO
               !print*, i, DNcat(i,1), DNcat(i,2), DNcat(i,3), DNcat(i,4)
            END DO
            print*, 'total', sum

            nred2 = nrows - nmiss
            ncat = nrows

            print*,'Number of clusters:', ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss

         end select

         ! now ACT

         select case (dim_switch_act)
         case(1)

            ! :D ACT 1D

            !file names
            cat_filename = 'data/ACT_SZ_cat.txt' ! :D !
            ylims_filename = 'data/ACT_SZ_ylims.txt'

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2_a = 987.5 ! ACT survey area (E-D56 field) in unit of deg2

            z0=0.0d0
            zmaxs_a=1.4    ! :D !
            dz=0.1d0

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2_a = 3.046174198d-4*surveypar%deg2_a

            Nz_a = dint((zmaxs_a-z0)/dz)+1

            print*,'Nz = ', Nz_a

            ! read catalogue
            filename=cat_filename
            nrows_a = 0
            print*,'Reading catalogue'

            open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
            IF (Reason > 0)  THEN
               print*,'Error opening file',filename
            endif

            print*,filename

            ncat_a  = 0 !number of clusters above s/n or ylim threshold

            DO
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4, col5
               !print*,Reason
               IF (Reason > 0)  THEN
                  print*,'Error in reading file'
                  print*,filename
                  stop
               ELSE IF (Reason < 0) THEN
                  exit
               ELSE
                  nrows_a = nrows_a + 1
                  if (col5 >= snr_act_cat) ncat_a = ncat_a + 1 !!!!!
               END IF
            END DO

            CLOSE (unit=iunit)
            print*,'done'

            Print*,'Catalogue Number of clusters=', nrows_a
            Print*,'Catalogue Number of clusters above the SNR_2.4 of', snr_act_cat, ' is ', ncat_a !!!!!
            ALLOCATE(SZcat_a(ncat_a,5), stat=iostat)! z, y
            if (iostat/=0) then
               print*,'allocation error'
            endif

            open (newunit=iunit,file=filename,status='old',form="formatted")
            ii=1
            DO i=1,nrows_a
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4, col5
               if (col5 >= snr_act_cat) then                   ! :D !
                  SZcat_a(ii,1)=col1  !redshift
                  SZcat_a(ii,2)=col2  !error on redshift
                  SZcat_a(ii,3)=col3  !detection S/N
                  SZcat_a(ii,4)=col4  !y0 !!!!!
                  SZcat_a(ii,5)=col5  !detection S/N_2.4
                  ii=ii+1
               endif
            ENDDO
            CLOSE (unit=iunit)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! load files describing ylim map

            if (ylim_switch_act == 2) then

               filename = ylims_filename
               nrows = File%TxtFileLines(filename)
               print*, 'Number of sky patches =', nrows

               surveypar%deg2_a = surveypar%deg2_a/nrows
               !          print*, surveypar%deg2_a

               if (allocated(y0lim_arr)) deallocate(y0lim_arr)
               allocate(y0lim_arr(nrows), stat=iostat)

               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit, file=filename, status='unknown', form="formatted")
               DO i = 1, nrows
                  READ(iunit,*,IOSTAT=Reason) col1
                  y0lim_arr(i) = col1
               ENDDO
               CLOSE (unit=iunit)

            endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! load files describing selection function
            if (allocated(z_a)) deallocate(z_a)
            allocate(z_a(Nz_a), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! z vector
            DO i = 0, Nz_a-1
               z_a(i+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) z_a(1) = z_a(1) + 1.e-8 ! for numerical problem when starting from 0.

            if (allocated(delNcat)) deallocate(delNcat)
            if (allocated(delN)) deallocate(delN)
            allocate(delNcat(Nz_a), delN(Nz_a), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            delNcat(:) = 0.

            nrows_a = size(SZcat_a(:,1))

            ! exclude clusters below s/n threshold q
            nmiss_a = 0
            DO ii=1,nrows_a
               if (SZcat_a(ii,1) <0.) nmiss=nmiss+1.
            enddo

            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz_a
               DO ii=1,nrows_a
                  if ((SZcat_a(ii,1) >= z_min) .and. (SZcat_a(ii,1) < z_max)) then
                     delNcat(I)=delNcat(I)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            sum=0.
            DO I=1,Nz_a
               sum=sum+delNcat(I)
            ENDDO

            nred2_a = nrows_a - nmiss_a
            ncat_a = nrows_a

            print*,'Number of clusters:', ncat_a
            print*,'Number of clusters with redshift:',nred2_a
            print*,'Number of clusters with no redshift:',nmiss_a
            print*,'Counts:',delNcat

            if (nmiss_a == 0) nmiss_switch=0

            SELECT CASE(nmiss_switch)
            CASE(0)
               print*,'Rescaling for missing redshifts',dble(nrows_a)/dble(nred2_a)
            CASE(1)
               print*,'Randomizing for missing redshifts',dble(nrows_a)/dble(nred2_a)
            END SELECT

         case(2)

            ! ACT 2D

            !file names
            cat_filename='data/ACT_SZ_cat.txt' ! :D !

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2_a = 987.5 ! ACT survey area (E-D56 field) in unit of deg2

            z0=0.0d0
            zmaxs_a=1.4    ! :D !
            dz=0.1d0

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2_a = 3.046174198d-4*surveypar%deg2_a

            Nz_a = dint((zmaxs_a-z0)/dz)+1
            Ny0 = y0binN

            print*,'Nz = ', Nz_a
            print*,'Ny0 = ', Ny0

            ! read catalogue
            filename=cat_filename
            nrows_a = 0
            print*,'Reading catalogue'

            open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
            IF (Reason > 0)  THEN
               print*,'Error opening file',filename
            endif

            print*,filename
            nrows_a = 1
            ncat_a  = 1 !number of clusters above s/n or ylim threshold

            DO
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4
               !print*,Reason
               IF (Reason > 0)  THEN
                  print*,'Error in reading file'
                  print*,filename
                  stop
               ELSE IF (Reason < 0) THEN
                  exit
               ELSE
                  nrows_a = nrows_a + 1
                  if (col4 >= ylim_act_cat) ncat_a = ncat_a + 1 !!!!!
               END IF
            END DO

            CLOSE (unit=iunit)
            print*,'done'

            Print*,'Catalogue Number of clusters=', nrows_a-1
            Print*,'Catalogue Number of clusters above the Ylim threshold of', ylim_act_cat, '*10^(-4)=', ncat_a !!!!!

            if (allocated(SZcat_a)) deallocate(SZcat_a)
            allocate(SZcat_a(ncat_a,4), stat=iostat)! z, y

            if (iostat/=0) then
               print*,'allocation error'
            endif

            open (newunit=iunit,file=filename,status='old',form="formatted")
            ii=1
            DO i=1,nrows_a
               READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4
               if (col4 >= ylim_act_cat) then                   ! :D !
                  SZcat_a(ii,1)=col1  !redshift
                  SZcat_a(ii,2)=col2  !error on redshift
                  SZcat_a(ii,3)=col3  !detection S/N
                  SZcat_a(ii,4)=col4  !y0 !!!!!
                  ii=ii+1
               endif
            ENDDO
            CLOSE (unit=iunit)

            ! load files describing selection function
            if (allocated(z_a)) deallocate(z_a)
            if (allocated(y0)) deallocate(y0)
            allocate(z_a(Nz_a), y0(Ny0), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! z vector
            DO i = 0, Nz_a-1
               z_a(i+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) z_a(1) = z_a(1) + 1.e-8 ! for numerical problem when starting from 0.

            ! y0 vector
            select case(y0bin)
            case(2)
               y0(1) = 0.8_dl
               y0(2) = 1.5_dl
            case(3)
               y0(1) = 0.8_dl
               y0(2) = 1.25_dl
               y0(3) = 1.7_dl
            end select
            print*, 'y0 = ', y0

            if (allocated(delN2Dcat)) deallocate(delN2Dcat)
            if (allocated(delN2D)) deallocate(delN2D)
            allocate(delN2Dcat(Nz_a, Ny0), delN2D(Nz_a, Ny0), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            delN2Dcat(:,:) = 0.

            nrows_a=size(SZcat_a(:,1))

            ! exclude clusters below s/n threshold q

            nmiss_a=0
            DO ii = 1, nrows_a
               if (SZcat_a(ii,1) <0.) nmiss_a = nmiss_a + 1.
            enddo

            z_min = z0
            z_max = z_min + dz

            select case(y0bin)
            case (2)
               DO i = 1, Nz_a
                  DO ii = 1, nrows_a
                     if ((SZcat_a(ii,1) >= z_min) .and. (SZcat_a(ii,1) < z_max) .and.&
                     (SZcat_a(ii,4) < y0(2)) .and. (SZcat_a(ii,4) >= y0(1))) then
                        delN2Dcat(i,1) = delN2Dcat(i,1) + 1.
                     end if
                  ENDDO
                  do ii = 1, nrows_a
                     if ((SZcat_a(ii,1) >= z_min) .and. (SZcat_a(ii,1) < z_max) .and. (SZcat_a(ii,4) >= y0(2))) then
                        delN2Dcat(i,2) = delN2Dcat(i,2) + 1.
                     end if
                  end do
                  z_min = z_min + dz
                  z_max = z_max + dz
               END DO

               print*, 'Catalogue Counts:'
               sum=0.
               DO i = 1, Nz_a
                  do j = 1, Ny0
                     sum = sum + delN2Dcat(i,j)
                  end do
                  print*, i, delN2Dcat(i,1), delN2Dcat(i,2)
                  ! print out predicted number counts
                  !                open(unit=1, file = "ACT_catN_2D_2bins_snr24.txt")
                  !                write(1, 101) i/10.0, delN2Dcat(i,1), delN2Dcat(i,2)
                  !101             format(F5.1, F5.1, F5.1)
               ENDDO
               !             close(unit=1)

               ! 2bins end ------------------------------------------------------

            case(3)
               do i = 1, Nz_a
                  do ii = 1, nrows_a
                     do j = 1, Ny0
                        if (j < Ny0) then
                           if ((SZcat_a(ii,1) >= z_min) .and. (SZcat_a(ii,1) < z_max) .and.&
                           (SZcat_a(ii,4) < y0(j+1)) .and. (SZcat_a(ii,4) >= y0(j))) then
                              delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                           end if
                        else if (j == Ny0) then
                           if ((SZcat_a(ii,1) >= z_min) .and. (SZcat_a(ii,1) < z_max) .and. (SZcat_a(ii,4) >= y0(j))) then
                              delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                           end if
                        end if
                     end do
                  end do
                  z_min = z_min + dz
                  z_max = z_max + dz
               end do

               print*, 'Catalogue Counts:'
               sum=0.
               DO i = 1, Nz_a
                  do j = 1, Ny0
                     sum = sum + delN2Dcat(i,j)
                  end do
                  print*, i, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3)
                  ! print out predicted number counts
                  !                   open(unit=1, file = "ACT_catN_2D_3bins.txt")
                  !                   write(1, 102) i/10.0, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3)
                  !102                format(F5.1, F5.1, F5.1, F5.1)
               ENDDO
               !                close(unit=1)

               ! 3bins end ------------------------------------------------------

            end select

            print*, 'Total Catalogue Number', sum

            do j = 1, Ny0
               sum = 0.
               do i = 1, Nz_a
                  sum = sum + delN2Dcat(i,j)
               end do
               print*, j, sum
            end do

            nred2_a = nrows_a - nmiss_a
            ncat_a = nrows_a

            print*,'Number of clusters:',ncat_a
            print*,'Number of clusters with redshift:',nred2_a
            print*,'Number of clusters with no redshift:',nmiss_a

         end select

      !------------------------------------------------------------------------!

      case(4) ! :D SPT

         !file names
         !cat_filename='data/SPT_SZ_cat.txt'                                    ! :D ! redshit / redshift_err / xi / Y_sz
         !cat_filename='data/SPT_SZ_cat_mass.txt'                               ! :D ! redshit / redshift_err / xi / Y_sz / mass (677)
         !cat_filename='data/SPT_Bocquet18_with_mass.txt'                        ! :D ! redshit / redshift_unc / xi / field / mass (677)
         cat_filename='data/SPT_Bocquet18_with_updated_redshift.txt'

         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         surveypar%deg2_s = 2500. ! SPT survey area in unit of deg2
         surveypar%deg2_s = 3.046174198d-4*surveypar%deg2_s


         z0 = 0.0d0
         zmaxs = 1.7    ! :D !
         dz = 0.1d0

         logximin = 0.6 ! snr = 4.5
         logximax = 1.7 ! snr = 42.35
         dlogxi = 0.3

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         nmiss_switch=0
         errorz_switch=0

         Nz = dint((zmaxs - z0)/dz) + 1
         Nxi = dint((logximax - logximin)/dlogxi) + 1
         print*,'Nz= ', Nz
         print*, 'Nxi= ', Nxi

         ! read catalogue
         filename = cat_filename

         nrows = 0
         print*,'Reading catalogue'

         open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
         IF (Reason > 0)  THEN
            print*,'Error opening file',filename
         endif

         print*,filename
         nrows = 0
         ncat = 0 !number of clusters above s/n or ylim threshold

         DO
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3, col4, col5
            !print*,Reason
            IF (Reason > 0)  THEN
               print*,'Error in reading file'
               print*,filename
               stop
            ELSE IF (Reason < 0) THEN
               exit
            ELSE
               nrows = nrows + 1
               if ((col3 >= snr_spt_cat) .and. (col1 > 0.25)) ncat = ncat + 1   ! :D ! cosmological sample
               !if ((col3 >= snr_spt_cat) .and. (col1 > 0.)) ncat = ncat + 1
            END IF
         END DO

         CLOSE (unit=iunit)
         print*,'done'

         Print*,'Catalogue Number of clusters=', nrows
         Print*,'Catalogue Number of clusters above the S/N threshold of', snr_spt_cat, '=', ncat ! :D !

         if (allocated(SZcat)) deallocate(SZcat)
         allocate(SZcat(ncat,5),stat=iostat)

         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit,file=filename,status='old',form="formatted")
         ii=1
         DO i=1,nrows
            READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4, col5
            if ((col3 >= snr_spt_cat)  .and. (col1 > 0.25))  then               ! :D ! cosmological sample
               !if ((col3 >= snr_spt_cat) .and. (col1 > 0.))  then
               SZcat(ii,1)=col1  !redshift
               SZcat(ii,2)=col2  !error on redshift
               SZcat(ii,3)=col3  !detection S/N : xi
               SZcat(ii,4)=col4  !field name
               SZcat(ii,5)=col5  !M500
               ii=ii+1
            endif
         ENDDO
         CLOSE (unit=iunit)

         ! load files describing selection function
         if (allocated(z)) deallocate(z)
         if (allocated(logxi)) deallocate(logxi)
         allocate(z(Nz), logxi(Nxi+1), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         ! z vector
         do i = 0, Nz-1
            z(i+1) = z0 + i*dz + 0.5_dl*dz
         enddo
         if (z0==0._dl) z(1) = z(1) + 1.e-8 ! for numerical problem when starting from 0.

         ! logxi vector
         xii = logximin + dlogxi/2.
         do i = 1, Nxi+1
            logxi(i) = xii
            xii = xii + dlogxi
         enddo

         print*, 'xi = ', 10.d0**logxi

         if (allocated(delN2Dcat)) deallocate(delN2Dcat)
         if (allocated(delNcat)) deallocate(delNcat)
         if (allocated(delN2D)) deallocate(delN2D)
         if (allocated(delN)) deallocate(delN)
         allocate(delN2Dcat(Nz, Nxi+1), delNcat(Nz), delN2D(Nz, Nxi+1), delN(Nz), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         delN2Dcat(:,:) = 0.
         delNcat(:) = 0.

         select case(dim_switch_spt)
         case(1) ! SPT 1D N(z)

             nrows=size(SZcat(:,1))

             ! exclude clusters below s/n threshold xi
             nmiss=0
             DO ii=1,nrows
                if (SZcat(ii,1) < 0.) nmiss=nmiss+1.
             enddo

             z_min=z0
             z_max=z_min+dz
             DO I=1,Nz
                DO ii=1,nrows
                   if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                      delNcat(I) = delNcat(I) + 1.
                   endif
                ENDDO
                z_min=z_min+dz
                z_max=z_max+dz
             END DO

             sum=0.
             DO I=1,Nz
                sum = sum + delNcat(I)
             ENDDO

             nred2 = nrows - nmiss
             !print*,nrows,nmiss

             ncat=nrows

             print*,'Catalogue S/N threshold:', snr_spt_cat
             print*,'Catalogue Number of clusters:', ncat
             print*,'Number of clusters with redshift:', nred2
             print*,'Number of clusters with no redshift:', nmiss
             print*,'Counts:', delNcat

             ! print out predicted number counts
             !       do i = 1, Nz
             !          open(unit=1, file = "SPT_catN.txt")
             !          write(1, 105) i/10.0, delNcat(i)
             !105       format(F5.1, F5.1)
             !       end do
             !       close(unit=1)

             if (nmiss==0) nmiss_switch=0

             SELECT CASE(nmiss_switch) ! it says 0 for simple scaling and 1 for MCMC
             CASE(0)
                print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
             CASE(1)
                print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
             END SELECT

         case(2) ! SPT 2D N(z,xi)

             nrows=size(SZcat(:,1))

             ! exclude clusters below s/n threshold xi
             nmiss=0
             DO ii=1,nrows
                if (SZcat(ii,1) < 0.) nmiss=nmiss+1.
             enddo

             z_min = z0
             z_max = z_min + dz
             DO I=1,Nz
                do j = 1, Nxi
                    xi_min = logxi(j) - dlogxi/2.
                    xi_max = logxi(j) + dlogxi/2.
                    xi_min = 10.d0**xi_min
                    xi_max = 10.d0**xi_max
                    do ii=1,nrows
                       if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                       (SZcat(ii,3) < xi_max) .and. (SZcat(ii,3) >= xi_min)) then
                          delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                       endif
                    enddo
                enddo
                j = Nxi + 1 ! the last bin constains all xi greater than what in the previous bin
                xi_min = xi_max
                do ii = 1, nrows
                    if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,3) >= xi_min)) then
                        delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                    endif
                enddo
                z_min = z_min + dz
                z_max = z_max + dz
             enddo

             do j = 1, Nxi + 1
                sum = 0.
                do i = 1, Nz
                    sum = sum + delN2Dcat(i,j)
                enddo
                print*, j, sum ! by signal-to-noise xi_cut
             enddo

             sum = 0.
             do i = 1, Nz
                do j = 1, Nxi + 1
                    sum = sum + delN2Dcat(i, j)
                enddo
                print*, i, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3), delN2Dcat(i,4), delN2Dcat(i,5)
             enddo
             print*, 'total = ', sum

             ! missing redshifts
             do j = 1, Nxi
                xi_min = logxi(j) - dlogxi/2.
                xi_max = logxi(j) + dlogxi/2.
                xi_min = 10.**xi_min
                xi_max = 10.**xi_max
                do ii = 1, nrows
                    if ((SZcat(ii,1) == 0) .and. (SZcat(ii,3) < xi_max) .and. (SZcat(ii,3) >= xi_min)) then
                        norm = 0.
                        do jj = 1, Nz
                            norm = norm + delN2Dcat(jj,j)
                        enddo
                        delN2Dcat(:,j) = delN2Dcat(:,j)*(norm + 1.d0)/norm
                    endif
                enddo
             enddo
             j = Nxi + 1
             xi_min = xi_max
             do ii = 1, nrows
                 if ((SZcat(ii,1) == 0) .and. (SZcat(ii,3) >= xi_min)) then
                     norm = 0.
                     do jj = 1, Nz
                         norm = norm + delN2Dcat(jj,j)
                     enddo
                     delN2Dcat(:,j) = delN2Dcat(:,j)*(norm + 1.d0)/norm
                 endif
             enddo

             nred2 = nrows - nmiss
             !print*,nrows,nmiss

             ncat = nrows

             print*,'Catalogue S/N threshold:', snr_spt_cat
             print*,'Catalogue Number of clusters:', ncat
             print*,'Number of clusters with redshift:', nred2
             print*,'Number of clusters with no redshift:', nmiss
             !print*,'Counts:', delNzcat

             ! print out predicted number counts
             !       do i = 1, Nz
             !          open(unit=1, file = "SPT_catN.txt")
             !          write(1, 105) i/10.0, delNcat(i)
             !105       format(F5.1, F5.1)
             !       end do
             !       close(unit=1)

             if (nmiss==0) nmiss_switch=0

             SELECT CASE(nmiss_switch) ! it says 0 for simple scaling and 1 for MCMC
             CASE(0)
                print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
             CASE(1)
                print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
             END SELECT

         end select


      !------------------------------------------------------------------------!

      case(5) ! :D survey_switch 5 = both Planck and SPT

         ! Planck first
         SELECT CASE(dim_switch_plc)
         CASE(1)

            ! planck N(z)

            !file names
            cat_filename='data/SZ_cat.txt'
            !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
            thetas_filename='data/SZ_thetas.txt'
            skyfracs_filename='data/SZ_skyfracs.txt'
            ylims_filename='data/SZ_ylims.txt'

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2 = 41253.0  !full sky (sky fraction handled in skyfracs file)

            z0 = 0.0d0
            zmaxs = 1.
            dz = 0.1d0

            surveypar%ylimin=1.e-3
            surveypar%ymaxin=1.e-1
            logymin=0.7 !s/n=6
            logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
            dlogy=0.25

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            ylim_switch=-1  ! >0=constant ylim, <0=variable ylim
            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2 = 3.046174198d-4*surveypar%deg2
            Nz = DINT((zmaxs - z0)/dz)+1
            Ny = DINT((logymax-logymin)/dlogy)+1
            print*,'Nz=', Nz


            !ylims file                    for the variable ylim case
            if (ylim_switch < 0) then

               ymin=-1.*ymin

               ! read catalogue
               filename=cat_filename

               nrows = 0
               print*,'Reading catalogue'

               open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
               IF (Reason > 0)  THEN
                  print*,'Error opening file',filename
               endif

               print*,filename
               nrows = 1
               ncat = 1 !number of clusters above s/n threshold

               DO
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  !print*,Reason
                  IF (Reason > 0)  THEN
                     print*,'Error in reading file'
                     print*,filename
                     stop
                  ELSE IF (Reason < 0) THEN
                     exit
                  ELSE
                     nrows = nrows + 1
                     if (col3 >= q) ncat = ncat + 1
                  END IF
               END DO

               CLOSE (unit=iunit)
               print*,'done'

               Print*,'Catalogue Number of clusters=', nrows
               Print*,'Catalogue Number of clusters above the S/N threshold of', q, '=', ncat

               if (allocated(SZcat)) deallocate(SZcat)
               allocate(SZcat(ncat, 3), stat=iostat)! z,y
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='old',form="formatted")
               ii=1
               DO i=1,nrows
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  if (col3 >= q) then
                     SZcat(ii,1)=col1  !redshift
                     SZcat(ii,2)=col2  !error on redshift
                     SZcat(ii,3)=col3  !detection S/N
                     ii=ii+1
                  endif
               ENDDO
               CLOSE (unit=iunit)

               ! load files describing selection function

               !file with theta
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(thetas_filename, thetas, n=nthetas)
               print*,'Number of size thetas=',nthetas

               !file with skyfracs
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
               print*,'Number of patches=',npatches

               filename=ylims_filename
               nrows=File%TxtFileLines(filename)
               print*,'Number of size y =', nrows
               if (nrows /= npatches*nthetas) then
                  print*,'Format error for ylims.txt:'
                  print*,'Expected rows:',npatches*nthetas
                  print*,'Actual rows:',nrows
                  stop
               endif
               allocate(ylims(npatches,nthetas),stat=iostat)
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='unknown',form="formatted")
               i=1
               j=1
               DO ii=1, nrows
                  READ(iunit,*,IOSTAT=Reason)  col1
                  ylims(i,j)=col1
                  i=i+1
                  if (i > npatches) then
                     i=1
                     j=j+1
                  endif
               ENDDO
               CLOSE (unit=iunit)
            endif

            if (allocated(z)) deallocate(z)
            if (allocated(logy)) deallocate(logy)
            allocate(z(Nz), logy(Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! logy vector
            yi=logymin+dlogy/2.
            DO I=1,Ny+1
               logy(I)=yi
               yi=yi+dlogy
            END DO

            print*,'q =',10.d0**logy

            ! z vector
            DO I=0,Nz-1
               Z(I+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) Z(1) = Z(1) + 1.e-8 ! for numerical problem when starting from 0.

            if (allocated(DNcat)) deallocate(DNcat)
            if (allocated(DNzcat)) deallocate(DNzcat)
            if (allocated(DNycat)) deallocate(DNycat)
            if (allocated(DN)) deallocate(DN)
            if (allocated(DNz)) deallocate(DNz)
            if (allocated(DNy)) deallocate(DNy)
            allocate(DNcat(Nz, Ny+1), DNzcat(Nz), DNycat(Ny), DN(Nz+1, Ny+1), DNz(Nz), DNy(Ny), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            DNcat(:,:)=0.
            DNzcat(:)=0.

            nrows = size(SZcat(:,1))
            ! exclude clusters below s/n threshold q
            nmiss = 0
            DO ii=1,nrows
               if (SZcat(ii,1) < 0.) nmiss = nmiss + 1.
            enddo

            z_min=z0
            z_max=z_min+dz
            DO I=1,Nz
               DO ii=1,nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                     DNzcat(I) = DNzcat(I)+1.
                  endif
               ENDDO
               z_min=z_min+dz
               z_max=z_max+dz
            END DO

            sum=0.
            DO I=1,Nz
               sum = sum + DNzcat(I)
            ENDDO

            nred2 = nrows - nmiss
            !print*, nrows, nmiss

            ncat = nrows

            print*,'Number of clusters:',ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss
            print*,'Counts:',DNzcat !!!!! array of catalogue cluster numbercount for each bin (10 except for 0)

            if (nmiss == 0) nmiss_switch=0

            SELECT CASE(nmiss_switch)
            CASE(0)
               print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
            CASE(1)
               print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
            END SELECT

         CASE(2) ! 2D

            ! Planck 2D

            ! compute catalogue counts in z and q
            ! compute P(q|qm) once and for all

            !file names
            cat_filename='data/SZ_cat.txt'
            !cat_filename='data/SZ_cat_intersection.txt' !intersection catalogue
            thetas_filename='data/SZ_thetas.txt'
            skyfracs_filename='data/SZ_skyfracs.txt'
            ylims_filename='data/SZ_ylims.txt'

            massfnpar%dso=500.
            surveypar%ab=0.
            surveypar%nb=0.
            surveypar%sfid=0.
            Nscat=0

            surveypar%deg2 = 41253.0  !full sky (sky fraction handled in skyfracs file)

            z0 = 0.0d0
            zmaxs = 1.
            dz = 0.1d0

            surveypar%ylimin=1.e-3
            surveypar%ymaxin=1.e-1
            logymin=0.7 !s/n=6
            logymax=1.5 !s/n=32. (higher s/n in the next to last bin)
            dlogy=0.25

            if (massfnpar%psind ==2) print*,'Using Tinker mass function'
            if (massfnpar%psind ==3) print*,'Using Watson mass function'

            ylim_switch=-1  ! >0=constant ylim, <0=variable ylim
            nmiss_switch=0
            errorz_switch=0

            surveypar%deg2 = 3.046174198d-4*surveypar%deg2
            Nz = DINT((zmaxs - z0)/dz)+1
            Ny = DINT((logymax-logymin)/dlogy)+1 !!!!! Nq
            print*,'Nz=',Nz !!!!!!!!!
            print*,'Ny=',Ny

            !ylims file                    for the variable ylim case
            if (ylim_switch < 0) then

               ymin=-1.*ymin

               ! read catalogue
               filename=cat_filename

               nrows = 0

               print*,'Reading catalogue'

               open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
               IF (Reason > 0)  THEN
                  print*,'Error opening file',filename
               endif

               print*,filename
               nrows = 1
               ncat  = 1 !number of clusters above s/n threshold

               DO
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  !print*,Reason
                  IF (Reason > 0)  THEN
                     print*,'Error in reading file'
                     print*,filename
                     stop
                  ELSE IF (Reason < 0) THEN
                     exit
                  ELSE
                     nrows = nrows + 1
                     if (col3 >= q) ncat = ncat + 1
                  END IF
               END DO

               CLOSE (unit=iunit)
               print*,'done'

               Print*,'Catalogue Number of clusters=', nrows
               Print*,'Catalogue Number of clusters above the S/N threshold of', q, '=', ncat

               if (allocated(SZcat)) deallocate(SZcat)
               allocate(SZcat(ncat, 3), stat=iostat)

               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='old',form="formatted")
               ii=1
               DO i=1,nrows
                  READ(iunit,*,IOSTAT=Reason)  col1, col2, col3
                  if (col3 >= q) then
                     SZcat(ii,1)=col1  !redshift
                     SZcat(ii,2)=col2  !error on redshift
                     SZcat(ii,3)=col3  !detection S/N
                     ii=ii+1
                  endif
               ENDDO
               CLOSE (unit=iunit)

               ! load files describing selection function

               !file with theta
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(thetas_filename, thetas, n=nthetas)
               print*,'Number of size thetas=',nthetas

               !file with skyfracs
               ! theta, nbins, first_index, last_index,first_index2, last_index2
               call File%LoadTxt(skyfracs_filename, skyfracs, n=npatches)
               print*,'Number of patches=',npatches

               filename=ylims_filename
               nrows=File%TxtFileLines(filename)
               print*,'Number of size y =', nrows
               if (nrows /= npatches*nthetas) then
                  print*,'Format error for ylims.txt:'
                  print*,'Expected rows:',npatches*nthetas
                  print*,'Actual rows:',nrows
                  stop
               endif
               allocate(ylims(npatches,nthetas),stat=iostat)
               if (iostat/=0) then
                  print*,'allocation error'
               endif

               open (newunit=iunit,file=filename,status='unknown',form="formatted")
               i=1
               j=1
               DO ii=1, nrows
                  READ(iunit,*,IOSTAT=Reason)  col1
                  ylims(i,j)=col1
                  i=i+1
                  if (i > npatches) then
                     i=1
                     j=j+1
                  endif
               ENDDO
               CLOSE (unit=iunit)
            endif

            if (allocated(z)) deallocate(z)
            if (allocated(logy)) deallocate(logy)
            allocate(z(Nz), logy(Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            ! logy vector
            yi=logymin+dlogy/2.
            DO I=1,Ny+1
               logy(I)=yi
               yi=yi+dlogy
            END DO

            print*,'q =',10.d0**logy

            ! z vector
            DO I=0,Nz-1
               Z(I+1) = z0 + I*dz + 0.5_dl*dz
            END DO
            if (z0==0._dl) Z(1) = Z(1) + 1.e-8 ! for numerical problem when starting from 0.

            if (allocated(DNcat)) deallocate(DNcat)
            if (allocated(DN)) deallocate(DN)
            allocate(DNcat(Nz,Ny+1), DN(Nz+1,Ny+1), stat=iostat)

            if (iostat /= 0) then
               print *, "Cannot allocate work arrays"
               stop
            endif

            DNcat(:,:)=0.
            nrows = size(SZcat(:,1))

            nmiss = 0
            DO ii = 1, nrows
               if (SZcat(ii,1) < 0.d0) nmiss = nmiss + 1.
            enddo

            z_min = z0
            z_max = z_min + dz
            DO I = 1, Nz
               DO J = 1, Ny
                  y_min = logY(J) - dlogy/2.
                  y_max = logY(J) + dlogy/2.
                  y_min = 10.d0**y_min
                  y_max = 10.d0**y_max
                  DO ii = 1, nrows
                     if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and.&
                     (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                        DNcat(I,J) = DNcat(I,J)+1.
                     endif
                  ENDDO
               ENDDO
               J = Ny + 1 ! the last bin contains all S/N greater than what in the previous bin
               y_min = y_max
               DO ii = 1, nrows
                  if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max) .and. (SZcat(ii,3) >= y_min)) then
                     DNcat(I,J) = DNcat(I,J) + 1.
                  endif
               ENDDO
               z_min = z_min + dz
               z_max = z_max + dz
            END DO

            !missing redshifts
            DO J = 1, Ny
               y_min = logY(J) - dlogy/2.
               y_max = logY(J) + dlogy/2.
               y_min = 10.**y_min
               y_max = 10.**y_max
               DO ii = 1, nrows
                  if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) < y_max) .and. (SZcat(ii,3) >= y_min)) then
                     norm = 0.
                     do jj=1, Nz
                        norm = norm + DNcat(jj,J)
                     enddo
                     DNcat(:,J) = DNcat(:,J)*(norm + 1.d0)/norm
                  endif
               ENDDO
            ENDDO
            J = Ny + 1 ! the last bin contains all S/N greater than what in the previous bin
            y_min = y_max
            DO ii = 1, nrows
               if ((SZcat(ii,1) == -1) .and. (SZcat(ii,3) >= y_min)) then
                  norm = 0.
                  do jj = 1, Nz
                     norm = norm + DNcat(jj,J)
                  enddo
                  DNcat(:,J) = DNcat(:,J)*(norm + 1.d0)/norm
               endif
            ENDDO
            !end missing redshifts

            sum = 0.
            !print*, 'Catalogue counts:'
            DO I = 1, Nz
               DO J = 1, Ny+1
                  sum = sum + DNcat(I,J)
               ENDDO
               !print*, i, DNcat(i,1), DNcat(i,2), DNcat(i,3), DNcat(i,4)
            END DO
            print*, 'total', sum

            nred2 = nrows - nmiss
            ncat = nrows

            print*,'Number of clusters:', ncat
            print*,'Number of clusters with redshift:',nred2
            print*,'Number of clusters with no redshift:',nmiss

         end select





         ! now :D SPT

         !file names
         !cat_filename='data/SPT_SZ_cat.txt'                                    ! :D ! redshit / redshift_err / xi / Y_sz
         !cat_filename='data/SPT_SZ_cat_mass.txt'                               ! :D ! redshit / redshift_err / xi / Y_sz / mass (677)
         !cat_filename='data/SPT_Bocquet18_with_mass.txt'                        ! :D ! redshit / redshift_unc / xi / field / mass (677)
         cat_filename='data/SPT_Bocquet18_with_updated_redshift.txt'

         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         surveypar%deg2_s = 2500. ! SPT survey area in unit of deg2
         surveypar%deg2_s = 3.046174198d-4*surveypar%deg2_s

         z0 = 0.0d0
         zmaxs_s = 1.7    ! :D !
         dz = 0.1d0

         logximin = 0.6 ! snr = 4.5
         logximax = 1.7 ! snr = 42.35
         dlogxi = 0.3

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         nmiss_switch=0
         errorz_switch=0



         Nz_s = dint((zmaxs_s - z0)/dz) + 1
         Nxi = dint((logximax - logximin)/dlogxi) + 1
         print*,'Nz= ', Nz_s
         print*, 'Nxi= ', Nxi

         ! read catalogue
         filename = cat_filename

         nrows_s = 0
         print*,'Reading catalogue'

         open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
         IF (Reason > 0)  THEN
            print*,'Error opening file',filename
         endif

         print*,filename
         nrows_s = 0
         ncat_s = 0 !number of clusters above s/n or ylim threshold

         DO
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3, col4, col5
            !print*,Reason
            IF (Reason > 0)  THEN
               print*,'Error in reading file'
               print*,filename
               stop
            ELSE IF (Reason < 0) THEN
               exit
            ELSE
               nrows_s = nrows_s + 1
               if ((col3 >= snr_spt_cat) .and. (col1 > 0.25)) ncat_s = ncat_s + 1   ! :D ! cosmological sample
               !if ((col3 >= snr_spt_cat) .and. (col1 > 0.)) ncat = ncat + 1
            END IF
         END DO

         CLOSE (unit=iunit)
         print*,'done'

         Print*,'Catalogue Number of clusters=', nrows_s
         Print*,'Catalogue Number of clusters above the S/N threshold of', snr_spt_cat, '=', ncat_s ! :D !

         if (allocated(SZcat_s)) deallocate(SZcat_s)
         allocate(SZcat_s(ncat_s,5),stat=iostat)

         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit,file=filename,status='old',form="formatted")
         ii=1
         DO i=1,nrows_s
            READ(iunit,*,IOSTAT=Reason)  col1, col2, col3, col4, col5
            if ((col3 >= snr_spt_cat)  .and. (col1 > 0.25))  then               ! :D ! cosmological sample
               !if ((col3 >= snr_spt_cat) .and. (col1 > 0.))  then
               SZcat_s(ii,1)=col1  !redshift
               SZcat_s(ii,2)=col2  !error on redshift
               SZcat_s(ii,3)=col3  !detection S/N : xi
               SZcat_s(ii,4)=col4  !field name
               SZcat_s(ii,5)=col5  !M500
               ii=ii+1
            endif
         ENDDO
         CLOSE (unit=iunit)

         ! load files describing selection function
         if (allocated(z_s)) deallocate(z_s)
         if (allocated(logxi)) deallocate(logxi)
         allocate(z_s(Nz_s), logxi(Nxi+1), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         ! z vector
         do i = 0, Nz_s-1
            z_s(i+1) = z0 + i*dz + 0.5_dl*dz
         enddo
         if (z0==0._dl) z_s(1) = z_s(1) + 1.e-8 ! for numerical problem when starting from 0.

         ! logxi vector
         xii = logximin + dlogxi/2.
         do i = 1, Nxi+1
            logxi(i) = xii
            xii = xii + dlogxi
         enddo

         print*, 'xi = ', 10.d0**logxi

         if (allocated(delN2Dcat)) deallocate(delN2Dcat)
         if (allocated(delNcat)) deallocate(delNcat)
         if (allocated(delN2D)) deallocate(delN2D)
         if (allocated(delN)) deallocate(delN)
         allocate(delN2Dcat(Nz_s, Nxi+1), delNcat(Nz_s), delN2D(Nz_s, Nxi+1), delN(Nz_s), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         delN2Dcat(:,:) = 0.
         delNcat(:) = 0.

         select case(dim_switch_spt)
         case(1) ! SPT 1D N(z)

             nrows_s=size(SZcat_s(:,1))

             ! exclude clusters below s/n threshold xi
             nmiss_s=0
             DO ii=1,nrows_s
                if (SZcat_s(ii,1) < 0.) nmiss_s=nmiss_s+1.
             enddo

             z_min=z0
             z_max=z_min+dz
             DO I=1,Nz_s
                DO ii=1,nrows_s
                   if ((SZcat_s(ii,1) >= z_min) .and. (SZcat_s(ii,1) < z_max)) then
                      delNcat(I) = delNcat(I) + 1.
                   endif
                ENDDO
                z_min=z_min+dz
                z_max=z_max+dz
             END DO

             sum=0.
             DO I=1,Nz_s
                sum = sum + delNcat(I)
             ENDDO

             nred2_s = nrows_s - nmiss_s
             !print*,nrows,nmiss

             ncat_s = nrows_s

             print*,'Catalogue S/N threshold:', snr_spt_cat
             print*,'Catalogue Number of clusters:', ncat_s
             print*,'Number of clusters with redshift:', nred2_s
             print*,'Number of clusters with no redshift:', nmiss_s
             print*,'Counts:', delNcat

             ! print out predicted number counts
             !       do i = 1, Nz
             !          open(unit=1, file = "SPT_catN.txt")
             !          write(1, 105) i/10.0, delNcat(i)
             !105       format(F5.1, F5.1)
             !       end do
             !       close(unit=1)

             if (nmiss_s==0) nmiss_switch=0

             SELECT CASE(nmiss_switch) ! it says 0 for simple scaling and 1 for MCMC
             CASE(0)
                print*,'Rescaling for missing redshifts',dble(nrows_s)/dble(nred2_s)
             CASE(1)
                print*,'Randomizing for missing redshifts',dble(nrows_s)/dble(nred2_s)
             END SELECT

         case(2) ! SPT 2D N(z,xi)

             nrows_s=size(SZcat_s(:,1))

             ! exclude clusters below s/n threshold xi
             nmiss_s=0
             DO ii=1,nrows_s
                if (SZcat_s(ii,1) < 0.) nmiss_s = nmiss_s + 1.
             enddo

             z_min = z0
             z_max = z_min + dz
             DO I = 1, Nz_s
                do j = 1, Nxi
                    xi_min = logxi(j) - dlogxi/2.
                    xi_max = logxi(j) + dlogxi/2.
                    xi_min = 10.d0**xi_min
                    xi_max = 10.d0**xi_max
                    do ii = 1, nrows_s
                       if ((SZcat_s(ii,1) >= z_min) .and. (SZcat_s(ii,1) < z_max) .and.&
                       (SZcat_s(ii,3) < xi_max) .and. (SZcat_s(ii,3) >= xi_min)) then
                          delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                       endif
                    enddo
                enddo
                j = Nxi + 1 ! the last bin constains all xi greater than what in the previous bin
                xi_min = xi_max
                do ii = 1, nrows_s
                    if ((SZcat_s(ii,1) >= z_min) .and. (SZcat_s(ii,1) < z_max) .and. (SZcat_s(ii,3) >= xi_min)) then
                        delN2Dcat(i,j) = delN2Dcat(i,j) + 1.
                    endif
                enddo
                z_min = z_min + dz
                z_max = z_max + dz
             enddo

             do j = 1, Nxi + 1
                sum = 0.
                do i = 1, Nz_s
                    sum = sum + delN2Dcat(i,j)
                enddo
                print*, j, sum ! by signal-to-noise xi_cut
             enddo

             sum = 0.
             do i = 1, Nz_s
                do j = 1, Nxi + 1
                    sum = sum + delN2Dcat(i, j)
                enddo
                print*, i, delN2Dcat(i,1), delN2Dcat(i,2), delN2Dcat(i,3), delN2Dcat(i,4), delN2Dcat(i,5)
             enddo
             print*, 'total = ', sum

             ! missing redshifts
             do j = 1, Nxi
                xi_min = logxi(j) - dlogxi/2.
                xi_max = logxi(j) + dlogxi/2.
                xi_min = 10.**xi_min
                xi_max = 10.**xi_max
                do ii = 1, nrows_s
                    if ((SZcat_s(ii,1) == 0) .and. (SZcat_s(ii,3) < xi_max) .and. (SZcat_s(ii,3) >= xi_min)) then
                        norm = 0.
                        do jj = 1, Nz_s
                            norm = norm + delN2Dcat(jj,j)
                        enddo
                        delN2Dcat(:,j) = delN2Dcat(:,j)*(norm + 1.d0)/norm
                    endif
                enddo
             enddo
             j = Nxi + 1
             xi_min = xi_max
             do ii = 1, nrows_s
                 if ((SZcat_s(ii,1) == 0) .and. (SZcat_s(ii,3) >= xi_min)) then
                     norm = 0.
                     do jj = 1, Nz_s
                         norm = norm + delN2Dcat(jj,j)
                     enddo
                     delN2Dcat(:,j) = delN2Dcat(:,j)*(norm + 1.d0)/norm
                 endif
             enddo

             nred2_s = nrows_s - nmiss_s
             !print*,nrows,nmiss

             ncat_s = nrows_s

             print*,'Catalogue S/N threshold:', snr_spt_cat
             print*,'Catalogue Number of clusters:', ncat_s
             print*,'Number of clusters with redshift:', nred2_s
             print*,'Number of clusters with no redshift:', nmiss_s
             !print*,'Counts:', delNzcat

             ! print out predicted number counts
             !       do i = 1, Nz
             !          open(unit=1, file = "SPT_catN.txt")
             !          write(1, 105) i/10.0, delNcat(i)
             !105       format(F5.1, F5.1)
             !       end do
             !       close(unit=1)

             if (nmiss_s==0) nmiss_switch=0

             SELECT CASE(nmiss_switch) ! it says 0 for simple scaling and 1 for MCMC
             CASE(0)
                print*,'Rescaling for missing redshifts',dble(nrows_s)/dble(nred2_s)
             CASE(1)
                print*,'Randomizing for missing redshifts',dble(nrows_s)/dble(nred2_s)
             END SELECT

         end select




      !------------------------------------------------------------------------!

      case(6) ! :D survey switch 6 = SO single tile test

         !file names
         !cat_filename='data/MFMF_SOSim_3freq_small_cat.txt'        ! :D ! (586) full catalogue with fixedSNR > 4   z / zErr / SNR
         cat_filename='data/MFMF_SOSim_3freq_tiles_cat.txt'        ! :D ! (29722) full catalogue with fixedSNR > 4   z / zErr / SNR    - fullmap test
         !cat_filename='data/SOSim_3freq_tiles_cat_SNR.txt'          ! :D ! (29722) full catalogue with SNR > 4         z / zErr / SNR   - fullmap test
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab.txt'          ! :D ! (5270) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab02_1070.txt'   ! :D ! (1070) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab03_495.txt'    ! :D ! (495) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab04_269.txt'    ! :D ! (269) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab08_75.txt'     ! :D ! (75) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab10_44.txt'     ! :D ! (44) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab12_24.txt'     ! :D ! (24) noise map with fixedSNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_005.txt'    ! :D ! (1046) stepsize = 0.05e-7 noise map with SNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_01.txt'     ! :D ! (632) stepsize = 0.1e-7 noise map with SNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_05.txt'     ! :D ! (150) stepsize = 0.5e-7 noise map with SNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_10.txt'     ! :D !  (79) stepsize = 1e-7 noise map with SNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_20.txt'     ! :D !  (41) stepsize = 2e-7 noise map with SNR > 5   skypatch / noise
         !ylims_filename='data/MFMF_SOSim_3freq_small_RMSTab_downsample_50.txt'     ! :D !  (17) stepsize = 5e-7 noise map with SNR > 5   skypatch / noise    before update
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_downsample_5.txt'         ! :D !  (2280) stepsize = 5e-7   skypatch / noise / tileName - fullmap test
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_10bin_downsample_10.txt' ! :D ! (1234)   stepsize = 10e-7  skypatch / noise / tileName - fullmap test
         ylims_filename='data/SOSim_3freq_tiles_RMSTab_downsample_50.txt'         ! :D ! (393)   stepsize = 50e-7  skypatch / noise / tileName - fullmap test
                                                                                        ! don't need tilename here anyway - it doesn't matter

         !Qfunc_filename='data/SOSim_3freq_small_Qfit.txt'            ! SO single tile test
         Qfunc_filename='data/SOSim_3freq_tiles_Qfit_mean.txt'        ! SO fullmap test with one average Q func

         print*, ylims_filename
         print*, Qfunc_filename

         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         !surveypar%deg2_SOtest = 599.35289
         surveypar%deg2_SOtest = 17826.8145007   ! fullmap test with one average Q func
         surveypar%deg2_SOtest = 3.046174198d-4*surveypar%deg2_SOtest

         z0 = 0.0d0
         zmaxs = 2.0   ! :D !
         dz = 0.1d0

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         nmiss_switch=0
         errorz_switch=0

         Nz = dint((zmaxs-z0)/dz)+1

         print*, 'Nz=', Nz

         ! read catalogue
         filename=cat_filename
         nrows=0
         print*,'Reading catalogue'

         open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
         IF (Reason > 0)  THEN
            print*,'Error opening file',filename
         endif

         print*, filename
         nrows=0
         ncat=0   !number of clusters above s/n threshold

         DO
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3
            !print*,Reason
            IF (Reason > 0)  THEN
               print*,'Error in reading file'
               print*,filename
               stop
            ELSE IF (Reason < 0) THEN
               exit
            ELSE
               nrows = nrows + 1
               if (col3 >= snr_so_cat) ncat = ncat + 1
            END IF
         END DO

         CLOSE (unit=iunit)
         print*,'done'
         Print*,'Catalogue Number of clusters = ', nrows
         Print*,'Catalogue Number of clusters above the SNR of', snr_so_cat, ' = ', ncat

         if (allocated(SZcat)) deallocate(SZcat)
         ALLOCATE(SZcat(ncat,3),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit,file=filename,status='old',form="formatted")
         ii=1
         DO i=1,nrows
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3
            if (col3 >= snr_so_cat) then
               SZcat(ii,1) = col1  !redshift
               SZcat(ii,2) = col2  !error on redshift
               SZcat(ii,3) = col3  !detection S/N
               ii=ii+1
            endif
         ENDDO
         CLOSE (unit=iunit)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! load files describing ylim map

         filename = ylims_filename
         nrows = File%TxtFileLines(filename)
         print*, 'Number of sky patches =', nrows

         if (allocated(noise_SOtest)) deallocate(noise_SOtest)
         if (allocated(skyfracs_SOtest)) deallocate(skyfracs_SOtest)
         allocate(noise_SOtest(nrows), skyfracs_SOtest(nrows), stat=iostat)

         ! surveypar%deg2_a = 3.046174198d-4*surveypar%deg2_a
         ! from deg2 to steradian (unit of solid angle)
         ! 3.046174198d-4 = (pi/180 deg)^2

         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit, file=filename, status='unknown', form="formatted")
         DO i = 1, nrows
            READ(iunit,*,IOSTAT=Reason) col1, col2
            skyfracs_SOtest(i) = col1*3.046174198d-4
            noise_SOtest(i)  = col2
         ENDDO
         CLOSE (unit=iunit)

         sum = 0.
         do i = 1, nrows
            sum = sum + skyfracs_SOtest(i)
         enddo

         print*, 'total survey area (deg2) = ', sum/3.046174198d-4

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         filename = Qfunc_filename
         nrows = File%TxtFileLines(filename)

         if (allocated(theta)) deallocate(theta)
         if (allocated(Q0)) deallocate(Q0)
         allocate(theta(nrows), Q0(nrows), stat=iostat)
         if (iostat/=0) then
            print*, 'allocation error'
         endif

         open (newunit=iunit, file=filename, status='unknown', form="formatted")
         do i = 1, nrows
            read(iunit, *, iostat=reason), col1, col2
            theta(i) = col1
            Q0(i) = col2
         enddo
         close(unit=iunit)




         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         if (allocated(z)) deallocate(z)
         allocate(z(Nz), stat=iostat)
         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         ! z vector
         DO i = 0, Nz-1
            z(i+1) = z0 + i*dz + 0.5_dl*dz
         END DO
         if (z0==0._dl) z(1) = z(1) + 1.e-8 ! for numerical problem when starting from 0.

         if (allocated(delNcat)) deallocate(delNcat)
         if (allocated(delN)) deallocate(delN)
         allocate(delNcat(Nz), delN(Nz), stat=iostat)

         if (iostat /= 0) then
            print *, "Cannot allocate work arrays"
            stop
         endif

         delNcat(:) = 0.

         nrows=size(SZcat(:,1))

         nmiss=0
         DO ii=1,nrows
            if (SZcat(ii,1) <0.) nmiss=nmiss+1.
         enddo

         z_min=z0
         z_max=z_min+dz
         DO I=1,Nz
            DO ii=1,nrows
               if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
                  delNcat(I)=delNcat(I)+1.
               endif
            ENDDO
            z_min=z_min+dz
            z_max=z_max+dz
         END DO

         sum=0.
         DO I=1,Nz
            sum=sum+delNcat(I)
         ENDDO
         nred2=nrows-nmiss
         ncat=nrows

         print*,'Number of clusters :',ncat
         print*,'Number of clusters with redshift:',nred2
         print*,'Number of clusters with no redshift:',nmiss
         print*,'Counts here:',delNcat

         if (nmiss==0) nmiss_switch=0

         SELECT CASE(nmiss_switch)
         CASE(0)
            print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
         CASE(1)
            print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
         END SELECT

      !------------------------------------------------------------------------!

      case(7) ! :D survey switch 7 = SO full map with 540 tiles - updated version has 346 tiles

         !file names
         cat_filename='data/MFMF_SOSim_3freq_tiles_cat.txt'                        ! :D ! (29722) full catalogue with fixedSNR > 4   z / zErr / SNR
         !cat_filename='data/SOSim_3freq_tiles_cat_SNR.txt'                        ! :D ! (29722) full catalogue with SNR > 4   z / zErr / SNR

         !ylims_filename='data/MFMF_SOSim_3freq_tiles_RMSTab_test.txt'             ! :D ! (202156) noise map with SNR > 5  skypatch / y0lim / tilename
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_downsample_5.txt'          ! :D ! (2280)   noise map with SNR > 5  skypatch / y0lim / tilename  ! updated v. 20190927 0.005e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_downsample_50.txt'         ! :D ! (393)    noise map with SNR > 5  skypatch / y0lim / tilename  ! updated v. 20190927 0.05e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_10bin_downsample_5.txt'  ! :D ! (2280) Qfunc is binned (10 bins) ! updated v. 20190927 0.005e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_10bin_downsample_10.txt' ! :D ! (1234) Qfunc is binned (10 bins) ! updated v. 20190927 0.01e-4
         ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_10bin_downsample_50.txt' ! :D ! (393)  Qfunc is binned (10 bins) ! updated v. 20190927 0.05e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_2bin_downsample_5.txt'   ! :D ! (2280) Qfunc is binned (2 bins)  ! updated v. 20190927 0.005e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_2bin_downsample_10.txt'  ! :D ! (1234) Qfunc is binned (2 bins)  ! updated v. 20190927 0.01e-4
         !ylims_filename='data/SOSim_3freq_tiles_RMSTab_Q_2bin_downsample_50.txt'  ! :D ! (393)  Qfunc is binned (2 bins)  ! updated v. 20190927 0.05e-4

         theta_filename='data/SOSim_3freq_tiles_theta_updated.txt'
         Qfunc_filename='data/SOSim_3freq_tiles_Qfit_10bins.txt'
         !Qfunc_filename='data/SOSim_3freq_tiles_Qfit_2bins.txt'

         print*, ylims_filename
         print*, Qfunc_filename

         massfnpar%dso=500.
         surveypar%ab=0.
         surveypar%nb=0.
         surveypar%sfid=0.
         Nscat=0

         surveypar%deg2_SOfull = 17826.8145007
         surveypar%deg2_SOfull = 3.046174198d-4*surveypar%deg2_SOfull

         z0 = 0.0d0
         zmaxs = 2.0   ! :D !
         dz = 0.1d0

         if (massfnpar%psind ==2) print*,'Using Tinker mass function'
         if (massfnpar%psind ==3) print*,'Using Watson mass function'

         nmiss_switch=0
         errorz_switch=0

         Nz = dint((zmaxs-z0)/dz)+1

         print*, 'Nz=', Nz

         ! read catalogue
         filename=cat_filename
         nrows=0
         print*,'Reading catalogue'

         open (newunit=iunit,file=filename,status='old',form="formatted",iostat=reason)
         IF (Reason > 0)  THEN
            print*,'Error opening file',filename
         endif

         print*,filename
         nrows=0
         ncat=0   !number of clusters above s/n threshold

         DO
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3
            !print*,Reason
            IF (Reason > 0)  THEN
               print*,'Error in reading file'
               print*,filename
               stop
            ELSE IF (Reason < 0) THEN
               exit
            ELSE
               nrows = nrows + 1
               if (col3 >= snr_so_cat) ncat = ncat + 1
            END IF
         END DO

         CLOSE (unit=iunit)
         print*,'done'
         Print*,'Catalogue Number of clusters = ', nrows
         Print*,'Catalogue Number of clusters above the SNR of', snr_so_cat, ' = ', ncat

         if (allocated(SZcat)) deallocate(SZcat)
         ALLOCATE(SZcat(ncat,3),stat=iostat)
         if (iostat/=0) then
            print*,'allocation error'
         endif

         open (newunit=iunit,file=filename,status='old',form="formatted")
         ii=1
         DO i=1,nrows
            READ(iunit,*,IOSTAT=Reason) col1, col2, col3
            if (col3 >= snr_so_cat) then
               SZcat(ii,1) = col1  !redshift
               SZcat(ii,2) = col2  !error on redshift
               SZcat(ii,3) = col3  !detection S/N
               ii=ii+1
            endif
         ENDDO
         CLOSE (unit=iunit)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! load files describing ylim map

         filename = ylims_filename
         nrows = File%TxtFileLines(filename)
         print*, 'Number of sky patches =', nrows

         if (allocated(noise_SOfull)) deallocate(noise_SOfull)
         if (allocated(skyfracs_SOfull)) deallocate(skyfracs_SOfull)
         if (allocated(tilenames)) deallocate(tilenames)

         allocate(noise_SOfull(nrows), stat=iostat)
         allocate(skyfracs_SOfull(nrows), stat=iostat)
         allocate(tilenames(nrows), stat=iostat)

         !surveypar%deg2_a = 3.046174198d-4*surveypar%deg2_a
         ! from deg2 to steradian (unit of solid angle)
         ! 3.046174198d-4 = (pi/180 deg)^2

         if (iostat/=0) then
         print*,'allocation error'
         endif

         open (newunit=iunit, file=filename, status='unknown', form="formatted")
         DO i = 1, nrows
         READ(iunit,*,IOSTAT=Reason) col1, col2, col3
         skyfracs_SOfull(i) = col1*3.046174198d-4
         noise_SOfull(i)  = col2
         tilenames(i) = col3
         ENDDO
         CLOSE (unit=iunit)

         sum = 0.
         do i = 1, nrows
         sum = sum + skyfracs_SOfull(i)
         end do

         print*,'total survey area (deg2) = ', sum/3.046174198d-4

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! load files describing Q func

         filename = theta_filename
         nrows = File%TxtFileLines(filename)

         if (allocated(theta)) deallocate(theta)
         allocate(theta(nrows), stat=iostat)
         if (iostat/=0) then
            print*, 'allocation error'
         endif

         open (newunit=iunit, file=filename, status='unknown', form="formatted")
         do i = 1, nrows
            read(iunit, *, iostat=reason), col1
            theta(i) = col1
         enddo
         close(unit=iunit)



         filename = Qfunc_filename
         nrows = File%TxtFileLines(filename)
         ncols = File%TxtFileColumns(filename)
         print*, 'Number of Qfuncs =', nrows, ' x ',  ncols

         if (allocated(Q00)) deallocate(Q00)
         allocate(Q00(nrows, ncols), stat=iostat)
         if (iostat/=0) then
            print*, 'allocation error'
         endif

         open (newunit=iunit, file=filename, status='unknown', form="formatted")
         do i = 1, nrows
            read(iunit, *, iostat=reason), (Q00(i,j), j=1,ncols)
         enddo
         close(unit=iunit)







         if (allocated(z)) deallocate(z)
         ALLOCATE(Z(Nz), stat=iostat)

         if (iostat /= 0) then
         print *, "Cannot allocate work arrays"
         stop
         endif

         ! z vector
         DO i = 0, Nz-1
         z(i+1) = z0 + i*dz + 0.5_dl*dz
         END DO
         if (z0==0._dl) z(1) = z(1) + 1.e-8 ! for numerical problem when starting from 0.


         if (allocated(delNcat)) deallocate(delNcat)
         if (allocated(delN)) deallocate(delN)
         allocate(delNcat(Nz), delN(Nz), stat=iostat)

         if (iostat /= 0) then
         print *, "Cannot allocate work arrays"
         stop
         endif

         delNcat(:) = 0.

         nrows=size(SZcat(:,1))

         nmiss=0
         DO ii=1,nrows
         if (SZcat(ii,1) <0.) nmiss=nmiss+1.
         enddo

         z_min=z0
         z_max=z_min+dz
         DO I=1,Nz
         DO ii=1,nrows
            if ((SZcat(ii,1) >= z_min) .and. (SZcat(ii,1) < z_max)) then
               delNcat(I)=delNcat(I)+1.
            endif
         ENDDO
         z_min=z_min+dz
         z_max=z_max+dz
         END DO

         sum=0.
         DO I=1,Nz
         sum=sum+delNcat(I)
         ENDDO
         nred2=nrows-nmiss
         ncat=nrows

         print*,'Number of clusters :',ncat
         print*,'Number of clusters with redshift:',nred2
         print*,'Number of clusters with no redshift:',nmiss
         print*,'Counts here:',delNcat

         if (nmiss==0) nmiss_switch=0

         SELECT CASE(nmiss_switch)
         CASE(0)
         print*,'Rescaling for missing redshifts',dble(nrows)/dble(nred2)
         CASE(1)
         print*,'Randomizing for missing redshifts',dble(nrows)/dble(nred2)
         END SELECT

      end select

      !------------------------------------------------------------------------!

      do_sz_init = .false.
      print*,'End SZ initialization'

   end subroutine SZ_init

   !---------------------------------------------------------------------------!

   function SZCC_Cash(this,CMB,Theory,DataParams)
      ! likelihood computation
      ! SZ nuisance in dataparams
      use cosmology
      use numbercounts
      use, intrinsic :: ieee_arithmetic
      Class(SZLikelihood) :: this
      Class(CMBParams):: CMB
      Class(TCosmoTheoryPredictions), target :: Theory
      real(mcp) DataParams(:)
      real(mcp) SZCC_Cash
      real(mcp) SZCC_Cash_p, SZCC_Cash_a, SZCC_Cash_s
      INTEGER :: N,i,j,Nf,ii
      REAL(DL) :: sum,factorial,ln_factorial,fact!,SZCC_Cash_exp

      SZCC_Cash   = logzero
      SZCC_Cash_p = logzero
      SZCC_Cash_a = logzero
      SZCC_Cash_s = logzero

      !Mapping of cosmo parameters
      cosmopar%H0=CMB%H0
      cosmopar%w0=CMB%w
      cosmopar%w1=0.0 ! Note, we do not support evolution of w right now
      cosmopar%omegam=CMB%omc+CMB%omb+CMB%omnu
      cosmopar%omegav=CMB%omv
      cosmopar%omegak=CMB%omk
      cosmopar%n=CMB%InitPower(ns_index)
      cosmopar%sig8=Theory%sigma_8
      cosmopar%omegabh2=CMB%ombh2
      cosmopar%gamma=-1
      !cosmopar%s8omegamp3=Theory%sigma_8*((CMB%omc+CMB%omb+CMB%omnu)/0.27)**0.3

      !Mapping of nuisance parameters
      cosmopar%alpha=DataParams(1)
      cosmopar%ystar=10.**(DataParams(2))/(2.**cosmopar%alpha)*0.00472724
      cosmopar%logystar=(DataParams(2))!to set the prior on logy
      cosmopar%bias=DataParams(3)
      cosmopar%biasinv=1./DataParams(3)
      cosmopar%sigmaM=DataParams(4)
      cosmopar%beta=DataParams(5)
      !-------------------------------------- :D
      cosmopar%ylim_act=DataParams(6)
      cosmopar%bias_act=DataParams(7)
      cosmopar%a_act=Dataparams(8)
      cosmopar%b_act=DataParams(9)
      cosmopar%d_act=DataParams(10)
      cosmopar%bias_spt=DataParams(11)
      cosmopar%a_spt=DataParams(12)
      cosmopar%b_spt=DataParams(13)
      cosmopar%c_spt=DataParams(14)
      cosmopar%d_spt=DataParams(15)
      !-------------------------------------- :D
      cosmopar%sigmaR=>Theory%sigma_R
      call INIGROWTH

      ! Planck
      DN(:,:)=0.
      DNz(:)=0.
      DNy(:)=0.

      ! :D  ACT and SPT and SO
      delN(:)=0.
      delN2D(:,:)=0.

      select case(survey_switch)
      !------------------------------------------------------------------------!
      case(1) ! Planck

         select case(dim_switch_plc)
         case(1) ! Planck N(z)

            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            DO I=1,Nz           ! 10 bins for z which is delta z = 0.1
               sum=sum+DNz(I)   ! integration dN(z) over z
               ! print out predicted number counts
               !             open(unit=1, file = "plc_dNdz_bias.txt")
               !             write(1, 200) I/10.0 - 0.1, DNz(I)
               !200          format(F5.1, F)
            ENDDO
            !          close(unit=1)

            if (print_counts==1) then
               print*,'prdicted counts'
               print*,DNz    !!!!!!!! N(z=1), N(z=2), N(z=3),...
               print*,'total predicted counts',sum
            endif

            if (ieee_is_nan(DNZ(1))) then
               print*,'NaN found in theory counts!'
               stop
            endif

            DNzcat = DNzcat/dble(nred2)*dble(ncat) !!!!! from catalogue
            !          print*,'catalogue counts =',  DNzcat

            !          DNz = DNz/0.9
            !          print*, 'number count divided by 0.9 =', DNz

            sum=0. !!!!! Planck 1D likelihood
            do i=1,Nz
               if (DNz(i) /= 0.) then
                  ln_factorial=0.
                  if (DNzcat(i) /= 0.) ln_factorial=0.918939+(DNzcat(i)+0.5)*dlog(DNzcat(i))-DNzcat(i) !Stirling
                     sum=sum-1.*(DNzcat(i)*dlog(DNz(i))-DNz(i)-ln_factorial)
               end if
            end do
            SZCC_Cash=sum

            DNzcat = DNzcat*dble(nred2)/dble(ncat)
            !          DNz = DNz*0.9
            !          print*, 'number count multiplied by 0.9 =', DNz

         case(2) ! Planck N(z,q)

            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            if (print_counts==1) print*,'predicted counts'
            DO I=1,Nz                                       !! 10 bins for z : delta z = 0.1
               if (print_counts==1) print*,I,DN(I,:)
               do J=1,Ny+1                                  !!  5 bins for q : delta log q = 0.25
                  sum=sum+DN(I,J)
               enddo
            ENDDO

            if (print_counts==1) print*,'total counts',sum

            !sum=0.
            !DO I=1,Nz
            !   sum=0.
            !   do J=1,Ny+1
            !      sum=sum+DN(I,J)
            !   enddo
            !   if (print_counts==1) print*, i, sum
            !ENDDO


            sum=0. !!!!! Planck 2D likelihood
            do i=1,Nz
               do j=1,Ny+1
                  if (DN(i,j) /= 0.) then
                     ln_factorial=0.
                     if (DNcat(i,j)/=0.) then
                        if (DNcat(i,j) >10.) then
                           ln_factorial=0.918939+(DNcat(i,j)+0.5)*dlog(DNcat(i,j))-DNcat(i,j)
                           !Stirling approximation only for more than 10 elements
                        else
                           !direct computation of factorial
                           fact=1.
                           do ii=1,int(DNcat(i,j))
                              fact=ii*fact
                           enddo
                           ln_factorial=dlog(fact)
                        endif
                     endif
                     sum=sum-1.*(DNcat(i,j)*dlog(DN(i,j))-DN(i,j)-ln_factorial)
                  end if
               enddo
            end do
            SZCC_Cash=sum

         end select

         print*,'SZCC_Cash', sum
         !Print*,'SZ lnlike = ',SZCC_Cash

         !priors for SZ nuisance params
         SZCC_Cash = SZCC_Cash + pystar*(cosmopar%logystar-(-0.186))**2./(2.*0.021**2.) + &
         palpha*(cosmopar%alpha-1.789)**2./(2.*0.084**2.) + psigma*(cosmopar%sigmaM-0.075)**2./(2.*0.01**2.) !cosmomc_sz
         !print*,'prior on ystar, alpha, scatter'

         SZCC_Cash = SZCC_Cash +pbeta*(cosmopar%beta-0.6666666)**2/(2.*0.5**2.) !prior beta evolution

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash = SZCC_Cash + lens*(cosmopar%biasinv-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias-0.780)**2./(2.*0.092**2.)!final

      !------------------------------------------------------------------------!

      case(2) ! ACT

         select case (dim_switch_act)
         case(1) ! ACT N(z)

            call deltaN_act(z, Nz, delN, y0lims)


!            if (ylim_switch_act == 1) then  ! constant ylim

!               call deltaN(Nz, z, delN) !!!

!            elseif (ylim_switch_act == 2) then ! ylim from the map

               !call Mlim_act_allocate(y0lim_arr)
!               call y0_act_allocate(y0_act_arr)

!!!               call deltaN_act_ymap(Nz, z, delN, y0lim_arr) !!!
!               call deltaN_act_ymap(Nz, z, delN, noise_act_arr, surveydeg2_act) !!!

               !call deltaN_mass(Nz, z, delN, SZcat) !! to calculate ACT mass estimates!!!! temp!!!

               !call Mlim_act_deallocate
!               call y0_act_deallocate

               !            sum=0.
               !            do i = 1, Nz
               !            	sum = sum + delN(i)
               !            end do
               !            print*, 'number counts', sum

               !             delN = delN/0.9

               !             sum=0.
               !             do i = 1, Nz
               !                sum = sum + delN(i)
               !             end do
               !             print*, 'considering 90 percent completeness', sum

!            endif

               !         delN = delN*0.9

            sum=0.
            do i = 1, Nz
               sum = sum + delN(i)
               !print out predicted number counts
               !              open(unit=1, file = "act_dNdz_B0.txt")
               !              write(1, 113) i/10.0 - 0.1, delN(i)
               !113           format(F5.1, F)
            end do
            !          close(unit=1)

            if (print_counts==1) then
               print*, 'prdicted counts'
               print*, delN
               print*, 'total counts = ', sum
            endif

            if (ieee_is_nan(delN(1))) then
               print*,'NaN found in theory counts!'
               stop
            endif

            !          delNcat = delNcat/dble(nred2)*dble(ncat) !!!!! from catalogue
            !          delNcat = delNcat/dble(0.9)

            ! considering 90% completeness !!!!!!!
            !          sum=0.
            !          do i = 1, Nz
            !             sum = sum + delN(i)/0.9
            !          end do
            !          print*, 'considering 90 percent completeness', sum

            !          delN = delN/0.9                       ! option 1

            sum=0. ! ACT 1D likelihood
            do i=1,Nz
               if (delN(i) /= 0.) then
                  ln_factorial=0.
                  if (delNcat(i) /= 0.) ln_factorial=0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
                  sum=sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
                  !                sum=sum-1.*(delNcat(i)*dlog(delN(i)/0.9)-(delN(i)/0.9)-ln_factorial)     ! option 2
               end if
            end do
            SZCC_Cash=sum

            !          delN = delN/0.9
            !          delNcat = delNcat*dble(0.9)

            print*,'SZCC_Cash ACT', sum
            !Print*,'SZ lnlike = ',SZCC_Cash


         case(2) ! ACT N(z, y0)

            call Mlim_act_allocate(y0lim_arr)
            call deltaN2D(Nz, z, Ny0, y0, delN2D)
            call Mlim_act_deallocate

            sum=0.
            if (print_counts == 1) print*, 'predicted counts'
            do i = 1, Nz
               !             if (print_counts == 1) print*, i, delN2D(i,:)
               do j = 1, Ny0
                  sum = sum + delN2D(i,j)
               end do
               ! print out predicted number counts *****
               !             open(unit=1, file = "ACT_theoryN_2D_SNR24.txt")
               !             write(1, 109) i/10.0, delN2D(i,1), delN2D(i,2)
               !109          format(F5.1, F10.5, F10.5)
            end do
            !          close(unit=1)

            ! just to see the counts
            do i = 1, Nz
               sum = 0.
               do j = 1, Ny0
                  sum = sum + delN2D(i,j)
               end do
               print*, i, sum
            end do

            if (print_counts==1) print*, ' in each y0 bin '

            do j = 1, Ny0
               sum = 0.
               do i = 1, Nz
                  sum = sum + delN2D(i,j)
               end do
               print*, j, sum
            end do

            sum = 0.
            do j = 1, Ny0
               do i = 1, Nz
                  sum = sum + delN2D(i,j)
               end do
            end do
            print*, 'total counts =', sum

            sum=0. ! ACT 2D likelihood
            do i = 1, Nz
               do j = 1, Ny0
                  if (delN2D(i,j) /= 0.) then
                     ln_factorial=0.
                     if (delN2Dcat(i,j) /= 0.) then
                        if (delN2Dcat(i,j) > 10.) then
                           ln_factorial = 0.918939 + (delN2Dcat(i,j) + 0.5)*dlog(delN2Dcat(i,j))-delN2Dcat(i,j) !Stirling approximation only for more than 10 elements
                        else
                           !direct computation of factorial
                           fact = 1.
                           do ii = 1, int(delN2Dcat(i,j))
                              fact = ii*fact
                           end do
                           ln_factorial = dlog(fact)
                        end if
                     end if
                     sum = sum - 1.*(delN2Dcat(i,j)*dlog(delN2D(i,j)) - delN2D(i,j) - ln_factorial)
                  end if
               end do
            end do
            SZCC_Cash=sum

            !print*,'SZCC_Cash', sum
            !Print*,'SZ lnlike = ',SZCC_Cash

         end select


         !priors for Planck SZ nuisance params
         !SZCC_Cash = SZCC_Cash + pystar*(cosmopar%logystar-(-0.186))**2./(2.*0.021**2.) + &
         !     palpha*(cosmopar%alpha-1.789)**2./(2.*0.084**2.) + psigma*(cosmopar%sigmaM-0.075)**2./(2.*0.01**2.)
         ! for planck coeff.
         !print*,'prior on ystar, alpha, sigmaM'

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash = SZCC_Cash + lens*((1./cosmopar%bias_act)-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias_act-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias_act-0.780)**2./(2.*0.092**2.)!final

         print*, 'SZCC_Cash final', SZCC_Cash

      !------------------------------------------------------------------------!

      case(3) ! Both Planck and ACT

         select case(dim_switch_plc)
         case(1) ! 1D Planck

            !Planck
            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            DO I=1,Nz             ! 10 bins for z which is delta z = 0.1
               sum=sum+DNz(I)     ! integration dN(z) over z
               ! print out predicted number counts
               !             open(unit=1, file = "d_numbers.txt")
               !             write(1, 200) I/10.0, DNz(I)
               !200          format(F5.1, F)
            ENDDO
            !          close(unit=1)

            if (print_counts==1) then
               print*,'prdicted counts of Planck 1D SZ'
               print*,DNz    !!!!!!!! N(z=1), N(z=2), N(z=3),...
               print*,'total counts',sum
            endif

            if (ieee_is_nan(DNZ(1))) then
               print*,'NaN found in theory counts!'
               stop
            endif

            DNzcat=DNzcat/dble(nred2)*dble(ncat)

            sum=0. ! Planck 1D likelihood
            do i=1,Nz
               if (DNz(i) /= 0.) then
                  ln_factorial=0.
                  if (DNzcat(i) /= 0.) ln_factorial=0.918939+(DNzcat(i)+0.5)*dlog(DNzcat(i))-DNzcat(i) !Stirling
                  sum=sum-1.*(DNzcat(i)*dlog(DNz(i))-DNz(i)-ln_factorial)
               end if
            end do
            SZCC_Cash_p = sum

            DNzcat = DNzcat*dble(nred2)/dble(ncat)

         case(2)

            ! Planck
            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            if (print_counts==1) print*,'predicted counts of Planck 2D SZ '
            DO I=1,Nz                                      !! 10 bins for z : delta z = 0.1
               if (print_counts==1) print*,I,DN(I,:)
               do J=1,Ny+1                                  !!  5 bins for q : delta log q = 0.25
                  sum=sum+DN(I,J)
               enddo
            ENDDO
            if (print_counts==1) print*,'total counts of Planck 2D SZ', sum

            sum=0. !!!!! Planck 2D likelihood
            do i=1,Nz
               do j=1,Ny+1
                  if (DN(i,j) /= 0.) then
                     ln_factorial=0.
                     if (DNcat(i,j)/=0.) then
                        if (DNcat(i,j) >10.) then
                           ln_factorial=0.918939+(DNcat(i,j)+0.5)*dlog(DNcat(i,j))-DNcat(i,j)
                           !Stirling approximation only for more than 10 elements
                        else
                           !direct computation of factorial
                           fact=1.
                           do ii=1,int(DNcat(i,j))
                              fact=ii*fact
                           enddo
                           ln_factorial=dlog(fact)
                        endif
                     endif
                     sum=sum-1.*(DNcat(i,j)*dlog(DN(i,j))-DN(i,j)-ln_factorial)
                  end if
               enddo
            end do
            SZCC_Cash_p = sum

         end select

         print*,'SZCC_Cash_Planck', sum
         !Print*,'SZ lnlike = ',SZCC_Cash_p

         !priors for SZ nuisance params
         SZCC_Cash_p = SZCC_Cash_p + pystar*(cosmopar%logystar-(-0.186))**2./(2.*0.021**2.) + &
         palpha*(cosmopar%alpha-1.789)**2./(2.*0.084**2.) + psigma*(cosmopar%sigmaM-0.075)**2./(2.*0.01**2.) !cosmomc_sz
         !print*,'prior on ystar, alpha, scatter'

         SZCC_Cash_p = SZCC_Cash_p +pbeta*(cosmopar%beta-0.6666666)**2/(2.*0.5**2.) !prior beta evolution

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash_p = SZCC_Cash_p + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash_p = SZCC_Cash_p + lens*(cosmopar%biasinv-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash_p = SZCC_Cash_p + wtg*(cosmopar%bias-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash_p = SZCC_Cash_p + cccp*(cosmopar%bias-0.780)**2./(2.*0.092**2.)!final


         select case(dim_switch_act)
         case(1)

            if (ylim_switch_act == 1) then  ! constant ylim

               call deltaN(Nz_a, z_a, delN) !!!

            elseif (ylim_switch_act == 2) then ! ylim from the map

               call Mlim_act_allocate(y0lim_arr)
               call y0_act_allocate(y0_act_arr)

!!!               call deltaN_act_ymap(Nz_a, z_a, delN, y0lim_arr) !
               call deltaN_act_ymap(Nz, z, delN, noise_act_arr, surveydeg2_act)

               !call deltaN_mass(Nz, z, delN, SZcat) !! to calculate ACT mass estimates!!!! temp!!!

               call Mlim_act_deallocate
               call y0_act_deallocate

            endif

            sum=0.
            do i = 1, Nz_a
               sum = sum + delN(i)
               ! print out predicted number counts
               !             open(unit=1, file = "d_redshift_number_act.txt")
               !             write(1, 100) i/10.0, delN(i)
               !100          format(F5.1, F)
            end do
            !          close(unit=1)

            if (print_counts==1) then
               print*, 'prdicted counts of ACT SZ 1D'
               print*, delN
               print*, 'total counts = ', sum
            endif

            if (ieee_is_nan(delN(1))) then
               print*,'NaN found in theory counts!'
               stop
            endif

            delNcat=delNcat/dble(nred2_a)*dble(ncat_a) !!!!! from catalogue

            sum=0.
            ! ACT 1D likelihood
            do i=1,Nz_a
               if (delN(i) /= 0.) then
                  ln_factorial=0.
                  if (delNcat(i) /= 0.) ln_factorial=0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
                  sum=sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
               end if
            end do
            SZCC_Cash_a = sum

            case(2)

               ! ACT 2D
               call Mlim_act_deallocate
               call deltaN2D(Nz_a, z_a, Ny0, y0, delN2D)
               call Mlim_act_deallocate

               sum=0.
               if (print_counts == 1) print*, 'predicted counts'
               do i = 1, Nz_a
                  if (print_counts == 1) print*, i, delN2D(i,:)
                  do j = 1, Ny0
                     sum = sum + delN2D(i,j)
                  end do
                  ! print out predicted number counts *****
                  !                 open(unit=1, file = "ACT_theoryN_2D.txt")
                  !                 write(1, 100) i/10.0, delN2D(i,1), delN2D(i,2), delN2D(i,3), delN2D(i,4), delN2D(i,5)
                  !100              format(F5.1, F10.5, F10.5, F10.5, F10.5, F10.5)
               end do
               !              close(unit=1)
               if (print_counts==1) print*, ' in each y0 bin '
               do j = 1, Ny0
                  sum = 0.
                  do i = 1, Nz_a
                     sum = sum + delN2D(i,j)
                  end do
                  print*, j, sum
               end do

               sum = 0.
               do j = 1, Ny0
                  do i = 1, Nz_a
                     sum = sum + delN2D(i,j)
                  end do
               end do
               print*, 'total counts =', sum

               sum=0.
               ! ACT 2D likelihood
               do i = 1, Nz_a
                  do j = 1, Ny0
                     if (delN2D(i,j) /= 0.) then
                        ln_factorial=0.
                        if (delN2Dcat(i,j) /= 0.) then
                           if (delN2Dcat(i,j) > 10.) then
                              ln_factorial = 0.918939 + (delN2Dcat(i,j) + 0.5)*dlog(delN2Dcat(i,j))-delN2Dcat(i,j) !Stirling approximation only for more than 10 elements
                           else
                              !direct computation of factorial
                              fact = 1.
                              do ii = 1, int(delN2Dcat(i,j))
                                 fact = ii*fact
                              end do
                              ln_factorial = dlog(fact)
                           end if
                        end if
                        sum = sum - 1.*(delN2Dcat(i,j)*dlog(delN2D(i,j)) - delN2D(i,j) - ln_factorial)
                     end if
                  end do
               end do
               SZCC_Cash_a = sum
            end select

            print*,'SZCC_Cash_ACT', sum
            !Print*,'SZ lnlike = ', SZCC_Cash_A

            !PRIORS ON NS AND OMEGABH2:
            SZCC_Cash_a = SZCC_Cash_a + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
            oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

            !PRIORS ON THE BIAS:
            SZCC_Cash_a = SZCC_Cash_a + lens*((1./cosmopar%bias_act)-0.99)**2./(2.*0.19**2.) !final
            !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

            !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
            !1-b=0.81+-0.08. !clash
            SZCC_Cash_a = SZCC_Cash_a + wtg*(cosmopar%bias_act-0.688)**2./(2.*0.072**2.)!final

            SZCC_Cash_a = SZCC_Cash_a + cccp*(cosmopar%bias_act-0.780)**2./(2.*0.092**2.)!final

            SZCC_Cash = SZCC_Cash_p + SZCC_Cash_a

            print*, 'SZCC_Cash Planck + ACT', SZCC_Cash

      !------------------------------------------------------------------------!

      case(4) ! SPT

         select case(dim_switch_spt)
         case(1)

             call deltaN_spt(z, Nz, logxi, Nxi, delN, delN2D, dim_switch_spt)

             sum=0.
             do i = 1, Nz
                sum = sum + delN(i)
                ! print out predicted number counts
                !          open(unit=1, file = "spt_dNdz_bias.txt")
                !          write(1, 100) i/10.0 - 0.1, delN(i)
                !100       format(F5.1, F)
             end do
             !       close(unit=1)

             if (print_counts==1) then
                print*, 'prdicted counts'
                print*, delN
                print*, 'total counts = ', sum
                print*, 'signal to noise of SPT =', snr_spt_cat
             endif

             if (ieee_is_nan(delN(1))) then
                print*,'NaN found in theory counts!'
                stop
             endif

             delNcat=delNcat/dble(nred2)*dble(ncat)
             delNcat=delNcat/365.*377.
             !print*, 'delNcat = ',  delNcat

             !sum=0.
             !DO I=1,Nz
             !   sum=sum+delNcat(I)
             !ENDDO
             !print*, 'Total = ', sum

             sum = 0.
             do i = 1, Nz
                if (delN(i) /= 0.) then
                   ln_factorial = 0.
                   if (delNcat(i) /= 0.) ln_factorial = 0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
                   sum = sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
                end if
             end do
             SZCC_Cash = sum

             delNcat = delNcat*dble(nred2)/dble(ncat)
             delNcat = delNcat*365./377.

         case(2)

            call deltaN_spt(z, Nz, logxi, Nxi, delN, delN2D, dim_switch_spt)

            sum = 0.
            if (print_counts == 1) print*, 'predicted counts'
            do i = 1, Nz
                if (print_counts == 1) print*, i, delN2D(i,:)
                do j = 1, Nxi+1
                    sum = sum + delN2D(i,j)
                enddo
            enddo

            if (print_counts == 1) print*, 'total counts', sum

            do j = 1, Nxi+1
               sum = 0.
               do i = 1, Nz
                  sum = sum + delN2D(i,j)
               end do
               print*, j, sum
            end do




            delNcat=delNcat/dble(nred2)*dble(ncat)
            delN2Dcat = delN2Dcat/365.*377.

            sum = 0.
            do i = 1, Nz
                do j = 1, Nxi+1
                    if (delN2D(i,j) /= 0.) then
                        ln_factorial = 0.
                        if (delN2Dcat(i,j) /= 0.) then
                            if (delN2Dcat(i,j) > 10.) then
                                ln_factorial = 0.918939 + (delN2Dcat(i,j)+0.5)*dlog(delN2Dcat(i,j)) - delN2Dcat(i,j)
                                !Stirling approximation only for more than 10 elements
                            else
                                !direct computation of factorial
                                fact = 1.
                                do ii = 1, int(delN2Dcat(i,j))
                                    fact = ii*fact
                                enddo
                                ln_factorial = dlog(fact)
                            endif
                        endif
                        sum = sum-1.*(delN2Dcat(i,j)*dlog(delN2D(i,j)) - delN2D(i,j) - ln_factorial)
                    endif
                enddo
            enddo
            SZCC_Cash = sum

            delNcat = delNcat*dble(nred2)/dble(ncat)
            delN2Dcat = delN2Dcat*365./377.

         end select



         print*,'SZCC_Cash', sum
         !Print*,'SZ lnlike = ',SZCC_Cash

         !priors on nusiance parameters:

!         SZCC_Cash = SZCC_Cash + pA_s*(cosmopar%a_spt - 5.38 )**2./(2.*1.61 **2.) + &
!                        pB_s*(cosmopar%b_spt - 1.340)**2./(2.*0.268**2.) + &
!                        pC_s*(cosmopar%c_spt - 0.49 )**2./(2.*0.49 **2.) + &
!                        pD_s*(cosmopar%d_spt - 0.13 )**2./(2.*0.13 **2.)

         SZCC_Cash = SZCC_Cash + pA_s*(cosmopar%a_spt - 4.842)**2./(2.*0.913**2.) + &
                        pB_s*(cosmopar%b_spt - 1.668)**2./(2.*0.083**2.) + &
                        pC_s*(cosmopar%c_spt - 0.550)**2./(2.*0.315**2.) + &
                        pD_s*(cosmopar%d_spt - 0.199)**2./(2.*0.069**2.)


         !print*, 'prior on A_s, B_s, C_s, D_s for SPT'

         !PRIORS ON NS AND OMEGABH2:
         !SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         !     oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         SZCC_Cash = SZCC_Cash + & ! pns*(CMB%InitPower(ns_index)-0.9619)**2./(2.*0.0073 **2.) +
                        oh2*(CMB%ombh2 - 0.02202)**2/(2*0.00045**2.) + &
                        pH0*(CMB%H0 - 73.8)**2./(2*2.4**2.) ! n_s and BBN : comparison with deHaan

         !print*, 'prior on n_s, omegabh2, H0 for SPT'

         !PRIORS ON THE BIAS:
         SZCC_Cash = SZCC_Cash + lens*((1./cosmopar%bias_spt)-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias_spt-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias_spt-0.780)**2./(2.*0.092**2.)!final

         !print*, 'SZCC_Cash final = ', SZCC_Cash



      !------------------------------------------------------------------------!

      case(5) ! Both Planck and SPT

         select case(dim_switch_plc)
         case(1) ! 1D Planck

            !Planck
            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            DO I=1,Nz             ! 10 bins for z which is delta z = 0.1
               sum=sum+DNz(I)     ! integration dN(z) over z
               ! print out predicted number counts
               !             open(unit=1, file = "d_numbers.txt")
               !             write(1, 200) I/10.0, DNz(I)
               !200          format(F5.1, F)
            ENDDO
            !          close(unit=1)

            if (print_counts==1) then
               print*,'prdicted counts of Planck 1D SZ'
               print*,DNz    !!!!!!!! N(z=1), N(z=2), N(z=3),...
               print*,'total counts',sum
            endif

            if (ieee_is_nan(DNZ(1))) then
               print*,'NaN found in theory counts!'
               stop
            endif

            DNzcat=DNzcat/dble(nred2)*dble(ncat)

            sum=0. ! Planck 1D likelihood
            do i=1,Nz
               if (DNz(i) /= 0.) then
                  ln_factorial=0.
                  if (DNzcat(i) /= 0.) ln_factorial=0.918939+(DNzcat(i)+0.5)*dlog(DNzcat(i))-DNzcat(i) !Stirling
                  sum=sum-1.*(DNzcat(i)*dlog(DNz(i))-DNz(i)-ln_factorial)
               end if
            end do
            SZCC_Cash_p = sum

            DNzcat = DNzcat*dble(nred2)/dble(ncat)

         case(2)

            ! Planck
            call deltaN_yz(Z,Nz,LOGY,Ny,DNZ,DN,skyfracs,thetas,ylims,dim_switch_plc,qa_ytot,erf_list)

            sum=0.
            if (print_counts==1) print*,'predicted counts of Planck 2D SZ '
            DO I=1,Nz                                      !! 10 bins for z : delta z = 0.1
               if (print_counts==1) print*,I,DN(I,:)
               do J=1,Ny+1                                  !!  5 bins for q : delta log q = 0.25
                  sum=sum+DN(I,J)
               enddo
            ENDDO
            if (print_counts==1) print*,'total counts of Planck 2D SZ', sum

            sum=0. !!!!! Planck 2D likelihood
            do i=1,Nz
               do j=1,Ny+1
                  if (DN(i,j) /= 0.) then
                     ln_factorial=0.
                     if (DNcat(i,j)/=0.) then
                        if (DNcat(i,j) >10.) then
                           ln_factorial=0.918939+(DNcat(i,j)+0.5)*dlog(DNcat(i,j))-DNcat(i,j)
                           !Stirling approximation only for more than 10 elements
                        else
                           !direct computation of factorial
                           fact=1.
                           do ii=1,int(DNcat(i,j))
                              fact=ii*fact
                           enddo
                           ln_factorial=dlog(fact)
                        endif
                     endif
                     sum=sum-1.*(DNcat(i,j)*dlog(DN(i,j))-DN(i,j)-ln_factorial)
                  end if
               enddo
            end do
            SZCC_Cash_p = sum

         end select

         print*,'SZCC_Cash_Planck', sum
         !Print*,'SZ lnlike = ',SZCC_Cash_p

         !priors for SZ nuisance params
         SZCC_Cash_p = SZCC_Cash_p + pystar*(cosmopar%logystar-(-0.186))**2./(2.*0.021**2.) + &
         palpha*(cosmopar%alpha-1.789)**2./(2.*0.084**2.) + psigma*(cosmopar%sigmaM-0.075)**2./(2.*0.01**2.) !cosmomc_sz
         !print*,'prior on ystar, alpha, scatter'

         SZCC_Cash_p = SZCC_Cash_p +pbeta*(cosmopar%beta-0.6666666)**2/(2.*0.5**2.) !prior beta evolution

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash_p = SZCC_Cash_p + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash_p = SZCC_Cash_p + lens*(cosmopar%biasinv-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash_p = SZCC_Cash_p + wtg*(cosmopar%bias-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash_p = SZCC_Cash_p + cccp*(cosmopar%bias-0.780)**2./(2.*0.092**2.)!final

         ! SPT

         select case(dim_switch_spt)
         case(1)

             call deltaN_spt(z_s, Nz_s, logxi, Nxi, delN, delN2D, dim_switch_spt)

             sum=0.
             do i = 1, Nz_s
                sum = sum + delN(i)
                ! print out predicted number counts
                !          open(unit=1, file = "spt_dNdz_bias.txt")
                !          write(1, 100) i/10.0 - 0.1, delN(i)
                !100       format(F5.1, F)
             end do
             !       close(unit=1)

             if (print_counts==1) then
                print*, 'prdicted counts'
                print*, delN
                print*, 'total counts = ', sum
                print*, 'signal to noise of SPT =', snr_spt_cat
             endif

             if (ieee_is_nan(delN(1))) then
                print*,'NaN found in theory counts!'
                stop
             endif

             delNcat=delNcat/dble(nred2_s)*dble(ncat_s)
             delNcat=delNcat/365.*377.
             !print*, 'delNcat = ',  delNcat

             !sum=0.
             !DO I=1,Nz
             !   sum=sum+delNcat(I)
             !ENDDO
             !print*, 'Total = ', sum

             sum = 0.
             do i = 1, Nz_s
                if (delN(i) /= 0.) then
                   ln_factorial = 0.
                   if (delNcat(i) /= 0.) ln_factorial = 0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
                   sum = sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
                end if
             end do
             SZCC_Cash_s = sum

             delNcat = delNcat*dble(nred2_s)/dble(ncat_s)
             delNcat = delNcat*365./377.

         case(2)

            call deltaN_spt(z_s, Nz_s, logxi, Nxi, delN, delN2D, dim_switch_spt)

            sum = 0.
            if (print_counts == 1) print*, 'predicted counts'
            do i = 1, Nz_s
                if (print_counts == 1) print*, i, delN2D(i,:)
                do j = 1, Nxi+1
                    sum = sum + delN2D(i,j)
                enddo
            enddo

            if (print_counts == 1) print*, 'total counts', sum

            delNcat=delNcat/dble(nred2_s)*dble(ncat_s)
            delN2Dcat = delN2Dcat/365.*377.

            sum = 0.
            do i = 1, Nz_s
                do j = 1, Nxi+1
                    if (delN2D(i,j) /= 0.) then
                        ln_factorial = 0.
                        if (delN2Dcat(i,j) /= 0.) then
                            if (delN2Dcat(i,j) > 10.) then
                                ln_factorial = 0.918939 + (delN2Dcat(i,j)+0.5)*dlog(delN2Dcat(i,j)) - delN2Dcat(i,j)
                                !Stirling approximation only for more than 10 elements
                            else
                                !direct computation of factorial
                                fact = 1.
                                do ii = 1, int(delN2Dcat(i,j))
                                    fact = ii*fact
                                enddo
                                ln_factorial = dlog(fact)
                            endif
                        endif
                        sum = sum-1.*(delN2Dcat(i,j)*dlog(delN2D(i,j)) - delN2D(i,j) - ln_factorial)
                    endif
                enddo
            enddo
            SZCC_Cash_s = sum

            delNcat = delNcat*dble(nred2_s)/dble(ncat_s)
            delN2Dcat = delN2Dcat*365./377.

         end select

         print*,'SZCC_Cash_SPT', sum
         !Print*,'SZ lnlike = ', SZCC_Cash_A


         SZCC_Cash = SZCC_Cash + pA_s*(cosmopar%a_spt - 4.842)**2./(2.*0.913**2.) + &
         pB_s*(cosmopar%b_spt - 1.668)**2./(2.*0.083**2.) + &
         pC_s*(cosmopar%c_spt - 0.550)**2./(2.*0.315**2.) + &
         pD_s*(cosmopar%d_spt - 0.199)**2./(2.*0.069**2.)

         !print*, 'prior on A_s, B_s, C_s, D_s for SPT'

         !PRIORS ON NS AND OMEGABH2:
         !SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         !     oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         SZCC_Cash = SZCC_Cash + & ! pns*(CMB%InitPower(ns_index)-0.9619)**2./(2.*0.0073 **2.) +
                  oh2*(CMB%ombh2 - 0.02202)**2/(2*0.00045**2.) + &
                  pH0*(CMB%H0 - 73.8)**2./(2*2.4**2.) ! n_s and BBN : comparison with deHaan

         !print*, 'prior on n_s, omegabh2, H0 for SPT'

         !PRIORS ON THE BIAS:
         SZCC_Cash_s = SZCC_Cash_s + lens*((1./cosmopar%bias_spt)-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash_s = SZCC_Cash_s + wtg*(cosmopar%bias_spt-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash_s = SZCC_Cash_s + cccp*(cosmopar%bias_spt-0.780)**2./(2.*0.092**2.)!final

         SZCC_Cash = SZCC_Cash_p + SZCC_Cash_s

         print*, 'SZCC_Cash Planck + SPT', SZCC_Cash

      !------------------------------------------------------------------------!

      case(6) ! SO single tile test

!         call y0_SO_allocate(y0_SO_arr)
!         call deltaN_SO(Nz, z, delN, noise_SOtest, surveydeg2_SO)
!         call y0_SO_deallocate

         call deltaN_SOtest(z, Nz, delN, skyfracs_SOtest, noise_SOtest, theta, Q0)


         !          sum=0.
         !          do i = 1, Nz
         !           	sum = sum + delN(i)
         !          end do
         !          print*, 'number counts', sum

         sum=0.
         do i = 1, Nz
            sum = sum + delN(i)
            !print out predicted number counts
            !           open(unit=1, file = "act_dNdz_B0.txt")
            !           write(1, 113) i/10.0 - 0.1, delN(i)
            !113        format(F5.1, F)
         end do
         !      close(unit=1)

         if (print_counts==1) then
            print*, 'prdicted counts'
            print*, delN
            print*, 'total counts = ', sum
         endif

         if (ieee_is_nan(delN(1))) then
            print*,'NaN found in theory counts!'
            stop
         endif

         sum=0. ! SO 1D likelihood
         do i=1,Nz
            if (delN(i) /= 0.) then
               ln_factorial=0.
               if (delNcat(i) /= 0.) ln_factorial=0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
               sum=sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
            end if
         end do
         SZCC_Cash=sum
         print*,'SZCC_Cash SO', sum
         !Print*,'SZ lnlike = ',SZCC_Cash

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash = SZCC_Cash + lens*((1./cosmopar%bias_act)-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias_act-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias_act-0.780)**2./(2.*0.092**2.)!final

         !print*, 'SZCC_Cash final', SZCC_Cash

      !------------------------------------------------------------------------!

      case(7) ! SO full map with 540 tiles - now with 346 tiles

!         call y0_SO_tiles_allocate(y0_SO_tiles_arr) !
!         call deltaN_SO_tiles(Nz, z, delN, noise_SOfull, skyfracs_SOfull, tilenames) !
!         call y0_SO_tiles_deallocate !

         !call precalc_y0_allocate(z, Nz, skyfracs_SOfull, tilenames, precalc_y0)
         call deltaN_SO_fullmap(z, Nz, delN, skyfracs_SOfull, noise_SOfull, tilenames, theta, Q00)
         !call precalc_y0_deallocate

         !          sum=0.
         !          do i = 1, Nz
         !           	sum = sum + delN(i)
         !          end do
         !          print*, 'number counts', sum

         sum=0.
         do i = 1, Nz
            sum = sum + delN(i)
            !print out predicted number counts
            !           open(unit=1, file = "act_dNdz_B0.txt")
            !           write(1, 113) i/10.0 - 0.1, delN(i)
            !113        format(F5.1, F)
         end do
         !      close(unit=1)

         if (print_counts==1) then
            print*, 'prdicted counts'
            print*, delN
            print*, 'total counts = ', sum
         endif

         if (ieee_is_nan(delN(1))) then
            print*,'NaN found in theory counts!'
            stop
         endif

         sum=0. ! SO 1D likelihood
         do i=1,Nz
            if (delN(i) /= 0.) then
               ln_factorial=0.
               if (delNcat(i) /= 0.) ln_factorial=0.918939+(delNcat(i)+0.5)*dlog(delNcat(i))-delNcat(i) !Stirling
               sum=sum-1.*(delNcat(i)*dlog(delN(i))-delN(i)-ln_factorial)
            end if
         end do
         SZCC_Cash=sum
         print*,'SZCC_Cash SO', sum
         !Print*,'SZ lnlike = ',SZCC_Cash

         !PRIORS ON NS AND OMEGABH2:
         SZCC_Cash = SZCC_Cash + pns*(CMB%InitPower(ns_index)-0.9624)**2./(2.*0.014 **2.) + &
         oh2*(CMB%ombh2 - 0.022)**2/(2*0.002**2.)  !BBN

         !PRIORS ON THE BIAS:
         SZCC_Cash = SZCC_Cash + lens*((1./cosmopar%bias_act)-0.99)**2./(2.*0.19**2.) !final
         !(1-b)^(-1)=0.99 pm 0.19. !Planck CMB lensing

         !SZCC_Cash = SZCC_Cash + clash*(cosmopar%bias-0.8)**2/(2.*0.11**2)!not final
         !1-b=0.81+-0.08. !clash
         SZCC_Cash = SZCC_Cash + wtg*(cosmopar%bias_act-0.688)**2./(2.*0.072**2.)!final

         SZCC_Cash = SZCC_Cash + cccp*(cosmopar%bias_act-0.780)**2./(2.*0.092**2.)!final

         !print*, 'SZCC_Cash final', SZCC_Cash

      end select

   end function SZCC_Cash
   !---------------------------------------------------------------------------!

end module szcounts
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

