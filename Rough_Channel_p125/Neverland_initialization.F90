module Neverland_initialization
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

public Neverland_initialize_topography
public Neverland_initialize_thickness

contains

! -----------------------------------------------------------------------------
!> This subroutine sets up the Neverland test case topography.
subroutine Neverland_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the Neverland test case topography
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: x, y
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "Neverland_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  real :: nl_roughness_amp
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_topography: setting topography", 5)

  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "NL_ROUGHNESS_AMP", nl_roughness_amp, &
                 "Amplitude of wavy signal in bathymetry.", default=0.05)

  PI = 4.0*atan(1.0)

!  Calculate the depth of the bottom.
  do i=is,ie
  do j=js,je
    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat
!  This sets topography that has a reentrant channel to the south.

    if ((j.eq.1) .or. (j.eq.jed)) then 
        D(i,j) = 0.0
    else  
        D(i,j) = 1.0  -  1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12)      !< The great northern wall and Antarctica   
                      -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y) &   !< roughness
                      -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)     !< roughness
   endif

!                - (1.2 * spike(x,0.2) + 1.2 * spike(x-1.0,0.2)) * spike(MIN(0.0,y-0.3),0.2) & !< South America
!                -  1.2 * spike(x-0.5,0.2) * spike(MIN(0.0,y-0.55),0.2)                       & !< Africa
!                -  1.1 * spike(y-1,0.12) - 1.1 * spike(y,0.12)                               & !< The great northern wall and Antarctica
!                -  1.2 * (spike(x,0.12)  + spike(x-1,0.12)) * spike(MAX(0.0,y-0.06),0.12)    & !< Antarctic Peninsula
!                -  0.1 * (cosbell(x,0.1) + cosbell(x-1,0.1))                                   & !< Drake Passage ridge
!                -  0.5 * cosbell(x-0.16,0.05) * (cosbell(y-0.18,0.13)**0.4)                  & !< Scotia Arc East
!                -  0.4 * (cosbell(x-0.09,0.08)**0.4) * cosbell(y-0.26,0.05)                  & !< Scotia Arc North
!                -  0.4 * (cosbell(x-0.08,0.08)**0.4) * cosbell(y-0.1,0.05)                   & !< Scotia Arc South
!                -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y)                            & !< roughness
!                -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)                              !< roughness


    if (D(i,j) < 0.0) D(i,j) = 0.0
    D(i,j) = D(i,j) * max_depth
  enddo
  enddo



end subroutine Neverland_initialize_topography
! -----------------------------------------------------------------------------

!> Returns the value of a cosine-bell function evaluated at x/L
 real function cosbell(x,L)

   real , intent(in) :: x       !< non-dimensional position
   real , intent(in) :: L       !< non-dimensional width
   real              :: PI      !< 3.1415926... calculated as 4*atan(1)

   PI      = 4.0*atan(1.0)
   cosbell = 0.5 * (1 + cos(PI*MIN(ABS(x/L),1.0)))
 end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L
 real function spike(x,L)

   real , intent(in) :: x       !< non-dimensional position
   real , intent(in) :: L       !< non-dimensional width
   real              :: PI      !< 3.1415926... calculated as 4*atan(1)

   PI    = 4.0*atan(1.0)
   spike = (1 - sin(PI*MIN(ABS(x/L),0.5)))
 end function spike


! -----------------------------------------------------------------------------
!> This subroutine initializes layer thicknesses for the Neverland test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverland_initialize_thickness(h, G, GV, param_file, eqn_of_state, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  real, intent(out), dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< The thickness that is being
                                                              !! initialized.
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  type(EOS_type),          pointer    :: eqn_of_state         !< integer that selects the
                                                              !! equation of state.
  real,                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure in Pa.
  ! Local variables
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                            ! negative because it is positive upward.      !
  real, dimension(SZK_(G)) :: h_profile ! Vector of initial thickness profile (m)
  real :: e_interface ! Current interface positoin (m)
  real :: phase
  real :: x,y
  real*8 :: RANDOM
  character(len=40)  :: mod = "Neverland_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, k1, is, ie, js, je, nz, itt, isd, ied, jsd, jed, junk, date_time(8)
  character*10 :: b(3)
   
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverland_initialization.F90, Neverland_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mod, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", fail_if_missing=.true.)

! e0 is the notional position of interfaces
  e0(1) = 0. ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie 
    x=(G%geoLonT(i,j)-G%west_lon)/G%len_lon
    y=(G%geoLatT(i,j)-G%south_lat)/G%len_lat 
    e_interface = -G%bathyT(i,j)
    call date_and_time(b(1), b(2), b(3), date_time)
    RANDOM = 100000 * (x * y)
!    print*, '            '
!    print*, '            ' 
!    print*,'RANDOM=',RANDOM
!    print*, '            '
!    print*, '            '
    junk = int(RANDOM * i * j) 
    do k=nz,1,-1
      phase = ranno(junk,0) 
!      print*, 'junk=',junk
!      print*,'phase=',phase
!      print*, '            '
!      print*, '            '
      h(i,j,k) = max( GV%Angstrom_z, e0(k) - e_interface ) 
      h(i,j,k) = h(i,j,k) + 0.0005*h(i,j,k)*(phase)
      e_interface = max( e0(k), e_interface - h(i,j,k) )
    enddo

  enddo ; enddo

end subroutine Neverland_initialize_thickness

!> Returns the value of a random number RANNO
real function ranno (junk,i)
! Controls random number generator.
!-----------------------------------------
! - If argument i.ne.0 it performs initialization with i=seed no.
! - If argument i.eq.0 it draws a random no.
!-----------------------------------------
  integer :: i,junk,ihold
  real :: twopi
!  save junk

  twopi =  4.*asin(1.)  
  if (i.ne.0) then
     if (i.gt.0) i = - i
     junk  = i
     ranno = (ran1(i)-0.5)*twopi
  else
     junk  = junk - 1
     ihold = junk
!     print*, '            '
!     print*, '            '
!     print*, 'JUNK2=',junk
!     print*, '            '
!     print*, '            '
     ranno = (ran1(ihold)-0.5)*twopi
  endif
  return
end function ranno

real function ran1(idum)
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ntab=32,ndiv=1+(im-1)/ntab
  real, parameter :: am=1./im,eps=1.2e-7,rnmx=1.-eps
  integer :: idum,j,k,iv(32),iy,junk
  save iv,iy
  data iv /ntab*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do j=ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        if (j.le.ntab) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum.lt.0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
end function ran1

! -----------------------------------------------------------------------------

!! \class Neverland_initialization
!!
!! The module configures the model for the Neverland experiment.
end module Neverland_initialization
