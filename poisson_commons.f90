module poisson_commons
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons

!#############################################
!QUMOND module
!Contains routines to initialise extra arrays
!and to add the phantom dark matter to the
!density. Also switches between Newtonian
!and MONDian quantities as appropriate for
!the solver.
!#############################################

module mond_commons
   use amr_parameters
   implicit none

   integer :: mond_n = 2         !Read in from namelist
   real(dp) :: imond_a0 = 1.2d-20 !Read in from namelist
   ! mond_cosmo is only for cosmo runs (obviously) -
   ! used to adjust MOND scale over cosmological time according
   ! to simple prescription: a0 = aexp**mond_cosmo*a0
   ! Note that a0 already has a factor of aexp from moving to
   ! (super)co-moving coordinates.
   integer :: mond_cosmo = 1     !Read in from namelist
   ! mond_omega_m is only for cosmo runs - sets particle
   ! mass to chosen value without changing background evolution
   real(dp) :: mond_omega_m = 0.3 !Read in from namelist

   integer :: tnbors = twondim*ndim

   integer, dimension(1:3,1:6,1:8) :: kkk, lll, mmm
   integer, dimension(1:6,1:8) :: uuu, yyy
   real(dp) :: mond_a0

   logical :: qumond_pot

   real(dp), allocatable, dimension(:) :: rho_pdm
   real(dp), allocatable, dimension(:) :: newtonian_phi
   real(dp), allocatable, dimension(:) :: mondian_phi
   integer :: rho_output_count=0,phi_output_count=0

contains   

subroutine init_mond_arrays
   use amr_parameters
   use amr_commons, only : myid
   implicit none

   real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

   kkk(1,1,1:8) = (/1,0,1,0,1,0,1,0/) ; lll(1,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(1,2,1:8) = (/0,2,0,2,0,2,0,2/) ; lll(1,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(1,3,1:8) = (/1,3,1,0,1,3,1,0/) ; lll(1,3,1:8) = (/3,0,0,0,3,0,0,0/)
   kkk(1,4,1:8) = (/0,2,4,2,0,2,4,2/) ; lll(1,4,1:8) = (/0,0,0,4,0,0,0,4/)
   kkk(1,5,1:8) = (/1,5,1,5,1,0,1,0/) ; lll(1,5,1:8) = (/5,0,5,0,0,0,0,0/)
   kkk(1,6,1:8) = (/0,2,0,2,6,2,6,2/) ; lll(1,6,1:8) = (/0,0,0,0,0,6,0,6/)

   kkk(2,1,1:8) = (/3,3,0,0,3,3,0,0/) ; lll(2,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(2,2,1:8) = (/0,0,4,4,0,0,4,4/) ; lll(2,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(2,3,1:8) = (/3,2,0,2,3,2,0,2/) ; lll(2,3,1:8) = (/0,3,0,0,0,3,0,0/)
   kkk(2,4,1:8) = (/1,0,1,4,1,0,1,4/) ; lll(2,4,1:8) = (/0,0,4,0,0,0,4,0/)
   kkk(2,5,1:8) = (/3,3,5,5,3,3,0,0/) ; lll(2,5,1:8) = (/5,5,0,0,0,0,0,0/)
   kkk(2,6,1:8) = (/0,0,4,4,6,6,4,4/) ; lll(2,6,1:8) = (/0,0,0,0,0,0,6,6/)

   kkk(3,1,1:8) = (/5,5,5,5,0,0,0,0/) ; lll(3,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(3,2,1:8) = (/0,0,0,0,6,6,6,6/) ; lll(3,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(3,3,1:8) = (/5,2,5,2,0,2,0,2/) ; lll(3,3,1:8) = (/0,5,0,5,0,0,0,0/)
   kkk(3,4,1:8) = (/1,0,1,0,1,6,1,6/) ; lll(3,4,1:8) = (/0,0,0,0,6,0,6,0/)
   kkk(3,5,1:8) = (/5,5,5,5,0,0,4,4/) ; lll(3,5,1:8) = (/0,0,4,4,0,0,0,0/)
   kkk(3,6,1:8) = (/3,3,0,0,6,6,6,6/) ; lll(3,6,1:8) = (/0,0,0,0,3,3,0,0/)

   mmm(1,1,1:8) = (/2,1,4,3,6,5,8,7/) ; uuu(1,1:8) = (/1,0,1,0,1,0,1,0/)
   mmm(1,2,1:8) = (/2,1,4,3,6,5,8,7/) ; uuu(2,1:8) = (/0,1,0,1,0,1,0,1/)
   mmm(1,3,1:8) = (/4,3,2,1,8,7,6,5/) ; uuu(3,1:8) = (/1,1,0,0,1,1,0,0/)
   mmm(1,4,1:8) = (/4,3,2,1,8,7,6,5/) ; uuu(4,1:8) = (/0,0,1,1,0,0,1,1/)
   mmm(1,5,1:8) = (/6,5,8,7,2,1,4,3/) ; uuu(5,1:8) = (/1,1,1,1,0,0,0,0/)
   mmm(1,6,1:8) = (/6,5,8,7,2,1,4,3/) ; uuu(6,1:8) = (/0,0,0,0,1,1,1,1/)

   mmm(2,1,1:8) = (/3,4,1,2,7,8,5,6/) ; yyy(1,1:8) = (/2,1,4,3,6,5,8,7/)
   mmm(2,2,1:8) = (/3,4,1,2,7,8,5,6/) ; yyy(2,1:8) = (/2,1,4,3,6,5,8,7/)
   mmm(2,3,1:8) = (/4,3,2,1,8,7,6,5/) ; yyy(3,1:8) = (/3,4,1,2,7,8,5,6/)
   mmm(2,4,1:8) = (/4,3,2,1,8,7,6,5/) ; yyy(4,1:8) = (/3,4,1,2,7,8,5,6/)
   mmm(2,5,1:8) = (/7,8,5,6,3,4,1,2/) ; yyy(5,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(2,6,1:8) = (/7,8,5,6,3,4,1,2/) ; yyy(6,1:8) = (/5,6,7,8,1,2,3,4/)

   mmm(3,1,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(3,2,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(3,3,1:8) = (/6,5,8,7,2,1,4,3/)
   mmm(3,4,1:8) = (/6,5,8,7,2,1,4,3/)
   mmm(3,5,1:8) = (/7,8,5,6,3,4,1,2/)
   mmm(3,6,1:8) = (/7,8,5,6,3,4,1,2/)

   ! Only valid for a non-cosmo run!
   ! Because the scale factors change with time in a cosmo run.
   if (.not. cosmo) then
       call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
       mond_a0 = imond_a0*100.0d0*scale_t**2/scale_l
   end if

end subroutine init_mond_arrays

subroutine make_phantom_halo(ilevel,icount)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel,icount

   real(dp), dimension(1:tnbors+1) :: phi_nbor
   real(dp), dimension(1:twondim) :: phi_vals
   real(dp), dimension(1:twondim) :: nu

   ! Arrays for vectorized interpol_phi
   real(dp), dimension(1:nvector,1:twotondim) :: phi_int
   integer,  dimension(1:nvector) :: ind_cell

   real(dp) :: dx, dx2, nb_sum, factor, grad_sqrd, phi_sum
   real(dp) :: oneoverdx, oneoverfourdx, oneoverfourpidx2
   real(dp) :: scale, fourpi, nb_phi, umond_a0

   integer  :: ngrid, nx_loc, ind
   integer  :: igrid_mg, idim, inbor, nbors_ind
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: kgshift, lgshift, igrid_nbor_amr
   integer  :: cpu_nbor_amr, icell_nbor_amr, ifathercell_nbor_amr
   integer  :: upper_cell_index, upper_cell_amr, upper_grid_amr
   integer  :: nbor_grid_amr, nbor_cell_index, nbor_cell_amr
   integer  :: rhosize

   real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! Set constants
   dx2  = (0.5d0**ilevel)**2
   dx   = 0.5d0**ilevel
   oneoverdx = 1.0d0/dx
   oneoverfourdx = 1.0d0/(4.0d0*dx)

   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) fourpi = 1.5D0*omega_m*aexp*scale
   oneoverfourpidx2 = 1.0d0/(fourpi*dx2)

   if(cosmo) then
       umond_a0 = (imond_a0*100.0d0*aexp*scale_t**2/scale_l)*(aexp**mond_cosmo)
       ! Including one factor of aexp for co-moving coords, then the additional
       ! factor is for an additional cosmo dependence.
       !if ((myid)==1) print *,"MOND_A0: ",umond_a0
   else
       umond_a0 = mond_a0
   endif

   ngrid=active(ilevel)%ngrid

   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         phi_nbor(tnbors+1) = phi(icell_amr)

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbour gridshift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid
                  if(kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                     ifathercell_nbor_amr = father(igrid_nbor_amr)
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                     ifathercell_nbor_amr = nbor(igrid_amr,kgshift)
                  else
                     igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                     ifathercell_nbor_amr = nbor(son(nbor(igrid_amr,kgshift)),lgshift)
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No DIAGONAL neighbor
                     ! Interpolate
                     ind_cell(1)=ifathercell_nbor_amr
                     call interpol_phi(ind_cell,phi_int,1,ilevel,icount)
                     nb_phi = phi_int(1,mmm(idim,inbor,ind))
                     phi_nbor(nbors_ind) = nb_phi
                     !phi_nbor(nbors_ind) = 0.5d0*nb_phi + 0.5d0*phi(icell_amr)
                     !phi_nbor(nbors_ind) = phi(icell_amr)
                  else
                     ! Get phi values on neighbouring cells
                     icell_nbor_amr  = igrid_nbor_amr + &
                        (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind) = phi(icell_nbor_amr)
                  endif
               end do
            end do
         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks
            if (f(icell_amr,3)<=0.0) cycle
            if (safe_mode(ilevel) .and. f(icell_amr,3)<1.0) cycle
            ! Use more complex mu calculation valid at boundaries

            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbor grid shift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid and its parent cell
                  if (kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                     ifathercell_nbor_amr = father(igrid_nbor_amr)
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                     ifathercell_nbor_amr = nbor(igrid_amr,kgshift)
                  else
                     if (son(nbor(igrid_amr,kgshift)) == 0) then
                     ! special case when we can't "leapfrog" to the diagonal cell because
                     ! the first shift takes us out of the AMR mesh.
                     ! First, integer division to find which number of cell we have in upper level grid
                     ! and then we calculate the amr index of that upper level grid
                        upper_cell_amr = nbor(igrid_amr,kgshift)
                        upper_cell_index = (upper_cell_amr/ngridmax)+1
                        ! Check if we need to move to the neighbour grid at the upper level or not
                        if (uuu(lgshift,upper_cell_index)==0) then
                            nbor_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                        else
                            upper_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                            nbor_grid_amr = son(nbor(upper_grid_amr,lgshift))
                        endif
                        ! Determine cell index of neighbour cell, depending on direction of 2nd shift
                        nbor_cell_index = yyy(lgshift,upper_cell_index)
                        ! Calculate amr index of neighbour cell
                        nbor_cell_amr = nbor_grid_amr + ncoarse + (nbor_cell_index - 1)*ngridmax
                        ! Find son grid in AMR mesh, if it exists
                        igrid_nbor_amr = son(nbor_cell_amr)
                        ifathercell_nbor_amr = nbor_cell_amr
                     else
                        igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                        ifathercell_nbor_amr = nbor(son(nbor(igrid_amr,kgshift)),lgshift)
                     endif
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No neighbor
                     ! Use interpolated values from upper level
                     ind_cell(1)=ifathercell_nbor_amr
                     call interpol_phi(ind_cell,phi_int,1,ilevel,icount)
                     nb_phi = phi_int(1,mmm(idim,inbor,ind))
                     phi_nbor(nbors_ind) = nb_phi
                     !phi_nbor(nbors_ind) = 0.5d0*nb_phi + 0.5d0*phi(icell_amr)
                     !phi_nbor(nbors_ind) = phi(icell_amr)
                  else
                     icell_nbor_amr = igrid_nbor_amr + (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind)  = phi(icell_nbor_amr)
                 end if
               end do
            end do
         end if

         grad_sqrd = ((phi_nbor(tnbors+1)-phi_nbor(1))*oneoverdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(11)-phi_nbor( 2)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(12)-phi_nbor( 3)-phi_nbor(13))*oneoverfourdx)**2
         nu(1) = nu_function(SQRT(grad_sqrd)/umond_a0)
         grad_sqrd = ((phi_nbor(4)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor(10)+phi_nbor( 5)-phi_nbor( 8)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor(16)+phi_nbor( 6)-phi_nbor( 9)-phi_nbor( 3))*oneoverfourdx)**2
         nu(4) = nu_function(SQRT(grad_sqrd)/umond_a0)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 8)-phi_nbor( 1)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor(2))*oneoverdx)**2 + &
                     ((phi_nbor(18)+phi_nbor( 6)-phi_nbor(14)-phi_nbor( 3))*oneoverfourdx)**2
         nu(2) = nu_function(SQRT(grad_sqrd)/umond_a0)
         grad_sqrd = ((phi_nbor(10)+phi_nbor( 4)-phi_nbor(11)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(17)-phi_nbor( 3)-phi_nbor(15))*oneoverfourdx)**2
         nu(5) = nu_function(SQRT(grad_sqrd)/umond_a0)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 9)-phi_nbor( 1)-phi_nbor(13))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(15)-phi_nbor( 2)-phi_nbor(14))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor( 3))*oneoverdx)**2
         nu(3) = nu_function(SQRT(grad_sqrd)/umond_a0)
         grad_sqrd = ((phi_nbor(16)+phi_nbor( 4)-phi_nbor(12)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor(17)+phi_nbor( 5)-phi_nbor(18)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)-phi_nbor(tnbors+1))*oneoverdx)**2
         nu(6) = nu_function(SQRT(grad_sqrd)/umond_a0)

         phi_vals = phi_nbor(1:twondim)
         nb_sum = nu(1)*phi_vals(1)+nu(2)*phi_vals(2)+nu(3)*phi_vals(3) + &
                  nu(4)*phi_vals(4)+nu(5)*phi_vals(5)+nu(6)*phi_vals(6)
         factor = nu(1)+nu(2)+nu(3)+nu(4)+nu(5)+nu(6)

         rho_pdm(icell_amr) = (nb_sum - factor*phi_nbor(tnbors+1))*oneoverfourpidx2

      end do
    end do

end subroutine make_phantom_halo

subroutine assign_newtonian
   use poisson_commons
   implicit none

   phi = newtonian_phi
   phi_old = newtonian_phi

end subroutine assign_newtonian

subroutine update_newtonian
   use poisson_commons
   implicit none

   newtonian_phi = phi

end subroutine update_newtonian

subroutine assign_mondian
   use poisson_commons
   implicit none

   phi = mondian_phi
   phi_old = mondian_phi

end subroutine assign_mondian

subroutine update_mondian
   use poisson_commons
   implicit none

   mondian_phi = phi

end subroutine update_mondian

! This nu function is not the same as the inverse mu function!
! This is the inverse mu function that appears in the standard formulation of QUMOND
! minus 1, so that the equation can be written as the addition of two density components
real(dp) function nu_function(y)
   implicit none

   real(dp), intent(in) :: y

   nu_function = -1.0d0 + (0.5d0 + 0.5d0*DSQRT(1.0d0 + 4.0d0/(y**mond_n)))**(1.0d0/mond_n)
   !nu_function = 0.5d0*DSQRT(1.0d0 + 4.0d0/y) - 0.5d0

end function nu_function

end module mond_commons
