module M_non_adiabatic_coupling

  use M_configuraciones
  use M_assemble_blocks  ! Using this to be able to use block and dblock structures

  ! Implicit none is a good idea, always
  implicit none
  type T_NAC_den
     type(T_assemble_neighbors), pointer :: tdenmat (:)
  end type T_NAC_den
  type T_NAC_vars
     integer ::  nac_inpfile = 123
     integer :: ntransitions
     real, allocatable :: c_old(:,:)
     real, allocatable :: suma(:,:)  !< dij.V analytical
     real, allocatable :: sumb(:,:)  !< dij.V numerical
     real, allocatable :: dnac(:,:)
     real, allocatable :: dnac_old(:,:)
     complex, allocatable :: c_na(:,:,:)
     type(T_forces), pointer :: dij(:,:,:)
     type(T_NAC_den), pointer :: band(:,:)
     type(T_vector), pointer :: ratom_old(:)
     type(T_vector), pointer :: vatom_old(:)

  end type T_NAC_vars

contains

  subroutine NAC_initialize(nac_vars, s)
    implicit none
    !Arguments
    type(T_NAC_vars), intent(inout) :: nac_vars
    type(T-structure), intent(in) :: s

    integer nac_inpfile

    integer istate
    integer stage
    integer iband
    integer :: ikpoint, iband
    integer :: nfermi
    integer numorb_max,norbitals

    complex a0,a1

    nac_inpfile = nac_vars%nac_inpfile
    open (unit = nac_inpfile, file = 'mdet.inp', status = 'old')
    ! Note : ntransitions is equal to nele in the old code
    read (nac_inpfile,*) nac_vars%ntransitions
    ntransitions = nac_vars%ntransitions
    stage = 1
    allocate(nac_vars%dij(ntransitions,ntransitions,s%natoms))
    ! Initialize NACs to zero
    do iatom = 1, s%natoms
       ! band-structure interactions
       nac_vars%dij(:,:,iatom)%kinetic = 0.0d0
       nac_vars%dij(:,:,iatom)%vna = 0.0d0
       nac_vars%dij(:,:,iatom)%vxc = 0.0d0
       nac_vars%dij(:,:,iatom)%vnl = 0.0d0
       nac_vars%dij(:,:,iatom)%ewald = 0.0d0
       
       ! corrections to the force
       nac_vars%dij(:,:,iatom)%usr = 0.0d0
       nac_vars%dij(:,:,iatom)%pulay = 0.0d0
       
       ! three-center interactions
       nac_vars%dij(:,:,iatom)%f3naa = 0.0d0
       nac_vars%dij(:,:,iatom)%f3nab = 0.0d0
       nac_vars%dij(:,:,iatom)%f3nac = 0.0d0
       nac_vars%dij(:,:,iatom)%f3xca = 0.0d0
       nac_vars%dij(:,:,iatom)%f3xcb = 0.0d0
       nac_vars%dij(:,:,iatom)%f3xcc = 0.0d0
       nac_vars%dij(:,:,iatom)%ftot  = 0.0d0
    end do
    allocate(nac_vars%band(ntransitions,ntransitions))
    allocate(nac_vars%c_na(ntransition,ntransitions,1))
! allocating map_ks, map_proj and iocc, these are in the nonadiabatic.f90 module in the old code, i have to create them here later
    allocate(nac_vars%map_ks(nac_vars%ntransitions))
    allocate(nac_vars%map_proj(nac_vars%ntransitions))
    allocate(nac_vars%iocc(nac_vars%ntransitions))
    ! Reading the transitions from mdet.inp file
    do istate = 1, nac_vars%ntransitions
       read (nac_inpfile,*) iband, nac_vars%iocc(istate)
       nac_vars%map_ks(istate) = iband
    end do
    close(nac_vars%nac_inpfile)
! allocate dnac
    allocate (nac_vars%dnac(nac_vars%ntransitions, nac_vars%ntransitions))
    allocate (nac_vars%dnac_old(nac_vars%ntransitions, nac_vars%ntransitions))
   
 ! Need allocation for imdet = 2, deal with it later
    allocate(nac_vars%c_na(nac_vars%ntransitions,nac_vars%ntransitions,nkpoints))
    allocate(nac_vars%ratom_old(natoms))
    allocate(nac_vars%vatom_old(natoms))
! These three need modification, and they might be allocated at some othe place in density, check them again
    allocate(s%kpoints(1)%eigen_old(s%norbitals))
    allocate(s%kpoints(1)%eigen(s%norbitals))
    allocate(s%kpoints(1)%c_Lowdin(s%norbitals,s%norbitals))
    allocate(nac_vars%eigen_1(nac_vars%ntransitions,nkpoints))
    allocate(nac_vars%eigen_0(nac_vars%ntransitions,nkpoints))
    ! MOved from M_density_MDET.f90 to here

    nfermi = int(s%ztot) / 2

   
    do iband = 1, nfermi
       nac_vars%foccupy_na (iband,ikpoint) = 1.0d0
       nac_vars%ioccupy_na (iband, ikpoint) = 2
       nac_vars%foccupy_na_TS (iband,ikpoint) = 1.0d0
       nac_vars%ioccupy_na_TS (iband, ikpoint) = 2
    end do
    if (int(qztot) .gt. 2*nfermi) then
       do ikpoint = 1, nkpoints
          nac_vars%ioccupy_na(nfermi+1,ikpoint) = 1
          nac_vars%foccupy_na (iband,ikpoint) = 0.5d0
          nac_vars%ioccupy_na_TS(nfermi+1,ikpoint) = 1
          nac_vars%foccupy_na_TS (iband,ikpoint) = 0.5d0
       end do
    end if

    call NAC_io(n, stage)

  end subroutine NAC_initialize

  subroutine NAC_finalize(n)
    implicit none
    type(T_NAC_vars), target :: n
    deallocate(nac_vars%dij)
  
    
  end subroutine NAC_finalize


  subroutine NAC_io(n, stage)

    implicit none
    type(T_NAC_vars), intent(inout) :: n
    integer :: stage

    if (stage == 1) then
       write(*,*) 'ilogfile', ilogfile

       write (ilogfile,*)
       write (ilogfile,*) ' Reading: mdet.inp '
       write (ilogfile,*) ' Number of transitions', n%ntransitions

    else if (stage == 'degenerate bands') then
      write(ilogfile,*)'TWO EIGENVALUES VERY CLOSE'
      write(ilogfile,*)'band', iband, eigen_k(map_ks(iband),ikpoint)
      write(ilogfile,*)'band', jband, eigen_k(map_ks(jband),ikpoint)
      write(ilogfile,*)'The nonadiabatic coupling is'
      write(ilogfile,*)'NOT CALCULATED'
    end if
    !        if (stage == 2) then
    !            write (ilogfile,*)
    !            write (ilogfile,*)'Call SCF_LOOP'
    !        end if
    !        if (stage == 3) then
    !            write (ilogfile,*)
    !            write (ilogfile,8)'Call getenery'





    !         write(*,*) 'The value of the elemnt (', i, i,') is:', this%m(i,i)


    ! Here you print on screen or ilogfile information about your MD

  end subroutine NAC_io


  subroutine NAC_fileio
    ! Here you should write info in those files  with massive data
    ! They should be written on HDF5, I will teach you this later
    ! Need to use the HDF5

  end subroutine NAC_fileio

  subroutine NAC_change_occupations(nac_vars)
    type(T_NAC_vars), intent(out) :: nac_vars

    integer ikpoint, iele, iband

    do ikpoint = 1, nac_vars%nkpoints
       do iele = 1, nac_vars%ntransitions
          iband = nac_vars%map_ks(iele)
          nac_vars%foccupy_na(iband,ikpoint) = nac_vars%iocc(iele)*0.5d0
          nac_vars%ioccupy_na(iband,ikpoint) = nac_vars%iocc(iele)
          nac_vars%foccupy_na_TS(iband,ikpoint) = nac_vars%iocc(iele)*0.5d0
          nac_vars%ioccupy_na_TS(iband,ikpoint) = nac_vars%iocc(iele)
       end do
    end do

  end subroutine NAC_change_occupations

  subroutine NAC_density_check(this, s)
    type(T_density_MDET), intent(out) :: this
    type(T_structure), intent(in) :: s
    real :: qcharge = 0.0d0
    integer :: ikpoint, iband

    do ikpoint = 1, this%nkpoints
       do iband = 1, this%norbitals
          if (ioccupy_na(iband,ikpoint) .ne. 0) then
             qcharge = qcharge + 2.0d0*foccupy_na(iband,ikpoint)*s%kpoints(ikpoint)%weight
          end if
       end do
    end do
    if (abs(qcharge - this%qztot) .gt. tol) then
       write (*,*) '          qcharge = ', qcharge
       write (*,*) '          qztot = ', this%qztot
       write (*,*) 'must stop in subroutine init_mdet 1'
       stop
    end if

  end subroutine NAC_density_check


  subroutine NAC_build_dij(nac_vars,s)
    type(T_NAC_vars), intent(inout) :: nac_vars
!    type(T_structure), intent(in) :: s

       implicit none
        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom                   !< counter over atoms/neighbors

! Allocate Arrays
! ===========================================================================
! Forces are stored in a Type with each piece, this makes acessing them and use
! pretty easy across the game.

! Procedure
! ===========================================================================

!!$
!!$! Format Statements
!!$! ===========================================================================
!!$! None
!!$
!!$! End Subroutine
!!$! ===========================================================================
!!$        return
!!$        end subroutine initialize_forces


! ===========================================================================
! build_forces
! ===========================================================================
! Subroutine Description
! ===========================================================================
!> This subroutine builds the total forces by adding contributions from the
!! kinetic, Vna, etc. and stores to a T_force variable called forces (:).
!
! ===========================================================================
! Code written by:
!> @author 
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        implicit none

        include '../include/constants.h'

! Argument Declaration and Description
! ===========================================================================
        type(T_structure), target :: s            !< the structure to be used

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, ineigh, matom, ialpha  !< counter over atoms/neighbors
        integer in1, in2, in3        !< species number
        integer jatom, num_neigh     !< counters over neighbors
        integer mbeta                !< the cell containing neighbor of iatom
        integer norb_mu, norb_nu     !< size of the (mu, nu) block for pair
        integer ix                   !< counter over dimensions
        integer imu, inu             !< counter over MEs
        integer mneigh
        integer ibeta, jbeta

        real sumT

        real, dimension (3) :: r1, r2, r3 !< position of atoms

        type(T_forces), pointer :: pdalpha
        type(T_forces), pointer :: pdi
        type(T_forces), pointer :: pdj

        ! band-structure interactions
        type(T_assemble_block), pointer :: pK_neighbors
        type(T_assemble_neighbors), pointer :: pkinetic
        type(T_assemble_block), pointer :: pvna_neighbors
        type(T_assemble_neighbors), pointer :: pvna
        type(T_assemble_block), pointer :: pvxc_neighbors
        type(T_assemble_neighbors), pointer :: pvxc
        type(T_assemble_block), pointer :: pSR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldsr
        type(T_assemble_block), pointer :: pLR_neighbors
        type(T_assemble_neighbors), pointer :: pewaldlr

        ! for overlap repulsive force
        type(T_assemble_block), pointer :: pCape_neighbors
        type(T_assemble_neighbors), pointer :: pcapemat
        type(T_assemble_block), pointer :: poverlap_neighbors
        type(T_assemble_neighbors), pointer :: poverlap

        ! Density matrix stuff
        type(T_assemble_neighbors), pointer :: pdenmat
        type(T_assemble_block), pointer :: pRho_neighbors
        type(T_assemble_block), pointer :: pRho_neighbors_matom

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
!       T W O - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************
! loop over atoms in central cell
        do iband, nac_vars%ntransitions
           do jband, nac_vars%ntransitions
              do iatom = 1, s%natoms
                 matom = s%neigh_self(iatom)
                 in1 = s%atom(iatom)%imass
                 norb_mu = species(in1)%norb_max
                 num_neigh = s%neighbors(iatom)%neighn
                 
                 ! cut some lengthy notation
                 pdi=>nac_vars%dij(iband,jband,iatom)
                 
                 ! density matrix
                 pdenmat=>nac_vars%band(iband,jband)%tdenmat(iatom)
                 pTRho_neighbors_matom=>pdenmat%neighbors(matom)
                 
                 ! interactions for each contribution
                 pkinetic=>s%kinetic(iatom)
                 pvna=>s%vna(iatom)
                 pvxc=>s%vxc(iatom)
                 pewaldsr=>s%ewaldsr(iatom)
                 pewaldlr=>s%ewaldlr(iatom)
                 
! allocate force terms and initialize to zero
                 allocate (pdi%vna_atom (3, num_neigh)); pdi%vna_atom = 0.0d0
                 !         allocate (pfi%vna_ontop (3, num_neigh)); pfi%vna_ontop = 0.0d0
                 allocate (pdi%vxc_off_site (3, num_neigh)); pdi%vxc_off_site = 0.0d0
                 allocate (pdi%vxc_on_site (3, num_neigh)); pdi%vxc_on_site = 0.0d0
                 allocate (pdi%ewaldsr (3, num_neigh)); pdi%ewaldsr = 0.0d0
                 allocate (pdi%ewaldlr (3, num_neigh)); pdi%ewaldlr = 0.0d0

! Now loop over all neighbors ineigh of iatom.
                 do ineigh = 1, num_neigh
                    mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                    jatom = s%neighbors(iatom)%neigh_j(ineigh)
                    in2 = s%atom(jatom)%imass
                    norb_nu = species(in2)%norb_max

                    ! cut some lengthy notation
                    pdj=>nac_vars%dij(iband,jband,jatom)

                    ! density matrix - neighbors
                    pTRho_neighbors=>pdenmat%neighbors(ineigh)
                    
                    ! interactions - neighbors
                    pK_neighbors=>pkinetic%neighbors(ineigh)
                    pvna_neighbors=>pvna%neighbors(ineigh)
                    
! KINETIC FORCES (TWO-CENTER)
! ***************************************************************************
! The derivatives are tpx and, where p means derivative and x means crytal
! coordinates. The derivative is a vector in crystal
! coordinates and is stored in pK_neighbors%Dblock. The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4).
                    do ix = 1, 3
                       sumT = 0.0d0
                       do inu = 1, norb_nu
                          do imu = 1, norb_mu
                             sumT = sumT                                                &
                                   pTRho_neighbors%block(imu,inu)*pK_neighbors%Dblock(ix,imu,inu)
                          end do
                       end do
                       
! Now add sum to appropriate force term. see notes "the total band structure
! The (-1.d0) makes it "force-like".
                       ! direct term
                       pdi%kinetic(ix) = pdi%kinetic(ix) + (-1.0d0)*sumT
                       ! cross term
                       pdj%kinetic(ix) = pdj%kinetic(ix) - (-1.0d0)*sumT
                    end do ! do ix


! ASSEMBLE HARTREE (TWO-CENTER) FORCES - ATOM CASE
! ***************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
!
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.

! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Now loop over all neighbors ineigh of iatom.
! Notice the explicit negative sign, this makes it force like.
                    do inu = 1, norb_mu
                       do imu = 1, norb_mu
                          pdi%vna_atom(:,ineigh) = pdi%vna_atom(:,ineigh)             &
                               &           - pTRho_neighbors_matom%block(imu,inu)*pvna_neighbors%Dblock(:,imu,inu)
                       end do
                    end do
                 end do ! end loop over neighbors


! ASSEMBLE HARTREE (TWO-CENTER) FORCES - ONTOP CASE
! ***************************************************************************
! Now loop over all neighbors ineigh of iatom.
                 do ineigh = 1, num_neigh
                    mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                    jatom = s%neighbors(iatom)%neigh_j(ineigh)
                    in2 = s%atom(jatom)%imass
                    norb_nu = species(in2)%norb_max
                    
                    ! cut some lengthy notation
                    pdj=>s%nac_vars%dij(iband,jband,jatom)

                    ! density matrix - neighbors
                    pTRho_neighbors=>ptdenmat%neighbors(ineigh)

                    ! cut some lengthy notation for vna
                    pvna_neighbors=>pvna%neighbors(ineigh)

! If r1 .ne. r2, then this is a case where the potential is located at one of
! the sites of a wavefunction (ontop case).
                    if (iatom .eq. jatom .and. mbeta .eq. 0) then
! Do nothing here - special case. Interaction already calculated in atm case.

                    else

! Notice the explicit negative sign, this makes it force like.
                       do inu = 1, norb_nu
                          do imu = 1, norb_mu
                             pdi%vna_ontop(:,ineigh) = pdi%vna_ontop(:,ineigh)          &
                                  &            - pTRho_neighbors%block(imu,inu)*pvna_neighbors%Dblocko(:,imu,inu)
                          end do
                       end do
                    end if
                 end do ! end loop over neighbors


! ASSEMBLE EXCHANGE-CORRELATION (TWO-CENTER) FORCE - ON-SITE CASE
! ***************************************************************************
! The vxc two-center forces are: vxc_on_site and vxc_off_site.

! First we calculate the on-site force contributions.
! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vxc_on_site term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.

! Now loop over all neighbors ineigh of iatom.
                 pvxc_neighbors=>pvxc%neighbors(matom)
                 do ineigh = 1, num_neigh            
                    do inu = 1, norb_mu
                       do imu = 1, norb_mu
                          pdi%vxc_on_site(:,ineigh) = pdi%vxc_on_site(:,ineigh)        &     
                               &           - pTRho_neighbors_matom%block(imu,inu)*pvxc_neighbors%Dblocko(:,imu,inu)
                       end do
                    end do
                 end do ! end loop over neighbors


! ASSEMBLE EXCHANGE-CORRELATION (TWO-CENTER) FORCE - OFF-SITE CASE
! ***************************************************************************
! Next, we calculate the off site force interaction terms.
! If r1 = r2, then this is a case of the self-interaction term or the
! one center term which has no force.
!
! Form the "force-like" derivative of the atom terms for vxc,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Now loop over all neighbors ineigh of iatom.
                 do ineigh = 1, num_neigh
                    mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                    jatom = s%neighbors(iatom)%neigh_j(ineigh)
                    in2 = s%atom(jatom)%imass
                    norb_nu = species(in2)%norb_max
                    
                    ! cut some lengthy notation
                    pdj=>nac_vars%dij(iband,jband,jatom)

                    ! density matrix - neighbors
                    pTRho_neighbors=>ptdenmat%neighbors(ineigh)

                    ! vxc interactions - neighbors
                    pvxc_neighbors=>pvxc%neighbors(ineigh)

! Notice the explicit negative sign, this makes it force like.
                    if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in atm case.

                    else
                       do inu = 1, norb_nu
                          do imu = 1, norb_mu
                             pdi%vxc_off_site(:,ineigh) = pdi%vxc_off_site(:,ineigh)    &
                                  &             - pTRho_neighbors%block(imu,inu)*pvxc_neighbors%Dblock(:,imu,inu)
                          end do
                       end do
                    end if
                 end do ! end loop over neighbors
                 

!  ASSEMBLE EWALD (TWO-CENTER) FORCES
! ***************************************************************************
! The Ewald two-center forces are: ewaldsr and ewaldlr.
!
! If r1 = r2, then this is a case of the self-interaction term or the
! one center term which has no force.
!
! Form the "force-like" derivative of the atom terms for vxc,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.

! Now loop over all neighbors ineigh of iatom.
                 do ineigh = 1, num_neigh
                    mbeta = s%neighbors(iatom)%neigh_b(ineigh)
                    jatom = s%neighbors(iatom)%neigh_j(ineigh)
                    in2 = s%atom(jatom)%imass
                    norb_nu = species(in2)%norb_max

                    ! cut some lengthy notation
                    pdj=>nac_vars%dij(iband,jband,jatom)

                    ! density matrix - neighbors
                    pTRho_neighbors=>ptdenmat%neighbors(ineigh)
                    
                    ! interactions - neighbors
                    pSR_neighbors=>pewaldsr%neighbors(ineigh)
                    pLR_neighbors=>pewaldlr%neighbors(ineigh)

! Notice the explicit negative sign, this makes it force like.
                    if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special self-interaction case.

                    else

! short-range part ewaldsr
                       do inu = 1, norb_nu
                          do imu = 1, norb_mu
                             pdi%ewaldsr(:,ineigh) = pdi%ewaldsr(:,ineigh)              &
                                  &             - 0.5d0*pTRho_neighbors%block(imu,inu)*pSR_neighbors%Dblock(:,imu,inu)
! Note - remove the 0.5d0 and make sure it gets into the Dassembler - I add it here
! because the 0.5d0 was here in the original assemble_F.f90 routine.
                          end do
                       end do

! long-range part ewaldsr
                       do inu = 1, norb_nu
                          do imu = 1, norb_mu
                             pdi%ewaldlr(:,ineigh) = pdi%ewaldlr(:,ineigh)              &
                                  &             - pTRho_neighbors%block(imu,inu)*pLR_neighbors%Dblock(:,imu,inu)
                          end do
                       end do
                    end if
                 end do ! end loop over neighbors
              end do ! end loop over atoms


! ADD CONTRIBUTIONS TO GET TOTAL BAND-STRUCTURE FORCE (TWO-CENTER)
! ***************************************************************************
! loop over atoms in central cell
              do iatom = 1, s%natoms
                 ! cut some lengthy notation
                 pdi=>nac_vars%dij(iband,jband,iatom)

! kinetic contribution to total force
! ****************************************************************************
                 pdi%ftot = pdi%ftot + pdi%kinetic

! Loop over all neighbors of iatom and add in the neighbor-contributed forces
! ****************************************************************************
                 num_neigh = s%neighbors(iatom)%neighn
                 do ineigh = 1, num_neigh

! cut some lengthy notation
                    jatom = s%neighbors(iatom)%neigh_j(ineigh)
                    pdj=>nac_vars%dij(iband,jband,jatom)

! vna contribution to total force
! ****************************************************************************
! Hartree forces - atom case
                    pdi%vna = pdi%vna + pdi%vna_atom(:,ineigh)
                    pdj%vna = pdj%vna - pdi%vna_atom(:,ineigh)

                    pdi%ftot = pdi%ftot + pdi%vna_atom(:,ineigh)
                    pdj%ftot = pdj%ftot - pdi%vna_atom(:,ineigh)

! Hartree forces - ontop terms
                    pdi%vna = pdi%vna + pdi%vna_ontop(:,ineigh)
                    pdj%vna = pdj%vna - pdi%vna_ontop(:,ineigh)

                    pdi%ftot = pdi%ftot + pdi%vna_ontop(:,ineigh)
                    pdj%ftot = pdj%ftot - pdi%vna_ontop(:,ineigh)
!
! vxc contribution to total force
! ****************************************************************************
! off site interactions
                    pdi%vxc = pdi%vxc + pdi%vxc_off_site(:,ineigh)
                    pdj%vxc = pdj%vxc - pdi%vxc_off_site(:,ineigh)

                    pdi%ftot = pdi%ftot + pdi%vxc_off_site(:,ineigh)
                    pdj%ftot = pdj%ftot - pdi%vxc_off_site(:,ineigh)

! on site interactions
                    pdi%vxc = pdi%vxc + pdi%vxc_on_site(:,ineigh)
                    pdj%vxc = pdj%vxc - pdi%vxc_on_site(:,ineigh)

                    pdi%ftot = pdi%ftot + pdi%vxc_on_site(:,ineigh)
                    pdj%ftot = pdj%ftot - pdi%vxc_on_site(:,ineigh)

! Ewald contribution to total force
! ****************************************************************************
! ewaldsr interactions
                    pdi%ewald = pdi%ewald - pdi%ewaldsr(:,ineigh)
                    pdj%ewald = pdj%ewald + pdi%ewaldsr(:,ineigh)

                    pdi%ftot = pdi%ftot - pdi%ewaldsr(:,ineigh)
                    pdj%ftot = pdj%ftot + pdi%ewaldsr(:,ineigh)

! ewaldlr interactions
                    pdi%ewald = pdi%ewald - pdi%ewaldlr(:,ineigh)
                    pdj%ewald = pdj%ewald + pdi%ewaldlr(:,ineigh)

                    pdi%ftot = pdi%ftot - pdi%ewaldlr(:,ineigh)
                    pdj%ftot = pdj%ftot + pdi%ewaldlr(:,ineigh)
                 end do ! end loop over neighbors
              end do ! end loop over atoms

! Vnl contribution to total force
! ****************************************************************************
! Loop over all atoms iatom in the central cell.
              do iatom = 1, s%natoms

! cut some lengthy notation
                 pdi=>nac_vars%dij(iband,jband,iatom)

! Loop over all neighbors of iatom and add in the neighbor-contributed forces
! Note - the neighbor mapping for vnl is different than neighbor mapping for
! other terms, so, we need to add in the contributions correctly.
                 do ineigh = 1, s%neighbors_PP(iatom)%neighn
                    jatom = s%neighbors_PP(iatom)%neigh_j(ineigh)

! cut some lengthy notation
                    pdj => nac_vars%dij(iband,jband,jatom)

! atom contribution
                    pdi%vnl = pdi%vnl + pdi%vnl_atom(:,ineigh)
                    pdj%vnl = pdj%vnl - pdi%vnl_atom(:,ineigh)

                    pdi%ftot = pdi%ftot + pdi%vnl_atom(:,ineigh)
                    pdj%ftot = pdj%ftot - pdi%vnl_atom(:,ineigh)
                 end do

! ontop left contribution
                 do ineigh = 1, s%neighbors_PPx(iatom)%neighn
                    jatom = s%neighbors_PPx(iatom)%neigh_j(ineigh)

! cut some lengthy notation
                    pdj=>nac_vars%dij(iband,jband,jatom)

                    pdi%vnl = pdi%vnl + 2.0d0*pdi%vnl_ontop(:,ineigh)
                    pdj%vnl = pdj%vnl - 2.0d0*pdi%vnl_ontop(:,ineigh)

                    pdi%ftot = pdi%ftot + 2.0d0*pdi%vnl_ontop(:,ineigh)
                    pdj%ftot = pdj%ftot - 2.0d0*pdi%vnl_ontop(:,ineigh)
                 end do ! end loop over neighbors
              end do  ! end loop over atoms
! ***************************************************************************
!                                   E N D
!       T W O - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************

! ADD CONTRIBUTIONS TO GET TOTAL BAND-STRUCTURE FORCE (THREE-CENTER)
! ***************************************************************************
! Loop over all atoms iatom in the central cell.
! Single-source loops (not dependent on neighbours)
              do iatom = 1, s%natoms
                 ! cut some lengthy notation
                 pdi=>nac_vars%dij(iband,jband,iatom)

! vna three-center contribution to the total force
! ****************************************************************************
                 pdi%ftot = pdi%ftot - pdi%f3naa - pdi%f3nab - pdi%f3nac

! vxc three-center contribution to the total force
! ****************************************************************************
                 pdi%ftot = pdi%ftot - pdi%f3xca - pdi%f3xcb - pdi%f3xcc
              end do ! end loop over atoms
! ***************************************************************************
!                                   E N D
!     T H R E E - C E N T E R   B A N D - S T R U C T U R E   F O R C E S
! ***************************************************************************

! ***************************************************************************
!
! Need to add the Ei<d/dr (phi)_mu|(phi)_nu>+Ej< (phi)_mu|d/dr(phi)_nu>
!            P U L A Y   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          in1 = s%atom(iatom)%imass
          norb_mu = species(in1)%norb_max
          num_neigh = s%neighbors(iatom)%neighn

          ! cut some lengthy notation
          pfi=>s%forces(iatom)

          ! density matrix with eigenvalues
          pcapemat=>s%capemat(iatom)

          ! interactions for each contribution
          poverlap=>s%overlap(iatom)

! Now loop over all neighbors ineigh of iatom.
          do ineigh = 1, num_neigh
            mbeta = s%neighbors(iatom)%neigh_b(ineigh)
            jatom = s%neighbors(iatom)%neigh_j(ineigh)
            in2 = s%atom(jatom)%imass
            norb_nu = species(in2)%norb_max

            ! cut some lengthy notation
            pfj=>s%forces(jatom)

            ! density matrix - neighbors
            pCape_neighbors=>pcapemat%neighbors(ineigh)

            ! interactions - neighbors
            poverlap_neighbors=>poverlap%neighbors(ineigh)

! The derivatives are tpx and, where p means derivative and x means crytal
! coordinates. The derivative is a vector in crystal
! coordinates and is stored in pK_neighbors%Dblock. The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4).
            do ix = 1, 3
              sumT = 0.0d0
              do inu = 1, norb_nu
                do imu = 1, norb_mu
                  sumT = sumT                                                &
                   + pCape_neighbors%block(imu,inu)*poverlap_neighbors%Dblock(ix,imu,inu)
                end do
              end do

! Now add sum to appropriate force term. see notes "the total band structure
! The (-1.d0) makes it "force-like".
              ! direct term
              pfi%pulay(ix) = pfi%pulay(ix) + (-1.0d0)*sumT
              ! cross term
              pfj%pulay(ix) = pfj%pulay(ix) - (-1.0d0)*sumT
            end do ! do ix
          end do ! end loop over neighbors
        end do ! end loop over atoms

! ADD CONTRIBUTIONS TO GET TOTAL FORCE AFTER PULAY CORRECTION
! ***************************************************************************
! loop over atoms in central cell
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)

! overlap repulsive contribution to total force
! ****************************************************************************
          pfi%ftot = pfi%ftot - pfi%pulay
        end do
! ***************************************************************************
!                                  E N D
!            P U L A Y   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************

! ***************************************************************************
!
!            U S R   C O R R E C T I O N S   (T W O - C E N T E R)
! ***************************************************************************
! Loop over all atoms iatom in the central cell.
        do iatom = 1, s%natoms
          ! cut some lengthy notation
          pfi=>s%forces(iatom)
          pfi%ftot = pfi%ftot + pfi%usr
        end do ! end loop over atoms

! Format Statements
! ===========================================================================
! None

! End Subroutine
! ===========================================================================
   return
 end subroutine NAC_build_dij
 
 subroutine overlap_sign (nac_vars, s, itime_step,species)
   integer,intent(in) :: itime_step 
   integer iband, iorbital 
   type(T_NAC_vars), intent(inout) :: nac_vars
   type(T_structure), intent(in) :: s
   type(T_species), intent(in) :: species(:)
   
   call ovelap_numeric(nac_vars,s,itime_step,species) ! I think this should not be in this subroutine it should be called from the main loop
   do iband = 1, nac_vars%ntransitions
      if (nac_vars(iband,iband) .lt. -0.1) then
         ! I have commented the write part this should be moved to the io subroutine
         !            write (*,*) 'The arbitrary sign of wavefunctions of this and previous time step are different'
         !            write (*,*) 'We will change the sign in order to have the same sign in all time steps!!'
         do iorbital = 1 , s%norbitals
            s%kpoints(1)%c(iorbital,nac_vars%map_ks(iband)) = - s%kpoints(1)%c(iorbital,nac_vars%map_ks(iband))
         end do
      end if
   end do
 end subroutine overlap_sign

 subroutine overlap_numeric(nac_vars, s , itime_step, species)
    implicit none 
    include '../include/interactions_2c.h'
    integer it
    integer imu, inu, jmu,jnu, iorbital !< counter over orbitals
    integer in1, in2, in3               !< species numbers
    integer iatom, jatom, ineigh        !< counter over atoms and neighbors
    integer isorp, interaction          !< which interaction and subtype
    integer ix
    integer ia
    integer ik,ij
    integer num_neigh                   !< number of neighbors
    integer mbeta                       !< the cell containing iatom's neighbor
    
    integer norb_mu, norb_nu         !< size of the block for the pair

    real z                           !< distance between r1 and r2
    real diff
    real delta

    real, dimension (3) :: eta        !< vector part of epsilon eps(:,3)
    real, dimension (3, 3) :: eps     !< the epsilon matrix
    real, dimension (3, 3, 3) :: deps !< derivative of epsilon matrix
    real, dimension (3) :: r1, r2, r3, v     !< positions of iatom and jatom
    real, dimension (3) :: sighat     !< unit vector along r2 - r1
    real, dimension (nac_vars%ntransitions, nac_vars%ntransitions) :: suma, sumb
    real y, rcutoff_i, rcutoff_j, range
    real, dimension(s%norbitals,s%norbitals) :: sover
    real, dimension(numorb_max,numorb_max) :: sm ! need to find the numorb_max
    real, dimension(3,numorb_max,numorb_max) :: sx ! need to find the numorb_max in the new fireball or change the whole method

    interface
       function distance (a, b)
         real distance
         real, intent(in), dimension (3) :: a, b
       end function distance
    end interface
    type(T_NAC_vars), intent(inout) :: nac_vars
    type(T_structure), intent(in) :: s
    type(T_species), intent(in) :: species(:)
    

    type(T_assemble_block), pointer :: pS_neighbors
    type(T_assemble_neighbors), pointer :: poverlap


    do iatom = 1, s%natoms
       r1 = nac_vars%ratom_old(iatom)%a
       rcutoff_i = 0.0d0
       in1 = s%atom(iatom)%imass
       norb_mu = species(in1)%norb_max
       do imu =1 , species(in1)%nssh
          if (species(in1)%shell(imu)%rcutoff .gt. rcutoff_i) rcutoff_i = species(in1)%shell(imu)%rcutoff
       end do! imu
       do jatom = 1, s%natoms
          r2 = s%atom(jatom)%ratom
          in2 = s%atom(iatom)%imass
          norb_nu = species(in2)%norb_max
          z = distance(r1,r2)
          rcutoff_j = 0.0d0
          do imu = 1 ,species(in2)%nssh  
             if (species(in2)%shell(imu)%rcutoff .gt. rcutoff_i) rcutoff_j = species(in2)%shell(imu)%rcutoff
          end do
          range = abs(rcutoff_i + rcutoff_j - 0.01d0)
          if (z .gt. range) then
             sover = 0
          else
             if (z .lt. 1.0d-05) then
                sighat(1) = 0.0d0
                sighat(2) = 0.0d0
                sighat(3) = 1.0d0
             else 
                sighat = (r2 - r1)/z
             end if
             call epsilon_function (r2, sighat, eps)
             isorp = 0
             interaction = P_overlap
             in3 = in2
             call getMEs_Fdata_2c (in1, in2, interaction, isorp, z,        &
                  &                            norb_mu, norb_nu, sm)
             call rotate (in1, in3, eps, norb_mu, norb_nu, sm, sx)
             do imu = 1, norb_mu
                jmu = imu +  s%iblock_slot(iatom)
                do inu = 1, norb_nu 
                   jnu = inu +  s%iblock_slot(jatom)
                   sover(jmu, jnu) = sx(imu,inu)
                end do ! end imu
             end do ! end inu
          end if
       end do ! jatom
    end do !iatom 
       
! ===========================================================================
! Calculate overlap between Kohn-Sham states at different
! time steps
! Non-adiabatic term: dot pruduct sum
    nac_vars%sumb = 0.0d0
    do ij = 1, nac_vars%ntransitions
       do ik = 1, nac_vars%ntransitions
          do imu = 1, s%norbitals
             do inu = 1, s%norbitals
                nac_vars%sumb(ik,ij) = nac_vars%sumb(ik,ij) +                                 & 
                     &      nac_vars%c_old(mmu,nac_vars%map_ks(ij))* s%kpoints(1)%c(inu,nac_vars%map_ks(ij))* sover(imu,inu)
             end do
          end do
       end do
    end do
    
    do ij = 1 , nac_vars%ntransitions
       do ik = 1, nac_vars%ntransitions
          nac_vars%dnac(ik,ij) = ((nac_vars%sumb(ik,ij)) - nac_vars%sumb(ij,ik))/(2.0d0*dt)
       end do
    end do


! analytical method, this can be added as an option, suma(ik,ij)= dij.V , 
    delta = 0.5d0
    nac_vars%suma = 0.0d0
    do ik = 1 , nac_vars%ntransitions
       do ij = 1 , nac_vars%ntransitions
          do iatom = 1 , s%natoms
             r1 = nac_vars%ratom_old(iatom)%a
             r2 = s%atom(iatom)%ratom
             v = (r2 - r1)/dt 
             do ix = 1, 3
                nac_vars%suma(ik,ij) = nac_vars%suma(ik,ij) + v(ix)*nac_vars%dij(ik,ij,iatom)%ftot(ix)
             end do
          end do
       end do
    end do    
 
  end subroutine overlap_numeric
  
  subroutine evolve_ks_state(nac_vars,den_vars,s)

  end subroutine evolve_ks_state  
  


  subroutine FSSH(nac_vars,s, itime_step)
    integer ij, ik, iswitch
    real xrand, aux, ajj, bkj
    complex akj
    real, dimansion(nac%vars) :: prob
    type(T_NAC_vars), intent(inout) :: nac_vars
    type(T_structure), intent(in) :: s
    
! ===========================================================================
! Calculate hopping probabilities for fewest switches
! map_ks(iband) gives back the corresponding adiabatic KS state
! we follow the possible transitions associated with states
! the reason using ij and ik to match the formula from the paper
    do ij = 1, nac_vars%ntransitions
!----------------------------------------------------------
! Random numbers for Monte-Carlo
       call random_number(xrand)
       ajj = real(conjg(nac_vars%c_na(ij,ij,1))*nac_vars%c_na(ij,ij,1))
       do ik = 1 , nac_vars%ntransitions
          akj = nac_vars%c_na(ij,ik,1)*conjg(nac_vars%c_na(ij,ik,1))
          bkj = -2.0d0*real(conjg(akj)*nac_vars%dnac(ik,ij))
! JOM-warning: may be later we can "imporve" this by using eq(29) in JCP
! 101 4657 (1994)
!----------------------------------------------------------
! probability of the jband ---> kband transition
          prob(ik) = bkj*dt/ajj
          write(*,*)'prob',ij,ik,prob(ik) ! move this to io subroutine
          if (prob(ik) .lt. 0.0d0) then
             prob(ik) = 0.0d0
          end if
       end do ! ik
!----------------------------------------------------------
! Monte-Carlo calculation for possible transitions
! JOM-warning : we should also allow transitions to states that are not
! fully occupied (ioccupy_na = 0, 1) [from states that are occupied
! ioccupy_na = 1, 2 ]. Use iocc for this (fix later)
!         iocc (ik) = ioccupy_na (ik, ikpoint)
!----------------------------------------------------------
       call mc_switch (xrand, nac_vars%ntransitions, prob, ij, iswitch)
       
       if (iswitch .ne. 0) then
          write(*,*)'SWITCH!!',ij, '--->',iswitch ! move it to NACio
!----------------------------------------------------------
! perform transition ij ---> iswitch
         call transition (itime_step, ij, iswitch, ikpoint)   
!----------------------------------------------------------
         return  ! we can only have one switch
      end if
   end do ! ij
 
 end subroutine FSSH

 subroutine mc_switch(xr, ntransitions, prob, ij, is) 
   integer ij,is,ntransitions
   integer ij, ik 
   integer jband, kband
   real xr
   real aux
   real, dimenstions(ntransitions), intent(in) :: prob
! Procedure
! ===========================================================================
! Check that sum of probabilities is smaller than 1
   aux = 0.0d0
   do ik = 1, ntransitions
! Consider only allowed transitions
      jband = nac_vars%map_ks(ij)
      kband = nac_vars%map_ks(ik)
      if (nac_vars%ioccupy_na(jband,1) .gt. 0) then
         if (nac_vars%ioccupy_na(kband,1) .lt. 2) then
            aux = aux + prob(ik)
         end if
      end if
   end do
   if (aux .gt. 1.0d0) then ! this might never happen and move the write to NACio
      write(*,*)'sum of probabilities greater than 1'
      write(*,*)'in mc_switches.f90'
      write(*,*)'total probabilty', aux
      write(*,*)'for state', ij
      do ik = 1, n
         write(*,*)'prob',ik,prob(k)
      end do
              write(*,*)'must stop' ! this was commented in the old code!! don't know why.
              stop
   end if
   is = 0 
   aux = 0.0d0
! Consider only allowed transitions
   do ik = 1, ntransition
      jband = map_ks(ij)
      kband = map_ks(ik)
      if (nac_vars%ioccupy_na(jband,1) .gt. 0) then     
          if (nac_vars%ioccupy_na(kband,1) .lt. 2) then    
             aux = aux + pr(ik)
             if (aux .gt. xr) then
                is = ik
                do ij = 1, ntransition
                   write(*,*)'prob',ij,prob(ij)
                end do
                exit
             end if
          end if
       end if
    end do
!----------------------------------------------------------
! If is = 0, no switch (hopping) between states
! Otherwise, switch from current state ( "ij" = map_ks(iele) in
! fewest_switches subroutine) to state "is"
!----------------------------------------------------------
    return 
  end subroutine mc_switch

subroutine transition(itime_step, ij,is,s,nac_vars) ! there are 3 subroutine in old code by this name I have used transition.f90 according to make file
  integer, intent(in) :: itime_step, ij, is
  
  integer jband, kband
  integer iatom,ix
  real ejump, ener
  real aa,bb,cc
  real alpha
  real xm
  real tkinetic
  real, parameter :: tolaa = 1.0d-8
  ! introduce the pointers


  jband = nac_vars%map_ks(ij)
  kband = nac_vars%map_ks(is)
  s%kpoints(1)%foccupy(jband) = s%kpoints(1)%ioccupy(jband) - 1
  s%kpoints(1)%foccupy(jband) = s%kpoints(1)%foccupy(jband) - 0.50d0
  s%kpoints(1)%ioccupy(kband) = s%kpoints(1)%ioccupy(kband) + 1
  s%kpoints(1)%foccupy(kband) = s%kpoints(1)%foccupy(kband) + 0.50d0

  ejump = s%kpoints(1)%eigen(kband) - s%kpoints(1)%eigen(jband)
  write(*,*)'ejump (eV) =', ejump  ! move this to NACio
  ener = ejump*P_fovermp
  write(*,*)'ener (dynamical units)=', ener ! move this to NACio
! ===========================================================================
! Find out if transition ij --> is is accesible (i.e. if there is enough
! kinetic energy 
!----------------------------------------------------------
 ! I skiped a step which prints the dij again
  aa = 0.0d0
  bb = 0.0d0
  do iatom = 1, s%natoms
     ! cut some lengthy notation
     pdij => nac_vars%dij(ij,is,iatom)%ftot
     pvatom => s%atom(iatom)%vatom
     in1 = s%atom(iatom)%imass
     xm = species(in1)%xmass
     do ix = 1,3
        aa = aa + 0.50d0*( pdij**2 )/xm
        bb = bb + pvatom(ix)*pdij(ix) 
     end do ! ix
  end do ! iatom
!----------------------------------------------------------
  write(*,*)'aa, 4*ener*aa =', aa, 4.0d0*aa*ener ! move to NACio
  write(*,*)'bb, bb**2=', bb, bb**2 ! move to NACio
  cc = bb**2 - 4.0d0*aa*ener
  write(*,*)'cc=', cc !Move to NACio
! ===========================================================================
  if (aa .gt. tolaa) then
!----------------------------------------------------------
! the transition is accepted
     if (cc .ge. 0.0d0) then
        write(*,*)'transition accepted'
!----------------------------------------------------------
        if (bb .ge. 0.0d0) then
           alfa = (bb - sqrt(cc))/(2.0d0*aa)
        else
           alfa = (bb + sqrt(cc))/(2.0d0*aa)
        end if
!----------------------------------------------------------
     else !(cc .ge. 0.0d0)
! the transition is NOT accepted
        write(*,*)'transition NOT accepted' ! Move to NACio
!----------------------------------------------------------
        s%kpoints(1)%foccupy(jband) = s%kpoints(1)%ioccupy(jband) + 1
        s%kpoints(1)%foccupy(jband) = s%kpoints(1)%foccupy(jband) + 0.50d0
        s%kpoints(1)%ioccupy(kband) = s%kpoints(1)%ioccupy(kband) - 1
        s%kpoints(1)%foccupy(kband) = s%kpoints(1)%foccupy(kband) - 0.50d0
! Define alfa for re-scaling velocities
! JOM-info : two possibilities: a) Do nothing; b) Reflection
! a) Do nothing
! b) Reflection (see JCP 101, 4657 (1994), pag. 4664)
        alfa = bb/aa
     end if  !(cc .ge. 0.0d0)
  else  !(aa .gt. tolaa)
     alfa = 0.0d0
  end if  !(aa .gt. tolaa)
  write(*,*)'alfa=', alfa ! move to NACio
! ===========================================================================
  do iatom = 1 , s%natoms
     write(*,*)'vatom-B',  (s%atom(iatom)%vatom(ix) , ix = 1,3)
  end do
!----------------------------------------------------------
  tkinetic = 0.0d0
  do iatom = 1, s%natoms
     in1 = s%atom(iatom)%imass
     xm = species(in1)%xmass
     ! cut some lengthy notation
     pvatom => s%atom(iatom)%vatom
     tkinetic = tkinetic                                            &
          &       + (0.5d0/P_fovermp)*xm                             &
          &      *(pvatom(1)**2 + pvatom(2)**2 + pvatom(3)**2)
  end do
  write(*,*)'KINETIC=',tkinetic ! move to NACio
!----------------------------------------------------------

!----------------------------------------------------------
! RESCALING VELOCITIES
  do iatom = 1, s%natoms
     do ix = 1, 3
        ! cut some lengthy notation
        pdij => nac_vars%dij(ij,is,iatom)%ftot
        pvatom => s%atom(iatom)%vatom
        in1 = s%atom(iatom)%imass
        xm = species(in1)%xmass
        pvatom (ix) = pvatom (ix) - alfa*pdij(ix)/xm
     end do ! ix
  end do ! iatom
!----------------------------------------------------------
!----------------------------------------------------------
  do iatom = 1, s%natoms
     write(*,*)'vatom-A',  (s%atom(iatom)%vatom (ix), ix = 1,3)
  end do
!----------------------------------------------------------
!----------------------------------------------------------
  tkinetic = 0.0d0
  do iatom = 1, s%natoms
     in1 = s%atom(iatom)%imass
     xm = species(in1)%xmass
     ! cut some lengthy notation
     pvatom => s%atom(iatom)%vatom
     tkinetic = tkinetic                                            &
          &       + (0.5d0/P_fovermp)*xm                             &
          &      *(pvatom(1)**2 + pvatom(2)**2 + pvatom(3)**2)
  end do
  write(*,*)'KINETIC=',tkinetic ! move to NACio




end subroutine transition



!          subroutine NAC_do_something1(this)

  ! Working routines, they act on the variables defined NAC_vars and
  ! operate changes
  !            type(NAC_vars), intent(inout) :: this
  !            integer :: i

  !            do i=1, this%size

  !               this%m(i,i) = i**2

  !            end do


  !         end subroutine NAC_do_something1


  !        subroutine NAC_do_something2

  !        end subroutine NAC_do_something2

  subroutine find_neigh_max_NAC(s,neigh_max, neighPP_max) !< this should not be here but we need it for allocating gh arrays, we can change it later

  type(T_structure), intent(in) :: s
  integer neigh_max,neighPP_max ,iatom

  neigh_max = -99
  do iatom = 1, s%natoms
    neigh_max = max(neigh_max,size(s%neighbors(iatom)%neigh_j))
  end do
  do iatom = 1, s%natoms
    neighPP_max = max(neigh_max, size(s%neighbors_PP(iatom)%neigh_j))
  end do

  end subroutine find_neigh_max_NAC


  subroutine NAC_normalization(this, nkpoints)
    type(T_NAC_vars), intent(inout) :: this

    complex :: a0 = cmplx(0.0d0,0.0d0)
    complex :: a1 = cmplx(1.0d0,0.0d0)
    complex :: cnorm, caux
    integer :: ikpoint, iele, jele, nkpoints
    real :: norm

    this%c_na = a0
    do ikpoint = 1, nkpoints
       do iele = 1, this%ntransitions
          this%c_na (iele, iele, ikpoint) = a1
          !       write(*,*)'c_na',iband,c_na (iband, iband, ikpoint)
       end do
    end do

    ! check Normalization !
    do ikpoint = 1, nkpoints
       do iele = 1, this%ntransitions
          cnorm = a0
          do jele = 1, this%ntransitions
             caux = this%c_na(iele,jele,ikpoint)
             cnorm = cnorm + caux*conjg(caux)
          end do
          norm = cabs (cnorm)
          write(*,*)'Norm of initial states',iele, norm
       end do
    end do

  end subroutine NAC_normalization


end module M_non_adiabatic_coupling
