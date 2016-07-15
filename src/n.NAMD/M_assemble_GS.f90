subroutine assemble_G_S(nac_vars,s)
  use M_configuraciones
  type(T_structure), intent(in) :: s
  type(T_NAC_vars), intent(inout) :: nac_vars
  integer iatom
  integer matom
  integer num_neigh
  integer in1
  
  real r1
  
  type(T_assemble_neighbors),pointer :: pover
  type(T_assemble_block), pointer :: pS_neighbors
  do iatom, s%natoms
     matom = s%neigh_self(iatom)
     r1 = s%atom(iatom)%ratom
     in1 = s%atom(iatom)%imass
     num_neigh = s%neighbors(iatom)%neighn
     pover => nac_vars%gover(iatom)
     do ineigh = 1, num_neigh
        mbeta = s%neighbors(iatom)%neigh_b(ineigh)
        jatom = s%neighbors(iatom)%neigh_j(ineigh)
        r2 = s%atom(jatom)%ratom + s%xl(mbeta)%a
        in2 = s%atom(jatom)%imass
        pS_neighbors=>pover%neighbors(ineigh)

! Allocate the block size
        norb_nu = species(in2)%norb_max
        allocate (pS_neighbors%Dblock(3, norb_mu, norb_nu))
        pS_neighbors%Dblock = 0.0d0
            
! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
        z = distance (r1, r2)
        ! unit vector in sigma direction.
        if (z .lt. 1.0d-05) then
           sighat(1) = 0.0d0
           sighat(2) = 0.0d0
           sighat(3) = 1.0d0
        else
           sighat = (r2 - r1)/z
        end if
        call epsilon_function (r2, sighat, eps)
        call Depsilon_2c (r1, r2, eps, deps)
        
! Get the matrix from the data files - which is the matrix in molecular
! coordinates (stored in sm). Rotate the matrix into crystal coordinates.
! The rotated  matrix elements are stored in sx, where x means crytal
! coordinates.
! For these interactions, there are no subtypes and isorp = 0
        if (iatom .eq. jatom .and. mbeta .eq. 0) then
           call build_gover1c(in1)
           do inu = 1, norb_nu
              do imu = 1, norb_mu
                 pS_neighbors%Dblock = gover1c(:,imu,inu)
              end do
           end do
        else
           isorp = 0
           intercation = P_overlap
           in3 = in2
           

  end do ! iatom
