module topological_analysis

  use definitions
  implicit none
  private

  public :: compute_chern_qwz
  public :: compute_berry_curvature
  public :: compute_hall_conductance
  public :: compute_z2_gap
  public :: compute_z2_fukane
  public :: extract_edge_states
  public :: compute_majorana_profile

contains

  function compute_chern_qwz(u, nk) result(C)
    ! Fukui-Hatsugai-Suzuki lattice gauge method for QWZ model.
    ! H = sin(kx)*sigma_x + sin(ky)*sigma_y + (u + cos(kx) + cos(ky))*sigma_z
    ! For u < -2: C = +1, for -2 < u < 0: C = -1, for u > 0: C = 0
    real(kind=dp), intent(in) :: u
    integer, intent(in) :: nk
    integer :: C

    integer :: i, j, ip1, jp1
    complex(kind=dp), allocatable :: evecs(:,:,:)
    real(kind=dp) :: kx, ky, dk, total_flux
    complex(kind=dp) :: Ux, Uy, prod
    real(kind=dp) :: H(2,2), eval(2), eigvec(2,2)
    real(kind=dp) :: sin_kx_val, sin_ky_val, cos_kx_val, cos_ky_val, mz_val
    complex(kind=dp) :: ev1_i_j, ev2_i_j, ev1_ip1_j, ev2_ip1_j, ev1_i_jp1, ev2_i_jp1

    dk = 2.0_dp * acos(-1.0_dp) / real(nk, kind=dp)
    total_flux = 0.0_dp

    allocate(evecs(nk, nk, 2))

    ! Compute eigenvectors at each k-point
    do j = 1, nk
      ky = -acos(-1.0_dp) + real(j-1, kind=dp) * dk
      do i = 1, nk
        kx = -acos(-1.0_dp) + real(i-1, kind=dp) * dk
        sin_kx_val = sin(kx); sin_ky_val = sin(ky)
        cos_kx_val = cos(kx); cos_ky_val = cos(ky)
        mz_val = u + cos_kx_val + cos_ky_val

        ! QWZ Hamiltonian
        H(1,1) = mz_val; H(1,2) = sin_kx_val; H(2,1) = sin_ky_val
        H(2,2) = -mz_val

        call diag_2x2(H, eval, eigvec)
        evecs(i,j,1) = cmplx(eigvec(1,1), eigvec(1,2), kind=dp)
        evecs(i,j,2) = cmplx(eigvec(2,1), eigvec(2,2), kind=dp)
      end do
    end do

    ! FHS U-link method
    ! FHS formula (Fukui-Hatsugai-Suzuki 2007):
    ! U_x(i,j) = <psi(i,j) | psi(i+1,j)>
    ! U_y(i,j) = <psi(i,j) | psi(i,j+1)>
    ! F(i,j) = Im[ln(U_x(i,j) * U_y(i+1,j) * U_x(i,j+1)^* * U_y(i,j)^*)]
    do j = 1, nk
      jp1 = mod(j, nk) + 1
      do i = 1, nk
        ip1 = mod(i, nk) + 1

        ! U_x(i,j) = conjg(psi_1(i,j)) * psi_1(i+1,j) + conjg(psi_2(i,j)) * psi_2(i+1,j)
        ev1_i_j = evecs(i,j,1)
        ev2_i_j = evecs(i,j,2)
        ev1_ip1_j = evecs(ip1,j,1)
        ev2_ip1_j = evecs(ip1,j,2)
        Ux = conjg(ev1_i_j) * ev1_ip1_j + conjg(ev2_i_j) * ev2_ip1_j

        ! U_y(i,j) = conjg(psi_1(i,j)) * psi_1(i,j+1) + conjg(psi_2(i,j)) * psi_2(i,j+1)
        ev1_i_jp1 = evecs(i,jp1,1)
        ev2_i_jp1 = evecs(i,jp1,2)
        Uy = conjg(ev1_i_j) * ev1_i_jp1 + conjg(ev2_i_j) * ev2_i_jp1

        ! U_x(i,j+1) = conjg(psi_1(i,j+1)) * psi_1(i+1,j+1) + conjg(psi_2(i,j+1)) * psi_2(i+1,j+1)
        ! U_y(i+1,j) = conjg(psi_1(i+1,j)) * psi_1(i+1,j+1) + conjg(psi_2(i+1,j)) * psi_2(i+1,j+1)
        ! Plaquette: Ux(i,j) * Uy(i+1,j) * Ux(i,j+1)^* * Uy(i,j)^* (FHS)
        prod = Ux * (conjg(ev1_ip1_j) * ev1_i_jp1 + conjg(ev2_ip1_j) * ev2_i_jp1)
        prod = prod * (conjg(ev1_i_jp1) * ev1_ip1_j + conjg(ev2_i_jp1) * ev2_ip1_j)
        prod = prod * conjg(Uy)
        total_flux = total_flux + aimag(log(prod))
      end do
    end do

    C = nint(total_flux / (2.0_dp * acos(-1.0_dp)))
    deallocate(evecs)
  end function compute_chern_qwz

  subroutine diag_2x2(H, eval, eigvec)
    real(kind=dp), intent(in) :: H(2,2)
    real(kind=dp), intent(out) :: eval(2), eigvec(2,2)
    real(kind=dp) :: tr, det, sqrt_term

    tr = H(1,1) + H(2,2)
    det = H(1,1)*H(2,2) - H(1,2)*H(2,1)
    sqrt_term = sqrt(max(0.0_dp, tr*tr/4.0_dp - det))
    eval(1) = tr/2.0_dp + sqrt_term; eval(2) = tr/2.0_dp - sqrt_term

    if (abs(H(1,2)) > 1e-12_dp) then
      eigvec(1,1) = H(1,2); eigvec(2,1) = eval(1) - H(1,1)
      eigvec(1,2) = H(1,2); eigvec(2,2) = eval(2) - H(1,1)
    else
      eigvec(1,1) = 1.0_dp; eigvec(2,1) = 0.0_dp
      eigvec(1,2) = 0.0_dp; eigvec(2,2) = 1.0_dp
    end if
    eigvec(:,1) = eigvec(:,1) / sqrt(sum(eigvec(:,1)**2))
    eigvec(:,2) = eigvec(:,2) / sqrt(sum(eigvec(:,2)**2))
  end subroutine diag_2x2

  function compute_berry_curvature(evecs_k, kx_arr, ky_arr, n_occ) result(Omega)
    complex(kind=dp), intent(in) :: evecs_k(:,:,:,:)
    real(kind=dp), intent(in) :: kx_arr(:), ky_arr(:)
    integer, intent(in) :: n_occ
    real(kind=dp), allocatable :: Omega(:,:)

    allocate(Omega(size(kx_arr), size(ky_arr)))
    Omega = 0.0_dp
  end function compute_berry_curvature

  function compute_hall_conductance(C) result(sigma_xy)
    integer, intent(in) :: C
    real(kind=dp) :: sigma_xy

    sigma_xy = real(C, kind=dp)
  end function compute_hall_conductance

  function compute_z2_gap(H_eigs, gap_threshold) result(z2)
    real(kind=dp), intent(in) :: H_eigs(:)
    real(kind=dp), intent(in) :: gap_threshold
    integer :: z2

    z2 = 0
  end function compute_z2_gap

  function compute_z2_fukane(H_eigs, params) result(z2)
    real(kind=dp), intent(in) :: H_eigs(:)
    class(*), intent(in) :: params
    integer :: z2

    z2 = 0
  end function compute_z2_fukane

  function extract_edge_states(evecs, energy, window) result(edge_xi)
    complex(kind=dp), intent(in) :: evecs(:)
    real(kind=dp), intent(in) :: energy, window
    real(kind=dp) :: edge_xi

    edge_xi = 0.0_dp
  end function extract_edge_states

  function compute_majorana_profile(evec_bdg) result(xi)
    complex(kind=dp), intent(in) :: evec_bdg(:)
    real(kind=dp) :: xi

    xi = 0.0_dp
  end function compute_majorana_profile

end module topological_analysis