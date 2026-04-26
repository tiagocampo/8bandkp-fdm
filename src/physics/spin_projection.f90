module spin_projection
  use definitions
  implicit none
  private
  public :: spin_weights

contains

  ! Returns w_up = sum of |<orbital_up|psi>|^2 for a single eigenstate.
  ! w_dw = 1 - w_up.
  ! Basis ordering: HH_UP(1), LH_UP(2), LH_DW(3), HH_DW(4),
  !                 SO_UP(5), SO_DW(6), EL_UP(7), EL_DW(8)
  subroutine spin_weights(psi, Ngrid, w_up, w_dw)
    complex(kind=dp), intent(in) :: psi(:)  ! (8*Ngrid)
    integer, intent(in) :: Ngrid
    real(kind=dp), intent(out) :: w_up, w_dw

    integer :: n, b, idx
    complex(kind=dp) :: c(8)
    real(kind=dp) :: px_up, py_up, pz_up, ps_up

    w_up = 0.0_dp
    do n = 1, Ngrid
      do b = 1, 8
        idx = (b - 1) * Ngrid + n
        c(b) = psi(idx)
      end do

      ! X_up = |HH_UP>/sqrt(2) + |LH_DW>/sqrt(6) - i|SO_DW>/sqrt(3)
      px_up = abs(c(1)/sqrt(2.0_dp) + c(3)/sqrt(6.0_dp) &
        & - cmplx(0,1,dp)*c(6)/sqrt(3.0_dp))**2
      ! Y_up = i|HH_UP>/sqrt(2) - i|LH_DW>/sqrt(6) - |SO_DW>/sqrt(3)
      py_up = abs(cmplx(0,1,dp)*c(1)/sqrt(2.0_dp) &
        & - cmplx(0,1,dp)*c(3)/sqrt(6.0_dp) - c(6)/sqrt(3.0_dp))**2
      ! Z_up = -i*sqrt(2)*|LH_UP>/sqrt(3) + |SO_UP>/sqrt(3)
      pz_up = abs(-cmplx(0,1,dp)*sqrt(2.0_dp)*c(2)/sqrt(3.0_dp) &
        & + c(5)/sqrt(3.0_dp))**2
      ! S_up = |EL_UP>
      ps_up = abs(c(7))**2

      w_up = w_up + px_up + py_up + pz_up + ps_up
    end do
    w_dw = 1.0_dp - w_up
  end subroutine spin_weights

end module spin_projection
