module pfaffian
  ! ==============================================================================
  ! Pfaffian for skew-symmetric matrices.
  !
  ! Algorithms:
  !   - skpfa_laplace / zskpfa_laplace: Laplace expansion along the first row,
  !     O(n!) work but bulletproof and sign-correct. Used for n <= ~12.
  !   - skpfa_reduction: Parlett-Reid tridiagonalization (Wimmer, arXiv:1102.3440)
  !     for general skew-symmetric A. Reduces A to skew-tridiagonal T via
  !     Householder reflections that preserve skew-symmetry, then evaluates
  !     Pf(T) = T(1,2) * T(3,4) * ... * T(n-1, n). Used for n > 12.
  !
  ! Skew-symmetrization guard: the input is symmetrized as (A - A^T)/2 before
  ! reduction; this guards against tiny roundoff that would otherwise cause
  ! the recurrence to produce a wrong sign.
  !
  ! Kitaev convention: Pf([[0, a], [-a, 0]]) = +a.
  !
  ! References:
  !   - M. Wimmer, arXiv:1102.3440 (2011).
  !   - A. Y. Kitaev, Phys.-Usp. 44, 131 (2001).
  ! ==============================================================================

  use definitions, only: dp
  implicit none
  private

  public :: real_pfaffian
  public :: complex_pfaffian
  public :: kitaev_majorana_number

contains

  ! ==============================================================================
  ! Real skew-symmetric Pfaffian. Dispatches to Laplace (n<=12) or Parlett-Reid.
  ! ==============================================================================
  function real_pfaffian(A) result(pf)
    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: pf
    real(kind=dp), allocatable :: Asym(:,:)
    integer :: n, i, j

    n = size(A, 1)
    if (n /= size(A, 2)) then
      pf = 0.0_dp
      return
    end if
    if (n == 0) then
      pf = 1.0_dp
      return
    end if
    if (mod(n, 2) /= 0) then
      pf = 0.0_dp
      return
    end if

    allocate(Asym(n, n))
    Asym = 0.0_dp
    do i = 1, n
      do j = i + 1, n
        Asym(i, j) = 0.5_dp * (A(i, j) - A(j, i))
        Asym(j, i) = -Asym(i, j)
      end do
    end do

    if (n <= 12) then
      pf = skpfa_laplace(Asym)
    else
      pf = skpfa_reduction(Asym, n)
    end if
    deallocate(Asym)
  end function real_pfaffian

  ! ==============================================================================
  ! Complex skew-symmetric Pfaffian.
  ! ==============================================================================
  function complex_pfaffian(A) result(pf)
    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: pf
    complex(kind=dp), allocatable :: Asym(:,:)
    integer :: n, i, j

    n = size(A, 1)
    if (n /= size(A, 2)) then
      pf = cmplx(0.0_dp, 0.0_dp, kind=dp)
      return
    end if
    if (n == 0) then
      pf = cmplx(1.0_dp, 0.0_dp, kind=dp)
      return
    end if
    if (mod(n, 2) /= 0) then
      pf = cmplx(0.0_dp, 0.0_dp, kind=dp)
      return
    end if

    allocate(Asym(n, n))
    Asym = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i = 1, n
      do j = i + 1, n
        Asym(i, j) = cmplx(0.5_dp, 0.0_dp, kind=dp) * (A(i, j) - A(j, i))
        Asym(j, i) = -Asym(i, j)
      end do
    end do

    if (n <= 12) then
      pf = zskpfa_laplace(Asym)
    else
      pf = zskpfa_reduction(Asym, n)
    end if
    deallocate(Asym)
  end function complex_pfaffian

  ! ==============================================================================
  ! Kitaev's Majorana number:
  !   M = sgn[Pf(H(k=0)*omega) * Pf(H(k=pi/a)*omega)]
  ! For class D BdG the antisymmetric structure matrix is omega = tau_y (x) I_n.
  ! Returns: +1 (trivial), -1 (topological), 0 (gap closure).
  ! ==============================================================================
  function kitaev_majorana_number(H_bdg, k_par_values, omega_struct) result(majorana_number)
    complex(kind=dp), intent(in) :: H_bdg(:,:)
    real(kind=dp), intent(in) :: k_par_values(:)
    complex(kind=dp), intent(in), optional :: omega_struct(:,:)
    integer :: majorana_number
    integer :: i, n_full
    real(kind=dp) :: prod
    complex(kind=dp), allocatable :: omega(:,:), A_work(:,:)
    complex(kind=dp) :: pf_val
    real(kind=dp) :: pf_abs, tol

    majorana_number = 0
    if (size(k_par_values) < 2) return
    n_full = size(H_bdg, 1)
    if (n_full /= size(H_bdg, 2) .or. mod(n_full, 2) /= 0) return

    if (present(omega_struct)) then
      allocate(omega(n_full, n_full))
      omega = omega_struct
    else
      call default_kitaev_omega(omega, n_full)
    end if

    tol = 1.0e-12_dp
    prod = 1.0_dp

    do i = 1, size(k_par_values)
      call assemble_skew_product(A_work, H_bdg, omega, n_full)
      pf_val = complex_pfaffian(A_work)
      pf_abs = real(sqrt(pf_val * conjg(pf_val)), kind=dp)
      if (pf_abs < tol) then
        prod = 0.0_dp
        if (allocated(A_work)) deallocate(A_work)
        exit
      end if
      prod = prod * real(pf_val, kind=dp)
      if (allocated(A_work)) deallocate(A_work)
    end do

    if (abs(prod) < tol) then
      majorana_number = 0
    else if (prod > 0.0_dp) then
      majorana_number = 1
    else
      majorana_number = -1
    end if

    if (allocated(omega)) deallocate(omega)
  end function kitaev_majorana_number

  subroutine default_kitaev_omega(omega, n_full)
    complex(kind=dp), allocatable, intent(out) :: omega(:,:)
    integer, intent(in) :: n_full
    integer :: half, i, j
    half = n_full / 2
    allocate(omega(n_full, n_full))
    omega = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i = 1, half
      do j = 1, half
        if (i == j) then
          omega(i, half + j) = cmplx(0.0_dp, -1.0_dp, kind=dp)
        else
          omega(i, half + j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        end if
        omega(half + j, i) = -omega(i, half + j)
      end do
    end do
  end subroutine default_kitaev_omega

  subroutine assemble_skew_product(A_work, H_bdg, omega, n_full)
    complex(kind=dp), allocatable, intent(out) :: A_work(:,:)
    complex(kind=dp), intent(in) :: H_bdg(:,:), omega(:,:)
    integer, intent(in) :: n_full
    complex(kind=dp), allocatable :: B(:,:)
    integer :: i, j

    allocate(B(n_full, n_full))
    B = matmul(H_bdg, omega)
    allocate(A_work(n_full, n_full))
    do i = 1, n_full
      do j = 1, n_full
        A_work(i, j) = cmplx(0.5_dp, 0.0_dp, kind=dp) * (B(i, j) - B(j, i))
      end do
    end do
    deallocate(B)
  end subroutine assemble_skew_product

  ! ==============================================================================
  ! Real Pfaffian via Laplace expansion (small n; O(n!)).
  ! Pf(A) = sum_{j=2..n} (-1)^j * A(1,j) * Pf(minor(A; 1, j))
  ! Convention: j=2 contributes +A(1,2)*Pf(minor).
  !   => Pf([[0,a],[-a,0]]) = +a (Kitaev).
  ! ==============================================================================
  recursive function skpfa_laplace(A) result(pf)
    real(kind=dp), intent(in) :: A(:,:)
    real(kind=dp) :: pf
    integer :: n, j, s
    real(kind=dp), allocatable :: sub(:,:)
    real(kind=dp) :: minor_val

    n = size(A, 1)
    if (n == 0) then
      pf = 1.0_dp
      return
    end if
    if (n == 2) then
      pf = A(1, 2)
      return
    end if

    pf = 0.0_dp
    s = 1
    do j = 2, n
      call extract_minor_sk(sub, A, 1, j)
      minor_val = skpfa_laplace(sub)
      deallocate(sub)
      if (s == 1) then
        pf = pf + A(1, j) * minor_val
      else
        pf = pf - A(1, j) * minor_val
      end if
      s = -s
    end do
  end function skpfa_laplace

  subroutine extract_minor_sk(sub, A, r, c)
    real(kind=dp), allocatable, intent(out) :: sub(:,:)
    real(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: r, c
    integer :: n, i_src, j_src, i_dst, j_dst

    n = size(A, 1)
    allocate(sub(n - 2, n - 2))
    i_dst = 0
    do i_src = 1, n
      if (i_src == r .or. i_src == c) cycle
      i_dst = i_dst + 1
      j_dst = 0
      do j_src = 1, n
        if (j_src == r .or. j_src == c) cycle
        j_dst = j_dst + 1
        sub(i_dst, j_dst) = A(i_src, j_src)
      end do
    end do
  end subroutine extract_minor_sk

  ! ==============================================================================
  ! Complex Pfaffian via Laplace expansion.
  ! ==============================================================================
  recursive function zskpfa_laplace(A) result(pf)
    complex(kind=dp), intent(in) :: A(:,:)
    complex(kind=dp) :: pf
    integer :: n, j, s
    complex(kind=dp), allocatable :: sub(:,:)
    complex(kind=dp) :: minor_val

    n = size(A, 1)
    if (n == 0) then
      pf = cmplx(1.0_dp, 0.0_dp, kind=dp)
      return
    end if
    if (n == 2) then
      pf = A(1, 2)
      return
    end if

    pf = cmplx(0.0_dp, 0.0_dp, kind=dp)
    s = 1
    do j = 2, n
      call extract_minor_zk(sub, A, 1, j)
      minor_val = zskpfa_laplace(sub)
      deallocate(sub)
      if (s == 1) then
        pf = pf + A(1, j) * minor_val
      else
        pf = pf - A(1, j) * minor_val
      end if
      s = -s
    end do
  end function zskpfa_laplace

  subroutine extract_minor_zk(sub, A, r, c)
    complex(kind=dp), allocatable, intent(out) :: sub(:,:)
    complex(kind=dp), intent(in) :: A(:,:)
    integer, intent(in) :: r, c
    integer :: n, i_src, j_src, i_dst, j_dst

    n = size(A, 1)
    allocate(sub(n - 2, n - 2))
    i_dst = 0
    do i_src = 1, n
      if (i_src == r .or. i_src == c) cycle
      i_dst = i_dst + 1
      j_dst = 0
      do j_src = 1, n
        if (j_src == r .or. j_src == c) cycle
        j_dst = j_dst + 1
        sub(i_dst, j_dst) = A(i_src, j_src)
      end do
    end do
  end subroutine extract_minor_zk

  ! ==============================================================================
  ! Real Parlett-Reid Pfaffian (Wimmer 2011, for n > 12).
  ! Reduces skew-symmetric A to skew-tridiagonal form, then evaluates the
  ! tridiagonal Pfaffian as T(1,2)*T(3,4)*...*T(n-1,n).
  ! ==============================================================================
  function skpfa_reduction(A_in, n) result(pf)
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: A_in(n, n)
    real(kind=dp) :: pf
    real(kind=dp), allocatable :: A(:,:)
    real(kind=dp) :: s, sn, sigma, alpha, tau, vkv
    real(kind=dp), allocatable :: v(:), w(:)
    integer :: k, i, j
    logical :: singular

    if (n == 0) then
      pf = 1.0_dp
      return
    end if
    if (n == 2) then
      pf = A_in(1, 2)
      return
    end if

    allocate(A(n, n))
    A = A_in
    pf = 1.0_dp
    singular = .false.

    k = 1
    do while (k <= n - 1)
      if (k + 2 > n) exit
      s = 0.0_dp
      do i = k + 2, n
        s = s + A(i, k) * A(i, k)
      end do
      if (s < tiny(1.0_dp)) then
        singular = .true.
        exit
      end if

      sn = sqrt(s)
      sigma = A(k + 1, k)
      if (sigma == 0.0_dp) then
        alpha = sn
        tau = 0.0_dp
      else
        alpha = -sign(sn, sigma)
      end if

      allocate(v(n))
      v = 0.0_dp
      v(k + 1) = A(k + 1, k) - alpha
      do i = k + 2, n
        v(i) = A(i, k)
      end do

      vkv = 0.0_dp
      do i = k + 1, n
        vkv = vkv + v(i) * v(i)
      end do
      if (vkv < tiny(1.0_dp)) then
        singular = .true.
        deallocate(v)
        exit
      end if
      tau = 2.0_dp / vkv

      allocate(w(n))
      w = 0.0_dp
      do i = 1, n
        do j = k + 1, n
          w(i) = w(i) + A(i, j) * v(j)
        end do
        w(i) = tau * w(i)
      end do

      do i = 1, n
        do j = 1, n
          A(i, j) = A(i, j) - v(i) * w(j) - w(i) * v(j)
        end do
      end do

      deallocate(v, w)
      k = k + 1
    end do

    if (singular) then
      pf = 0.0_dp
    else
      do i = 1, n - 1, 2
        pf = pf * A(i, i + 1)
      end do
    end if
    deallocate(A)
  end function skpfa_reduction

  ! ==============================================================================
  ! Complex Parlett-Reid Pfaffian.
  ! ==============================================================================
  function zskpfa_reduction(A_in, n) result(pf)
    integer, intent(in) :: n
    complex(kind=dp), intent(in) :: A_in(n, n)
    complex(kind=dp) :: pf
    complex(kind=dp), allocatable :: A(:,:)
    complex(kind=dp) :: s, sn, sigma, alpha, vkv, tau
    complex(kind=dp), allocatable :: v(:), w(:)
    real(kind=dp) :: s_norm
    integer :: k, i, j
    logical :: singular

    if (n == 0) then
      pf = cmplx(1.0_dp, 0.0_dp, kind=dp)
      return
    end if
    if (n == 2) then
      pf = A_in(1, 2)
      return
    end if

    allocate(A(n, n))
    A = A_in
    pf = cmplx(1.0_dp, 0.0_dp, kind=dp)
    singular = .false.

    k = 1
    do while (k <= n - 1)
      if (k + 2 > n) exit
      s_norm = 0.0_dp
      do i = k + 2, n
        s_norm = s_norm + real(A(i, k) * conjg(A(i, k)), kind=dp)
      end do
      if (s_norm < tiny(1.0_dp)) then
        singular = .true.
        exit
      end if

      s = cmplx(s_norm, 0.0_dp, kind=dp)
      sn = cmplx(sqrt(s_norm), 0.0_dp, kind=dp)
      sigma = A(k + 1, k)
      if (real(sigma * conjg(sigma), kind=dp) < tiny(1.0_dp)) then
        alpha = sn
        tau = cmplx(0.0_dp, 0.0_dp, kind=dp)
      else
        alpha = -sn * sigma / abs(sigma)
      end if

      allocate(v(n))
      v = cmplx(0.0_dp, 0.0_dp, kind=dp)
      v(k + 1) = sigma - alpha
      do i = k + 2, n
        v(i) = A(i, k)
      end do

      vkv = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do i = k + 1, n
        vkv = vkv + v(i) * v(i)
      end do
      if (real(vkv * conjg(vkv), kind=dp) < tiny(1.0_dp)) then
        singular = .true.
        deallocate(v)
        exit
      end if
      tau = cmplx(2.0_dp, 0.0_dp, kind=dp) / vkv

      allocate(w(n))
      w = cmplx(0.0_dp, 0.0_dp, kind=dp)
      do i = 1, n
        do j = k + 1, n
          w(i) = w(i) + A(i, j) * v(j)
        end do
        w(i) = tau * w(i)
      end do

      do i = 1, n
        do j = 1, n
          A(i, j) = A(i, j) - v(i) * w(j) - w(i) * v(j)
        end do
      end do

      deallocate(v, w)
      k = k + 1
    end do

    if (singular) then
      pf = cmplx(0.0_dp, 0.0_dp, kind=dp)
    else
      do i = 1, n - 1, 2
        pf = pf * A(i, i + 1)
      end do
    end if
    deallocate(A)
  end function zskpfa_reduction

end module pfaffian