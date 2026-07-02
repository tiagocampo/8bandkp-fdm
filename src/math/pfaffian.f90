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
  public :: kitaev_bdg_fixture_2band
  public :: polar_decomposition_sign

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
  ! Kitaev Majorana number via Lutchyn-Oreg sign-of-det (ADR 0008 section 1):
  !   M = sign( prod_k det(U_odd(k)) ), where U = sign(H_bdg(k)).
  ! =========================-=================================================
  function kitaev_majorana_number(H_k_array, k_par_values, omega_struct) result(majorana_number)
    complex(kind=dp), intent(in) :: H_k_array(:,:,:)
    real(kind=dp), intent(in) :: k_par_values(:)
    complex(kind=dp), intent(in), optional :: omega_struct(:,:)
    integer :: majorana_number
    integer :: i, n_full, n_k, n_odd, info
    real(kind=dp) :: sign_prod
    complex(kind=dp), allocatable :: U_odd(:,:)
    complex(kind=dp) :: det_val

    majorana_number = 0
    n_k = size(H_k_array, 3)
    if (n_k < 2) return
    if (size(k_par_values) < 2) return
    n_full = size(H_k_array, 1)
    if (n_full /= size(H_k_array, 2) .or. mod(n_full, 2) /= 0) return
    n_odd = n_full / 4
    if (n_odd < 1) return

    allocate(U_odd(n_odd, n_odd))
    sign_prod = 1.0_dp

    do i = 1, n_k
      call polar_decomposition_sign(H_k_array(:,:,i), U_odd, info)
      if (info == 1) then
        error stop 'kitaev_majorana_number: H_bdg not Hermitian'
      else if (info == 2) then
        error stop 'kitaev_majorana_number: zheev eigendecomposition failed'
      end if
      ! Compute det(U_odd) directly (small matrix).
      if (n_odd == 1) then
        det_val = U_odd(1, 1)
      else
        ! 2x2 det for n_odd=2 (typical case).
        det_val = U_odd(1,1) * U_odd(2,2) - U_odd(1,2) * U_odd(2,1)
      end if
      if (abs(det_val) < 1.0e-12_dp) then
        majorana_number = 0
        deallocate(U_odd)
        return
      end if
      sign_prod = sign_prod * sign(1.0_dp, real(det_val, kind=dp))
    end do

    if (sign_prod > 0.0_dp) then
      majorana_number = 1
    else
      majorana_number = -1
    end if

    if (allocated(U_odd)) deallocate(U_odd)
  end function kitaev_majorana_number

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
      s = A(k + 1, k) * A(k + 1, k)
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
      end do
      do i = 1, n
        do j = 1, n
          A(i, j) = A(i, j) + tau * (v(i) * w(j) - w(i) * v(j))
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
      s_norm = real(A(k + 1, k) * conjg(A(k + 1, k)), kind=dp)
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
      end do
      do i = 1, n
        do j = 1, n
          A(i, j) = A(i, j) + tau * (v(i) * w(j) - w(i) * v(j))
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

  ! ==============================================================================
  ! 2-band spinless p-wave Kitaev chain fixture at PHS-invariant momentum k.
  ! Nambu basis (c, c^dagger) block-diagonal: two identical 2x2 blocks.
  ! Each block: [[eps - mu, delta_0], [delta_0, -(eps - mu)]] with eps = -2 t cos k.
  ! Gap closes at |mu| = 2t.
  ! PHS: C H(k) C^{-1} = -H(-k) with C = tau_x K.
  ! ==============================================================================
  function kitaev_bdg_fixture_2band(mu, t, delta_0, k) result(H)
    real(kind=dp), intent(in) :: mu, t, delta_0, k
    complex(kind=dp) :: H(4, 4)
    real(kind=dp) :: eps
    eps = -2.0_dp * t * cos(k)
    H = cmplx(0.0_dp, 0.0_dp, kind=dp)
    H(1,1) = cmplx(eps - mu, 0.0_dp, kind=dp)
    H(2,2) = cmplx(-(eps - mu), 0.0_dp, kind=dp)
    H(1,2) = cmplx(delta_0, 0.0_dp, kind=dp)
    H(2,1) = cmplx(delta_0, 0.0_dp, kind=dp)
    H(3,3) = cmplx(eps - mu, 0.0_dp, kind=dp)
    H(4,4) = cmplx(-(eps - mu), 0.0_dp, kind=dp)
    H(3,4) = cmplx(delta_0, 0.0_dp, kind=dp)
    H(4,3) = cmplx(delta_0, 0.0_dp, kind=dp)
  end function kitaev_bdg_fixture_2band

  ! ==============================================================================
  ! Polar-decomposition sign (Lutchyn-Oreg majority-rule algorithm per ADR 0008
  ! section 1). For Hermitian H_bdg(k), sign(H) = sum_j sign(eigval_j) |v_j><v_j|
  ! where {eigval_j, v_j} are eigendecomposition pairs (zheev). The odd-parity
  ! Nambu subspace is rows/cols 1..n_sp/2, so U_odd is the restricted sign(H).
  ! det(U_odd) is purely real (since sign(H) is involutory and Hermitian) and
  ! its sign maps onto the BdG Z2 invariant.
  !
  ! info = 0: success; U_odd is n_odd x n_odd.
  ! info = 1: H_bdg shape wrong or not Hermitian.
  ! info = 2: LAPACK zheev failed.
  ! ==============================================================================
  subroutine polar_decomposition_sign(H_bdg, U_odd, info)
    use linalg, only: zheev
    complex(kind=dp), intent(in) :: H_bdg(:,:)
    complex(kind=dp), allocatable, intent(out) :: U_odd(:,:)
    integer, intent(out) :: info
    integer :: n_full, n_sp, n_odd, lwork, lrwork, info_lapack, i, j, k
    complex(kind=dp), allocatable :: work(:)
    real(kind=dp), allocatable :: rwork(:)
    complex(kind=dp), allocatable :: eigenvectors(:,:)
    real(kind=dp), allocatable :: eigvals(:)

    info = 0
    n_full = size(H_bdg, 1)
    if (n_full /= size(H_bdg, 2)) then; info = 1; return; end if
    if (mod(n_full, 2) /= 0) then; info = 1; return; end if
    n_sp = n_full / 2
    n_odd = n_sp / 2
    if (n_odd < 1) then; info = 1; return; end if

    ! Hermiticity check (strict).
    do i = 1, n_full
      do j = i + 1, n_full
        if (abs(H_bdg(i,j) - conjg(H_bdg(j,i))) > 1.0e-12_dp) then
          info = 1
          return
        end if
      end do
    end do

    ! Eigendecomposition via zheev.
    allocate(eigvals(n_full), eigenvectors(n_full, n_full))
    eigenvectors = H_bdg
    lwork = max(1, 2*n_full - 1)
    lrwork = max(1, 3*n_full - 2)
    allocate(work(max(lwork, n_full*n_full)), rwork(lrwork))
    call zheev('V', 'U', n_full, eigenvectors, n_full, eigvals, work, size(work), &
               rwork, info_lapack)
    if (info_lapack /= 0) then
      info = 2
      deallocate(work, rwork, eigvals, eigenvectors)
      return
    end if

    ! sign(H) restricted to odd-parity block.
    allocate(U_odd(n_odd, n_odd))
    U_odd = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do i = 1, n_odd
      do j = 1, n_odd
        do k = 1, n_full
          U_odd(i,j) = U_odd(i,j) + sign(1.0_dp, eigvals(k)) * &
                       eigenvectors(i, k) * conjg(eigenvectors(j, k))
        end do
      end do
    end do

    deallocate(work, rwork, eigvals, eigenvectors)
  end subroutine polar_decomposition_sign

end module pfaffian
