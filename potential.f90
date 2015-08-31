module potential
  use types, only: dp
  use iotools, only: loadtxt
  use interp, only: cubic_coeff, cubic_eval
  implicit none 
  
contains
  subroutine phi_coeffs_from_grfile(grfile, rr, phi, bb, cc, dd)
    character(*), intent(in) :: grfile
    real(dp), dimension(:), allocatable, intent(out) :: rr, phi, bb, cc, dd

    real(dp), dimension(:, :), allocatable :: gr_raw
    integer :: npts
  
    call loadtxt(grfile, 2, gr_raw)
    npts = size(rr)
    allocate(rr(npts), phi(npts), bb(npts), cc(npts), dd(npts))
    rr = gr_raw(:, 1)
    phi = -log(gr_raw(:, 2))
    call cubic_coeff(rr, phi, bb, cc, dd, npts)
  end subroutine phi_coeffs_from_grfile
  
  real(dp) function spline_phi_at(r, rr, phi, bb, cc, dd, n) result(s)
    integer, intent(in) :: n
    real(dp), intent(in) :: r
    real(dp), dimension(n), intent(in) :: rr, phi, bb, cc, dd
    
    if (r.lt.rr(1)) then
       s = phi(1)*rr(1)/r ! Extrapolate assuming phi ~ 1/r
    else if (r.gt.rr(n)) then
       s = 0.d0
    else
       s = cubic_eval(r, rr, phi, bb, cc, dd, n)
    end if
  end function spline_phi_at
end module potential
    

