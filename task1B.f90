program task1B

    use mod_param, only: rp, pi
    use mod_rk3, only: rk3
    use mod_save, only: save

    implicit none

    ! Parameters

    integer, parameter :: nx = 50
    real(rp), parameter :: L = 1.0_rp, c = 0.1_rp, v = 0.0_rp
    real(rp), parameter :: dt = 0.25_rp, t_end = 100.0_rp
    real(rp), parameter :: x0 = L/2.0_rp, sigma = 0.1_rp

    ! Schemes

    character(len=*), parameter :: cscheme = 'CS2', dscheme = 'none'
    character(len=*), parameter :: label = '1BCS2'

    ! Variables 

    real(rp) :: dx, t, CFL, dist
    real(rp), dimension(0:nx-1) :: x, phi, phi_exact
    integer :: i, n, nt

    dx = L/ real(nx, rp)
    do i=0, nx-1
        x(i) = i * dx
    end do

    CFL = c * dt / dx
    write(*,'(A, F8.4)') 'CFL = ', CFL

    do i=0, nx-1
        phi(i) = sin(2 * pi * x(i) / L)
    end do

    call save(label, 0, 0.0_rp, x, phi, phi, nx, L)

    nt = int(t_end / dt) 

    t = 0.0_rp

    do n=1, nt
        t = n * dt

        call rk3(phi, nx, dx, dt, c, v, cscheme, dscheme)

        do i=0, nx-1
            phi_exact(i) = exp(-v*(2*pi/L)**2*t)*sin(2*pi/L * (x(i) -c*t))
        end do

        if (n <= 50 .and. mod(n, 5) == 0) then
            call save(label, n, t, x, phi, phi_exact, nx, L)
        end if

    end do

    write(*, '(A)') 'Files saved: task1B_t*.txt'

end program task1B