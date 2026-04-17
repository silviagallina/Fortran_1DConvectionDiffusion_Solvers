program task1A

    use mod_param, only: rp
    use mod_rk3, only: rk3
    use mod_save, only: save

    implicit none

    ! Parameters

    integer, parameter :: nx = 100
    real(rp), parameter :: L = 1.0_rp, c = 0.1_rp, v = 0.0_rp
    real(rp), parameter :: dt = 1.0e-2_rp, t_end = 1.0_rp
    real(rp), parameter :: x0 = L/2.0_rp, sigma = 0.1_rp

    ! Schemes

    character(len=*), parameter :: cscheme = 'CS2', dscheme = 'none'
    character(len=*), parameter :: label = 'task1A'

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
        phi(i) = exp(-((x(i) - x0)**2 / (2.0_rp * sigma**2)))
    end do

    call save(label, 0, 0.0_rp, x, phi, phi, nx, L)

    nt = int(t_end / dt) 

    t = 0.0_rp

    do n=1, nt
        t = n * dt

        call rk3(phi, nx, dx, dt, c, v, cscheme, dscheme)

        do i=0, nx-1
            
            dist = x(i) - x0 - c*t
            
            dist = mod(dist, L)
            
            if (dist > L / 2.0_rp) then
                dist = dist - L
            else if (dist < -L / 2.0_rp) then
                dist = dist + L
            end if

            phi_exact(i) = (sigma / sqrt(sigma**2 + 2.0_rp * v * t)) * &
                           exp(-(dist**2) / (2.0_rp * (sigma**2 + 2.0_rp * v * t)))
        end do

        if (mod(n, 25) == 0) then
            call save(label, n, t, x, phi, phi_exact, nx, L)
        end if

    end do

    write(*, '(A)') 'Files saved: task1A_t*.txt'

end program task1A