module mod_schemes

    use mod_param, only: rp

    implicit none

contains

    subroutine CS2(phi, rhs, nx, dx, c)

        integer, intent(in) :: nx
        real(rp), intent(in) :: phi(0:nx-1), dx, c
        real(rp), intent(out) :: rhs(0:nx-1)

        integer :: i, ip1, im1

        do i=0, nx-1
            if (i == 0) then
                ip1 = i + 1
                im1 = nx - 1
            else if (i == nx - 1) then
                ip1 = 0
                im1 = i - 1
            else
                ip1 = i + 1
                im1 = i - 1
            end if

            rhs(i) = - c * (phi(ip1)-phi(im1)) / (2.0_rp * dx)
        
        end do

    end subroutine CS2

    subroutine FW1(phi, rhs, nx, dx, c)

        integer, intent(in) :: nx
        real(rp), intent(in) :: phi(0:nx-1), dx, c
        real(rp), intent(out) :: rhs(0:nx-1)

        integer :: i, ip1

        do i=0, nx-1
            if (i == nx - 1) then
                ip1 = 0
            else
                ip1 = i + 1
            end if

            rhs(i) = - c * (phi(ip1)-phi(i)) / dx
        
        end do

    end subroutine FW1

    subroutine BW1(phi, rhs, nx, dx, c)

        integer, intent(in) :: nx
        real(rp), intent(in) :: phi(0:nx-1), dx, c
        real(rp), intent(out) :: rhs(0:nx-1)

        integer :: i, im1

        do i=0, nx-1
            if (i == 0) then
                im1 = nx - 1
            else
                im1 = i - 1
            end if

            rhs(i) = - c * (phi(i)-phi(im1)) / dx
        
        end do

    end subroutine BW1

    subroutine CS4(phi, rhs, nx, dx, v)

        integer, intent(in) :: nx
        real(rp), intent(in) :: phi(0:nx-1), dx, v
        real(rp), intent(out) :: rhs(0:nx-1)

        integer :: i, ip1, ip2, im1, im2

        do i=0, nx-1
            if (i == 0) then
                ip1 = i + 1
                ip2 = i + 2
                im1 = nx - 1
                im2 = nx - 2
            else if (i == 1) then
                ip1 = i + 1
                ip2 = i + 2
                im1 = i - 1
                im2 = nx - 1
            else if (i == nx - 2) then
                ip1 = i + 1
                ip2 = 0
                im1 = i - 1
                im2 = i - 2
            else if (i == nx - 1) then
                ip1 = 0
                ip2 = 1
                im1 = i - 1
                im2 = i - 2
            else 
                ip1 = i + 1
                ip2 = i + 2
                im1 = i - 1
                im2 = i - 2
            end if

            rhs(i) = v * (- phi(ip2) + 16.0_rp*phi(ip1) - 30.0_rp*phi(i) + 16.0_rp*phi(im1) - phi(im2)) / (12.0_rp * dx**2)
        
        end do

    end subroutine CS4

    subroutine computerhs(phi, rhs, nx, dx, c, v, cscheme, dscheme)

        integer, intent(in) :: nx
        real(rp), intent(in) :: phi(0:nx-1), dx, c, v
        real(rp), intent(out) :: rhs(0:nx-1)

        character(len=*), intent(in) :: cscheme, dscheme

        real(rp):: crhs(0:nx-1), drhs(0:nx-1)

        crhs = 0.0_rp
        drhs = 0.0_rp

        select case (cscheme)
            case ("CS2")
                call CS2(phi, crhs, nx, dx, c)
            case ("FW1")
                call FW1(phi, crhs, nx, dx, c)
            case ("BW1")
                call BW1(phi, crhs, nx, dx, c)
            case ("none")
                crhs = 0.0_rp
            case default
                write(*, '(A, A)') 'WARNING: unknown convection scheme: ', trim(cscheme)
        end select

        select case (dscheme)
            case ("CS4")
                call CS4(phi, drhs, nx, dx, v)
            case ("none")
                drhs = 0.0_rp
            case default
                write(*, '(A, A)') 'WARNING: unknown diffusion scheme: ', trim(dscheme)
        end select

        rhs = crhs + drhs

    end subroutine computerhs

end module mod_schemes