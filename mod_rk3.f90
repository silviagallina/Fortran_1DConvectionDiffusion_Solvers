module mod_rk3

    use mod_param, only: rp
    use mod_schemes, only: computerhs

    implicit none

    contains 

        subroutine rk3(phi, nx, dx, dt, c, v, cscheme, dscheme)

            integer, intent(in) :: nx
            real(rp), intent(in) :: dx, dt, c, v
            real(rp), intent(inout) :: phi(0:nx-1)

            character(len=*), intent(in) :: cscheme, dscheme

            real(rp) :: rhs1(0:nx-1), rhs2(0:nx-1), rhs3(0:nx-1), phi1(0:nx-1), phi2(0:nx-1)

            ! Step 1

            call computerhs(phi, rhs1, nx, dx, c, v, cscheme, dscheme)
            phi1 = phi + 1.0_rp / 3.0_rp * dt * rhs1

            ! Step 2

            call computerhs(phi1, rhs2, nx, dx, c, v, cscheme, dscheme)
            phi2 = phi + 2.0_rp / 3.0_rp * dt * rhs2

            ! Step 3
            call computerhs(phi2, rhs3, nx, dx, c, v, cscheme, dscheme)

            ! Final update

            phi = phi + dt * (1.0_rp / 4.0_rp * rhs1 + 3.0_rp / 4.0_rp * rhs3)

        end subroutine rk3

end module mod_rk3