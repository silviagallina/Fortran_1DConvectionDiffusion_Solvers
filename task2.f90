program task2

    use mod_param, only: rp, pi
    use mod_rk3, only: rk3
    use mod_save, only: save

    implicit none

    ! Parameters

    integer, parameter :: nx = 100
    real(rp), parameter :: L = 1.0_rp, c = 0.0_rp, v = 0.01_rp
    real(rp), parameter :: dt = 1.0e-03_rp, t_end = 3.0_rp
    real(rp), parameter :: x0 = L/2.0_rp, sigma = 0.1_rp

    ! Schemes

    character(len=*), parameter :: cscheme = 'none', dscheme = 'CS4'
    character(len=*), parameter :: label = '2CS4'

    ! Variables 

    real(rp) :: dx, t, D, dist
    real(rp), dimension(0:nx-1) :: x, phi, phi_exact
    integer :: i, n, nt, k

    dx = L/ real(nx, rp)
    do i=0, nx-1
        x(i) = i * dx
    end do

    D = 4*v*dt/dx**2
    write(*,'(A, F8.4)') 'D = ', D

    do i=0, nx-1
        phi(i) = exp(-(x(i)-x0)**2 / (2.0_rp * sigma**2))
    end do

    call save(label, 0, 0.0_rp, x, phi, phi, nx, L)

    nt = int(t_end / dt) 

    t = 0.0_rp

    do n=1, nt
        t = n * dt

        call rk3(phi, nx, dx, dt, c, v, cscheme, dscheme)

       do i = 0, nx-1
            phi_exact(i) = 0.0_rp

            ! La soluzione esatta deve considerare la periodicità ai bordi. Sto aggiungendo delle gaussiane "immaginarie" a distanza L, 2L, ecc. 
            ! dal centro x = 0.5 per tenere conto di questo effetto. 
            ! La larghezza effettiva della gaussiana al tempo t è sqrt(sigma**2 + 2*v*t) a causa della diffusione. A t = 3.0, questa larghezza è 
            ! circa 0.265. La regola pratica è che la gaussiana è "praticamente tutta" dentro ±3σ dal centro (x=0.5). Quindi la campana al tempo
            ! finale si estende da circa 0.5 - 3*0.265 = -0.295 (che è equivalente a 0.705 a causa della periodicità) a 
            ! 0.5 + 3*0.265 = 0.5 + 0.795 = 1.295 (che è equivalente a 0.295). La campana occupa circa 0.8 di larghezza totale (6σ_eff), che è 
            ! l'80% del dominio e le code escono di 0.3 da entrambi i lati.  
            ! Se la distanza tra la copia e il dominio è meno di ~3σ_eff, quella copia contribuisce, mentre se e è oltre 3σ_eff, posso ignorarla.
            ! Ad esempio, le copie a ±2L (0.5±2.0) sono a distanza 1.5, che è più di 3σ_eff, quindi non contribuiscono in modo significativo e 
            ! quindi mi fermo a k = ±2.
                do k = -2, 2
                    dist = x(i) - x0 - c*t + real(k, rp) * L
                    
                    phi_exact(i) = phi_exact(i) + (sigma / sqrt(sigma**2 + 2.0_rp*v*t)) * &
                    exp(-(dist**2) / (2.0_rp*(sigma**2 + 2.0_rp*v*t)))
                end do

        end do


        if (mod(n, 250) == 0) then
            call save(label, n, t, x, phi, phi_exact, nx, L)
        end if

    end do

    write(*, '(A)') 'Files saved: task2_t*.txt'

end program task2