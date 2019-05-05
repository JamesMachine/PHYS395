
! compile with: gfortran -O3 -fdefault-real-8 problem1.f90 -o a && ./a

program simple_harmonic
    implicit none
    
    
    ! iterations for bisection cycle
    integer, parameter :: iterations = 1000
    real a, guess,x, dx, neg , fa, fb, n
    integer i ,alpha

    a = 0.0; guess=1.0;
    x=5.0; dx=0.01;
    neg=-1.0;
    n = floor(x/dx)


    fa=integrate(a, guess, neg*x, neg*dx,1) 
    fb=integrate(a, guess, x, dx,2) 


    contains
    
    ! potential
    pure function V(x); intent(in) x
        real V, x
        V = 0.5*(x*x)
        ! V = (1.0/4.0)*(x**4.0)
    end function
    
    
    subroutine evalf(y, dydx)
        real y(3), dydx(3) , E
        E=  1.5;
        ! x and psi and psi_prime and solving schrodinger equation
        associate( x => y(1), y1 => y(2), dy1dx => y(3) )
        dydx(1) = 1.0 
        dydx(2) = dy1dx 
        dydx(3) = 2.0*(V(x)-E)*y1
    
        end associate
    end subroutine evalf
    
    
    
    ! 10th order implicit Gauss-Legendre integrator
    subroutine gl10(y, dx)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dx; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dx, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dx
    end subroutine gl10
    
    function integrate(y1, dy1dx, x, dx, alpha)
        real y1,x, dx, integrate ,dy1dx
        real y(3); integer i, n
        integer alpha

        ! ! start from zero at a given gradient
        y = [0.0, y1, dy1dx]; !E0 = E(u)
        
        ! number of time steps needed
        n = floor(x/dx)
        if (alpha==1) then
        open(unit=alpha,file="data_q1_1.dat")
        elseif (alpha==2)then
        open(unit=alpha,file="data_q1_2.dat")
        end if

        do i = 1,n
            call gl10(y, dx); if (abs(y(1)) > 100.0) exit
            write (alpha,'(4g24.16)') i*dx, y(1), y(2), y(3)
        end do
        close(alpha)
    

        call gl10(y,x-n*dx)
        
        ! return state at time t
        integrate = y(3)
    end function
    
    
    end program simple_harmonic
    
    
    
    
    
    
    
    