
! compile with: gfortran -O3 -fdefault-real-8 prob3_odd.f90 -o a && ./a

program simple_harmonic
    implicit none
    
    
    ! iterations for bisection cycle
    integer, parameter :: iterations = 1000
    real a, guess,x, dx, neg , fa, fb, n, E1,E2,E3,EE1(10),EE2;
    real initial_even,slope_even,initial_odd,slope_odd,initial,slope
    integer i, j,alpha
    real ,dimension(10) :: E1_range ,E2_range

    initial_even=1.0; slope_even=0.0;
    initial_odd=0.0; slope_odd=1.0;

    ! a = 2.0; guess=0.0;
    x=8.0; dx=0.01;
    neg=-1.0;
    n = floor(x/dx)

    ! E1=0.2;
    ! E2=0.8;

    ! do i=1,10
    ! E1_range(i)=0.2+1.0*(i-1);
    ! E2_range(i)=0.8+1.0*(i-1);
    ! end do

    EE1=[1.0,2.0,3.0,4.0,5.0,6.0,6.0,6.0,6.0,6.0];

    do j=1,10;
        if ( mod(j,2)==1 ) then 
            initial= initial_even; slope=slope_even;
            elseif ( mod(j,2)==0 ) then 
             initial= initial_odd; slope=slope_odd; 
        endif 
        !  E1=E1_range(j); E2=E2_range(j);
        
        if (j==1) then 
        E1=0.2; E2=0.8;
        else
        E1=E3+sqrt(EE1(j-1)-(EE1(j-1)/2.0))/2.0
        E2=E3+(EE1(j));
        end if
        ! write(*,*) j, E1, E2
        E3=bisection(E1,E2,initial,slope);
        write(*,*) E3
        ! fa=integrate(a, guess, neg*x, neg*dx,j,E3) 
        ! fb=integrate(a, guess, x, dx,(j+1),E3) 
         fa=integrate(initial, slope, neg*x, neg*dx,3,E3) 
         fb=integrate(initial, slope, x, dx,4,E3) 
    end do 

    contains
    
    !Bisection_function
    function bisection(a,b,initial,slope); !intent(in) a,b,E
        real a,b,E,c,fa,fb,fc, bisection ,initial,slope
            fa = f(initial,slope,a)
            fb = f(initial,slope,b)    
    
            
            if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."
    
            ! bisect the interval
        do i = 1,iterations
            c = (a+b)/2.0; fc = f(initial,slope,c); if (fc == 0.0) exit
            
            ! Ridder's variation on the basic method
            c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(initial,slope,c)
            
            ! write (*,*) c, fc
            
            if (fa*fc < 0.0) then; b = c; fb = fc; end if
            if (fc*fb < 0.0) then; a = c; fa = fc; end if
        end do
        bisection = c
        end function bisection


    ! potential
        pure function V(x); intent(in) x
        real V, x
        V = (1.0/4.0)*(x**4.0)
    end function
    
    
    subroutine evalf(y, dydx,E)
        real y(3), dydx(3) , E
        ! E=0.5;
        ! x and psi and psi_prime and solving schrodinger equation
        associate( x => y(1), y1 => y(2), dy1dx => y(3) )
        dydx(1) = 1.0 
        dydx(2) = dy1dx 
        dydx(3) = 2.0*(V(x)-E)*y1
    
        end associate
    end subroutine evalf
    
    
    
    ! 10th order implicit Gauss-Legendre integrator
    subroutine gl10(y, dx,E)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dx,E; integer i, k
        
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
                        call evalf(y + g(:,i)*dx, g(:,i),E)
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dx
    end subroutine gl10
    
    function integrate(y1, dy1dx, x, dx, alpha,E)
        real y1,x, dx, integrate ,dy1dx,E
        real y(3); integer i, n,j
        integer alpha
        ! character(len=1024) :: filename1
        ! character(len=1024) :: filename2

        ! filename1='data_q2_1_'//trim(j)//'.dat'
        ! filename2='data_q2_2_'//trim(j+1)//'.dat'


        ! ! start from zero at a given gradient
        y = [0.0, y1, dy1dx]; !E0 = E(u)
        
        ! number of time steps needed
        n = floor(x/dx)

        !  if (alpha==3) then
        ! open(unit=alpha,file="filename1")
        ! elseif (alpha==4)then
        ! open(unit=alpha,file="filename2")
        ! end if
        if (alpha==3) then
        open(unit=alpha,file="data_q2_1.dat")
        elseif (alpha==4)then
        open(unit=alpha,file="data_q2_2.dat")
        end if

        do i = 1,n
            call gl10(y, dx,E); if (abs(y(1)) > 100.0) exit
            write (alpha,'(4g24.16)') i*dx, y(1), y(2), y(3)
        end do
        close(alpha)
    

        call gl10(y,x-n*dx,E)
        
        ! return state at time t
        integrate = y(2)
    end function


    ! function to find a root of...
    function f(initial,slope,E); intent(in) E,initial,slope
        real f, E,initial,slope
        f = integrate(initial,slope, x, dx,2,E) 
    end function
    
    
    end program simple_harmonic
    
    
    
    
    
    
    
    