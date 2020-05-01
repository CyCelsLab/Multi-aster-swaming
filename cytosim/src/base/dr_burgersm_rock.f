c * * * * * * * * * * * * * * * * * * * * * * * * *
c    Driver for ROCK4 (or ROCK2) at Burger problem
c * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use ROCK4 (or ROCK2). It solves a
c    system of ODEs resulting from discretization of the 
c    Burger equation (u=u(x,y)):
c     
c    u_t+(u^2/2)_x=m*u_{xx}  for 0<= t <=2.5, 0<= x <= 1
c                            where m=0.0003
c
c    with initial condition 
c
c    u(x,0)=1.5*x*(1-x)^2
c
c    and  boundary conditions
c
c    u(0,t)=u(1,t)=0 
c
c
c    We discretize the space variables with
c    x_i=i/(N+1) for i=0,1,...,N, with N=500.
c    We obtain a system of 500 equations. 
c    The spectral radius of the Jacobian can be  
c    estimated with the Gerhgorin theorem. Thus  
c    we provide an external function RHO, 
c    giving the spectral radius of the Jacobian 
c    matrix. As output point we choose t_out=2.5
c 
c
c--------------------------------------------------------
c ----- to integrate with rock2.f ----- 
c     include 'rock2.f'  
c
      subroutine cytosimf(dt,neqn,y,work,tol)   
      implicit double precision (a-h,o-z)
      integer neqn
      integer iwork(12)
      double precision dt,tol
      dimension y(neqn),work(5*neqn)
      external fburg
c --- common parameters for the problem -----
      common/trans/alf,ns,nssq,nsnsm1,nsm1sq
      
      
c--------------------------------------------------------
c     Initialise iwork: 
c      iwork(1)=1  RHO returns an upper bound for the spectral radius
c      iwork(2)=1  The Jacobian is constant (RHO is called once)
c      iwork(3)=0  Return and solution at tend
c      iwork(4)=0  Atol and rtol are scalar
c--------------------------------------------------------
      iwork(1)=0
      iwork(2)=0
      iwork(3)=0
      iwork(4)=0
c ----- initial and end point of integration -----
      t=0.0d0
      tend=dt
c ----- required tolerance -----
      rtol=tol
      atol=tol
c ----- initial step size -----
      h=dt
c ----- integration -----
c      write(6,*) 'Integration of the Burger problem'     
c ----- to integrate with rock2.f ----- 
      call rock2(neqn,t,tend,h,y,fburg,atol,rtol,work,
     &           iwork,idid) 
c
c ----- call of the subroutine rock4 -----
c       call rock4(neqn,t,tend,h,y,fburg,atol,rtol,work,
c     &           iwork,idid)
c ----- print solution -----
c      do j=1,500
c        write (6,*) y(j)
c      end do
c ----- print statistics -----
c      write(6,*) 'Solution is tabuled in file sol.out'
c      write(6,*) 'The value of IDID is',idid
      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write(6,91) iwork(5),iwork(6),iwork(7),iwork(8)
 91   format(' Number of f evaluation=',i5,' steps=',i4,
     &        ' accpt=',i4,' rejct=',i3)
c--------------------------------------------------------
c     End of main program
c--------------------------------------------------------
      end      
c--------------------------------------------------------
c     The subroutine RHO gives an estimation of the spectral 
c     radius of the Jacobian matrix of the problem. This
c     is a bound for the whole interval and thus RHO is called
c     once.
c--------------------------------------------------------
      double precision function rho(neqn,t,y)
      implicit double precision (a-h,o-z)
      rho=2400000.d0
      return
      end
c--------------------------------------------------------
c     The subroutine FBURG compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fburg(neqn,x,y,f)
      implicit double precision (a-h,o-z)
      dimension x(neqn),y(neqn),f(neqn)
      external fburgc
      call fburgc(y,f)
      return
      end
c
