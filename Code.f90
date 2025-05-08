program slinky

implicit none
double precision, allocatable :: x(:),y(:),m(:),ay(:),vy(:),ae(:),l(:),ax(:),vx(:)
double precision::g,dt,l0,tmax,t,c,k,tol,y0,omega,m_tot,y_com,x_com,length,tc,q,E,ke,rij
double precision::del_x,d
integer ::N,i,Nmax,j

!Nmax = 100  !for Section 2.1.1
N = 10

!do while(N <= Nmax)  !for Section 2.1.1
allocate(x(0:N),y(0:N),m(0:N),vx(0:N),ax(0:N),ay(0:N),vy(0:N),ae(0:N),l(0:N-1))

!Parameters
dt = 1.0d-4
g = 9.8d0
l0 = 0.0d0
tmax = 1.0d0
c = 5.d-3  !the other sections
d = 1.0d0  !for Section 2.1.2
del_x = d/N  !for Section 2.1.2
!c = 1.0d-2  !for Section 2.1.2
tol = 1.0d-6
omega = 30.0d0
ke = 9.0d9
q = 1.0d-12
E = 9.0d9

open(unit = 7, file = "results.txt")

!Initialization
i = 0
t = 0.0d0
y_com = 0.0d0
x_com = 0.0d0
k = 20.0d0
!k=9800    !true catenary curve
!m(0:N)=1.27d-4  !true catenary curve
m(0:N) = 1.0d-3
l(0:N-1) = 0.0d0
vy(0:N) = 0.0d0
vx(0:N) = 0.0d0
y(0:N) = 0.0d0
x(0:N) = 0.0d0
ay(N) = -9.8d0 
!ay(0:N) = 0.0d0  !for Section 2.1.2
ax(0:N) = 0.0d0
y_com = 0.0d0
m_tot = ((N+1)*m(1))

!!initial state of the slinky
!do i = 0,N                        !Section 2.1.2:
!x(i) = i*del_x
!end do


!do i=0,N "task 2.1.2 variable masses"
!if (mod(i, 2) == 1) then
!m(i)=2.0d-3
!else                       !Section 2.1.2:
!m(i)=1.0d-3
!end if
!end do 


!optional task (section 2.1.2) 
!y(0:N-1) = 0.0d0      
!y(N)=0.05d0
!x(N)=1.5d0
!do i = 0,N
!y_com = y_com + y(i)*m(i)   
!x_com = x_com + x(i)*m(i)
!end do

!y_com = y_com/m_tot
!x_com = x_com/m_tot


!do while(x(0) <= x(N))  !optional task codition (Section 2.1.2)

!do while (t <= tmax)  !for Section 2.1.2

do while (abs(ay(N)) >= tol) !for Section 2.1.1

y(1:N) = y(1:N) + vy(1:N)*dt  !for Section 2.1.1

!x(1:N-1) = x(1:N-1) + vx(1:N-1)*dt  !for Section 2.1.2
!y(1:N-1) = y(1:N-1) + vy(1:N-1)*dt  !for Section 2.1.2

!calculating Li (Section 2.1.2):
!do i = 0,N-1
!l(i) = sqrt((y(i+1) - y(i))**2 + (x(i+1) - x(i))**2)
!end do

do i = 0,N
ay(i) = 0.0d0

if(i == N)then
ay(i) = k*((y(i-1) - y(i))-l0) - m(i)*g - c*vy(i)
ay(i) = ay(i)/m(i)

!for section 2.1.2:
!ax(i) = k*((l(N-1)-l0)/l(N-1))*(x(i-1)-x(i)) - c*vx(i)
!ax(i) = ax(i)/m(i)
!ay(i) = k*((l(N-1)-l0)/l(N-1))*(y(i-1)-y(i)) - c*vy(i) - m(i)*g
!ay(i) = ay(i)/m(i)

else if (i == 0)then
ay(i) = k*((y(i+1) - y(i))-l0) - m(i)*g - c*vy(i)
ay(i) = ay(i)/m(i)

!for section 2.1.2:
!ax(i) = k*((l(i)-l0)/l(i))*(x(i+1)-x(i)) - c*vx(i)
!ax(i) = ax(i)/m(i)
!ay(i) = k*((l(i)-l0)/l(i))*(y(i+1)-y(i)) - c*vy(i) - m(i)*g
!ay(i) = ay(i)/m(i)

else
ay(i) = k*((y(i-1) - y(i))-l0) + k*((y(i+1) - y(i))-l0) - m(i)*g - c*vy(i)
ay(i) = ay(i)/m(i)

!for section 2.1.2:
!ax(i) = k*((l(i-1)-l0)/l(i-1))*(x(i-1)-x(i)) + k*((l(i)-l0)/l(i))*(x(i+1)-x(i)) - c*vx(i)
!ax(i) = ax(i)/m(i)
!ay(i) = k*((l(i-1)-l0)/l(i-1))*(y(i-1)-y(i)) + k*((l(i)-l0)/l(i))*(y(i+1)-y(i)) - c*vy(i) - m(i)*g
!ay(i) = ay(i)/m(i)
end if
end do

!updating velocity (Section 2.1.2):
!vx(1:N-1) = vx(1:N-1) + ax(1:N-1)*dt
!vy(1:N-1) = vy(1:N-1) + ay(1:N-1)*dt

!write(7,*)x(0:N),y(0:N)     "Optional task (section 2.1.2)"

!The center of mass
!y_com = 0.0d0
!x_com = 0.0d0
!do i = 0,N
!y_com = y_com + y(i)*m(i)   !"Optional task (section 2.1.2)"
!x_com = x_com + x(i)*m(i)
!end do
!y_com = y_com/m_tot
!x_com = x_com/m_tot

vy(1:N) = vy(1:N) + ay(1:N)*dt

t = t + dt

!write(7,*)x(0:N),y(0:N),x_com,y_com
end do

!for task 2.1.2:
!do i=0,N
!write(7,*)(x(i)-x(N)/2),(y(i)+1) !(shifting the data for catenary fitting)
!end do 
!write(6,*)(y(N/2)+1) !to get the catenary constant a 



!Initializations for Section 2.3:
length = abs(y(N) - y(0))
t = 0.0d0
c = 0.0d0
y_com = 0.0d0
tc = sqrt((2.0d0*length)/(3.0d0*g))  !From theory (time for mN to move)
ay(0:N) = 0.0d0
vy(0:N) = 0.0d0


!The center of mass of the slinky at t = 0
do i = 0,N
y_com = y_com + y(i)*m(i)
end do
y_com = y_com/m_tot


do while(abs(vy(N)) <= 1.20d-4)  !for Section 2.3.1

!do while(t <= tmax)   !for Section 2.3.2

!update the position
y(0:N) = y(0:N) + vy(0:N)*dt  !for Section 2.3.1 & for Section 2.3.2 (charged slinky)

!y(1:N) = y(1:N) + vy(1:N)*dt   !for Section 2.3.2 (osscilating v0)

!do i = 0,N
!ae(i) = 0.0d0  !for Section 2.3.2
!ae(i) = 0.0d0

!do j = 0,N
!if (j == i) cycle
!rij = y(j)-y(i)
!if (i > j)ae(i) = ae(i) + (ke*q**2)/rij**2  !for Section 2.3.2
!if (i < j)ae(i) = ae(i) - (ke*q**2)/rij**2
!end do
!end do

!update oscillate i particle costantly
!y(0) = (length*0.2d0)*cos(omega*t) -(length*0.2)  !for Section 2.3.2

do i = 0,N
ay(i) = 0.0d0

if(i == N)then
ay(i) = k*((y(i-1) - y(i))-l0) - m(i)*g - c*vy(i) !+ ae(i) + q*E  !for Section 2.3.2
ay(i) = ay(i)/m(i)

else if (i == 0)then
ay(i) = k*((y(i+1) - y(i))-l0) - m(i)*g - c*vy(i) !+ ae(i) + q*E  !for Section 2.3.2
ay(i) = ay(i)/m(i)

else
ay(i) = k*((y(i-1) - y(i))-l0) + k*((y(i+1) - y(i))-l0) - m(i)*g - c*vy(i) !+ ae(i) + q*E  
ay(i) = ay(i)/m(i)
end if
end do

vy(0:N) = vy(0:N) + ay(0:N)*dt  !for Section 2.3.1

!vy(1:N) = vy(1:N) + ay(1:N)*dt  !for Section 2.3.2

t = t + dt

!The center of mass
y_com = 0.0d0
do i = 0,N
y_com = y_com + y(i)*m(i)
end do
y_com = y_com/m_tot

write(7,*)t,0.0d0,y(0:N),y_com

end do

write (6,*)"t = ",t,"tc = ",tc


close(7)
deallocate(x,y,m,ay,vy,ae,l,ax,vx)
end program slinky

!!!!!!!!!!
!check if m is total m of slinky or m of particle in the article for com