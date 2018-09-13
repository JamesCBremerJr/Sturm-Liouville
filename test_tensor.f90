module test_tensor_functions

implicit double precision (a-h,o-z)

contains

subroutine testfun(x,y,val)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: x,y
double precision, intent(out) :: val

val = (1+y**2-y)*cos(x**2+y**2+1)

end subroutine

subroutine testfun2(x,y,val,derx,dery)
implicit double precision (a-h,o-z)
double precision, intent(in)  :: x,y
double precision, intent(out) :: val
val = y**2*x
derx = y**2
dery = 2*y*x
end subroutine

end module

program test_tensor

use utils
use chebyshev
use tensor
use test_tensor_functions

implicit double precision (a-h,o-z)

double precision, allocatable :: xs(:),ys(:),u(:,:)
double precision, allocatable :: vals(:,:),coefs(:)

!double precision, allocatable :: ab(:,:),vals2(:,:,:),coefs2(:,:)
integer*1, allocatable        :: compis(:),compjs(:)
double precision, allocatable :: compcoefs(:)

!
!  Construct the Chebyshev nodes and u matrix
!

n = 21
call chebpts(n,xs)
call elapsed(t1)
call tensor_umatrix(n,u)
call elapsed(t2)

call prin2("tensor_umatrix time = ",t2-t1)

!
!  Test an expansion
!

allocate(vals(n,n),coefs(n*n))

a =  0.11d0
b =  1.21d0

c = -0.23d0
d =  1.00d0

do i=1,n
x = xs(i) * (b-a)/2 + (b+a)/2
do j=1,n
y = xs(j) * (d-c)/2 + (c+d)/2
call testfun(x,y,val)
vals(i,j) = val
end do
end do

call tensor_coefs(n,u,vals,coefs)

!call prin2("vals = ",vals)
!call prin2("coefs = ",coefs)


x = 0.1d0
y = 0.5d0
call tensor_evalder(n,a,b,c,d,coefs,x,y,val,derx,dery)
call testfun(x,y,val0)

call prin2("evaluation error = ",val0-val)

stop


x = 0.412112d0
y = 0.53232d0
call tensor_eval(n,a,b,c,d,coefs,x,y,val)
call testfun(x,y,val0)

call prin2("evaluation error = ",val0-val)
call tensor_compressed_coefs(0.0d0,n,coefs,ncompcoefs,compcoefs)

call prini("ncoefs = ",n*(n+1)/2)
call prini("ncompcoefs = ",ncompcoefs)
call prin2("compcoefs = ",compcoefs)


val = 0
x = 0.412112d0
y = 0.53232d0
call testfun(x,y,val0)

call tensor_eval(n,a,b,c,d,coefs,x,y,val)
call  tensor_compressed_eval(compcoefs,a,b,c,d,x,y,val)

print *,""

print *,val
print *,val0
print *,(val-val0)/val0

stop


call tensor_compressed2_coefs(0.0d0,n,coefs,ncompcoefs,compcoefs,compis,compjs)

call tensor_compressed2_eval(eps,n,ncompcoefs,compcoefs,compis,compjs,a,b,c,d,x,y,val)
call prini("ncompcoefs = ",ncompcoefs)

x = 0.1d0
y = 0.5d0
call prin2("compressed evaluation error = ",val0-val)

stop

end program
