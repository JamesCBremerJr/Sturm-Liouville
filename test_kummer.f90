module test_kummer_subroutines

use utils
use chebyshev
use odesolve
use kummer
use iso_c_binding



implicit double precision (a-h,o-z)

type      qfun_data
double precision    :: dlambda,dgamma
end type  qfun_data


integer :: nintsbeta,k

double precision, allocatable :: xscheb(:),whtscheb(:),chebintl(:,:),chebintr(:,:), &
   ucheb(:,:),vcheb(:,:)

double precision, allocatable :: ab(:,:),alpha(:,:),alphap(:,:),alphapp(:,:),alphappp(:,:)
double precision, allocatable :: alphainv(:,:),alphainvp(:,:),abinv(:,:),xs(:,:)


double precision, allocatable :: abbeta(:,:),beta(:,:),betap(:,:),betapp(:,:),abin(:,:)

double precision, allocatable :: alpha_coefs(:,:),alphainv_coefs(:,:),alphap_coefs(:,:)
double precision, allocatable :: ts(:),vals(:),vals0(:),errs(:)
type(c_ptr)                   :: userptr
type(qfun_data), pointer      :: userdata


contains

subroutine qfun(t,val,userptr0)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: t
double precision, intent(out)  :: val
type(c_ptr)                    :: userptr0
type(qfun_data), pointer       :: userdata0

call c_f_pointer(userptr0,userdata0)
dlambda = userdata0%dlambda
dgamma  = userdata0%dgamma

val = dlambda**(2d0)  * sin(2*t)**2 / ( .01d0 + (t - 0.5d0)**4)
val = val + exp(-1/t)*dgamma**2

! val  = dlambda**(2d0)  * sin(2*t)**2 / ( .01d0 + (t - 0.5d0)**4) * exp(-1/t)

end subroutine




end module



program test_kummer

use utils
use odesolve
use kummer
use test_kummer_subroutines
use iso_c_binding

implicit double precision (a-h,o-z)

eps0 = epsilon(0.0d0)
pi   = acos(-1.0d0)


eps      = 1.0d-13
k        = 30
ifleft   = 0
dlambda  = 1000.0d0

call chebexps(k,xscheb,whtscheb,ucheb,vcheb,chebintl,chebintr)


!
!  Construct the phase function and its inverse
!

allocate(userdata)
userdata%dlambda = dlambda
userdata%dgamma  = dgamma

userptr          = c_loc(userdata)


call elapsed(t1)
call kummer_adap(eps,a,b,qfun,k,xscheb,chebintl,chebintr,ucheb, &
  nints,ab,alphap,alphapp,userptr)
call elapsed(t2)
t_phase = t2-t1


end program
