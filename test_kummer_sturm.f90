


program test_kummer_jacobi
use utils
use chebyshev
use odesolve
use kummer
use kummer_sturm

implicit double precision (a-h,o-z)
type(chebexps_data)             :: chebdata
type(sturm_phase)               :: phase
double precision, allocatable   :: avals(:,:),apvals(:,:),avals0(:),apvals0(:),dnus(:),ts(:)
double precision, allocatable   :: ab(:,:),dlambdas(:)
type(sturm_phase), allocatable  :: phases(:,:)


pi          = acos(-1.0d0)
eps         = 1.0d-13
k           = 30

call chebexps2(k,chebdata)

!
!  First build an expansion
!

call elapsed(t1)
call kummer_strum_expansion(chebdata,nints,ab,phases)
call elapsed(t2)
call prin2("expansion time = ",t2-t1)

call prini("nints = ",nints)
call prin2("ab = ",ab)

do int=1,nints
call  kummer_sturm_error(chebdata,nints,ab,phases,int,derr)
print *,int,derr
end do


call kummer_sturm_eigenvalues(chebdata,nints,ab,phases,ndlambda,dlambdas)

call prini("ndlambda = ",ndlambda)

end program
