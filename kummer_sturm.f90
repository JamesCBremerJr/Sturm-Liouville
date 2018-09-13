!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kummer_sturm
use utils
use chebyshev
use odesolve
use kummer


implicit double precision (a-h,o-z)


double precision, private               :: pi, eulergamma,sqrt2overpi,piover2,twopi,sqrtpiover2

data pi          / 3.14159265358979323846264338327950288d0  /
data eulergamma  / 0.577215664901532860606512090082402431d0 /
data sqrt2overpi / 0.797884560802865355879892119868763737d0 /
data piover2     / 1.57079632679489661923132169163975144d0  /
data twopi       / 6.28318530717958647692528676655900577d0  /
data sqrtpiover2 / 1.25331413731550025120788264240552263d0  /

type      sturm_phase
double precision              :: dlambda

integer                       :: k
double precision, allocatable :: xscheb(:)

! data for the phase function
integer                       :: nints
double precision, allocatable :: ab(:,:)
double precision, allocatable :: alpha(:,:),alphap(:,:),alphapp(:,:)

! data for the inverse phase function
integer                       :: nintsinv
double precision, allocatable :: abinv(:,:)
double precision, allocatable :: alphainv(:,:),alphainvp(:,:)

! phase and amplitude to make evaluating functions faster 
!integer                       :: nints
!double precision, allocatable :: ab(:,:)
!double precision, allocatable :: avals(:,:),psivals(:,:)

end type  sturm_phase


type       sturm_qdata
double precision              :: dlambda
end type   sturm_qdata

contains


subroutine kummer_sturm_phase(eps,dlambda,chebdata,phase)
implicit double precision (a-h,p-z)

double precision, intent(in)         :: dlambda
type(chebexps_data), intent(in)      :: chebdata
type(sturm_phase), intent(out)       :: phase

!
!
!

type(sturm_qdata), pointer           :: qdata
type(c_ptr)                          :: userptr


allocate(qdata)
qdata%dlambda = dlambda
userptr   = c_loc(qdata)

a        = 0.0d0
b        = 1.0d0

ifleft   = 1
k        = chebdata%k

allocate(phase%xscheb(k))

phase%k       = k
phase%dlambda = dlambda
phase%xscheb  = chebdata%xs

!
!  Construct the phase functions for P_dnu^(da,db) 
!


call kummer_adap(eps,a,b,sturmq,k,chebdata%xs,chebdata%aintl,chebdata%aintr,chebdata%u, &
   phase%nints,phase%ab,phase%alphap,phase%alphapp,userptr)


call kummer_phase(ifleft,k,chebdata%xs,chebdata%aintl,chebdata%aintr,chebdata%u, &
  phase%nints,phase%ab,phase%alpha,phase%alphap,phase%alphapp)

!
!  Make sure cos( alpha(a)) = 0
!

phase%alpha = phase%alpha + pi/2

end subroutine



subroutine kummer_sturm_phase_and_inverse(eps,dlambda,chebdata,phase)
implicit double precision (a-h,p-z)

double precision, intent(in)         :: dlambda
type(chebexps_data), intent(in)      :: chebdata
type(sturm_phase), intent(out)       :: phase

!
!
!

type(sturm_qdata), pointer           :: qdata
type(c_ptr)                          :: userptr


allocate(qdata)
qdata%dlambda = dlambda
userptr   = c_loc(qdata)

a        = 0.0d0
b        = 1.0d0

ifleft   = 1
k        = chebdata%k

allocate(phase%xscheb(k))

phase%k       = k
phase%dlambda = dlambda
phase%xscheb  = chebdata%xs

!
!  Construct the phase functions for P_dnu^(da,db) 
!


call kummer_adap(eps,a,b,sturmq,k,chebdata%xs,chebdata%aintl,chebdata%aintr,chebdata%u, &
   phase%nints,phase%ab,phase%alphap,phase%alphapp,userptr)


call kummer_phase(ifleft,k,chebdata%xs,chebdata%aintl,chebdata%aintr,chebdata%u, &
  phase%nints,phase%ab,phase%alpha,phase%alphap,phase%alphapp)

!
!  Make sure cos( alpha(a)) = 0
!

phase%alpha = phase%alpha + pi/2

!
!  Compute the inverse
!

call kummer_phase_inverse(phase%nints,phase%ab,k,phase%xscheb,chebdata%aintl,chebdata%u, &
    phase%alpha,phase%alphap,phase%nintsinv,phase%abinv,phase%alphainv,phase%alphainvp)


end subroutine


subroutine kummer_zeros()
implicit double precision (a-h,o-z)

!
!  Compute zeros of the phase 
!

end subroutine


subroutine sturmq(t,val,userptr)
implicit double precision (a-h,o-z)
double precision, intent(in)   :: t
double precision, intent(out)  :: val
type(c_ptr)                    :: userptr
type(sturm_qdata), pointer     :: qdata

call c_f_pointer(userptr,qdata)
dlambda = qdata%dlambda

val = dlambda**2 / (1+t**2) * cos(t)**4 + (3+sin(3*t**2)) / (.0001d0+t**2) 

end subroutine





subroutine kummer_strum_expansion(chebdata,nints,ab,phases)
implicit double precision (a-h,o-z)

type(sturm_phase), allocatable, intent(out)     :: phases(:,:)
integer, intent(out)                            :: nints
double precision, allocatable , intent(out)     :: ab(:,:)
type(chebexps_data), intent(in)                 :: chebdata

double precision, allocatable :: dlambdas(:)

!
!
!
!
!
type(sturm_phase), allocatable                  :: phases0(:,:),phases00(:)
integer                                         :: nints0
double precision, allocatable                   :: ab0(:,:),ab00(:)
integer, allocatable                            :: istack(:),ipivs(:)

k       = chebdata%k
eps     = 1.0d-13
epserr  = 1.0d-13
maxints = 1000

allocate(phases0(k,maxints),ab0(2,maxints),istack(maxints),phases00(k))


N      = 2**(20)
dN     = N
dd     = 6.0d0
nints0 = ceiling(log(N+0.0d0)/log(dd))

if (dd**(nints0) .lt. dn) nints0 = nints0+1

do i=1,nints0
ab0(1,i) = dd**(i-1)
ab0(2,i) = min(dN,dd**(i))
end do

nstack = nints0
do i=1,nints0
istack(i) = i
end do

call prin2("in sturm_expansion, initial ab = ",ab0(:,1:nints0))

write ( *,*) ""
write (13,*) ""

write ( *,*) "---[ Building phase ]-----------------------------------------------------------"
write (13,*) "---[ Building phase ]-----------------------------------------------------------"

write ( *,*) ""
write (13,*) ""

do while (nstack .gt. 0)

int    = istack(nstack)
nstack = nstack-1

!
!  Construct the phase functions for this interval
!
a        = ab0(1,int)
b        = ab0(2,int)
dlambdas = (b-a)/2 * chebdata%xs + (b+a)/2

do i=1,k
call kummer_sturm_phase(eps,dlambdas(i),chebdata,phases00(i))
end do

!
!  Test the error
!


derr = 0
call  kummer_sturm_error0(chebdata,k,a,b,phases00,derr)


ifsplit = 0
if (derr .gt. epserr) ifsplit = 1

write (*,*)    "   ",int,a,b,derr
write (13,*)   "   ",int,a,b,derr

if (ifsplit .eq. 1) then
c = (a+b)/2

ab0(1,int) = a
ab0(2,int) = c

nints0 = nints0 + 1
ab0(1,nints0) = c
ab0(2,nints0) = b

nstack         = nstack + 1
istack(nstack) = int

nstack = nstack + 1
istack(nstack) = nints0

else

do i=1,k
phases0(i,int) = phases00(i)
end do

endif
end do



nints = nints0

allocate(ab00(nints))
allocate(ab(2,nints))
allocate(ipivs(nints))

do i= 1,nints
ipivs(i)   = i
ab00(i)    = ab0(1,i)
end do

call insort2(nints,ab00,ipivs)

do i=1,nints
int0    = ipivs(i)
ab(1,i) = ab0(1,int0)
ab(2,i) = ab0(2,int0)
end do


write ( *,*) "--------------------------------------------------------------------------------"
write (13,*) "--------------------------------------------------------------------------------"


write ( *,*) ""
write (13,*) ""



!
!  Copy out the phase functions
!


allocate(phases(k,nints))


do int=1,nints


a        = ab(1,int)
b        = ab(2,int)
int0     = ipivs(int)



do i=1,k
 ! call kummer_sturm_phase(eps,dlambdas(i),chebdata,phases(i,int))
phases(i,int) = phases0(i,int0)
end do
end do

! call prin2("second pass time = ",t2-t1)
call prin2("in kummer_sturm_expansion, final ab = ",ab)

end subroutine



subroutine kummer_sturm_phase_eval(phase,nts,ts,avals,apvals)
implicit double precision (a-h,o-z)

type(sturm_phase), intent(in)               :: phase
integer, intent(in)                         :: nts
double precision, intent(in)                :: ts(nts)
double precision, intent(out)               :: avals(nts),apvals(nts)

double precision, allocatable :: ts0(:)

int0 = 1

do i=1,nts
t = ts(i)
do int=int0,phase%nints
b = phase%ab(2,int)
if (t .le. b) exit
end do
a = phase%ab(1,int)

call chebeval_two(a,b,phase%k,phase%xscheb,phase%alpha(1,int),phase%alphap(1,int),&
  t,avals(i),apvals(i))

end do

end subroutine




subroutine kummer_sturm_eigenvalues(chebdata,nints,ab,phases,ndlambda,dlambdas)
implicit double precision (a-h,o-z)

type(sturm_phase), allocatable, intent(in)      :: phases(:,:)
integer, intent(in)                             :: nints
double precision, allocatable , intent(in)      :: ab(:,:)
type(chebexps_data), intent(in)                 :: chebdata
double precision, allocatable, intent(out)      :: dlambdas(:)

!
!
!

double precision, allocatable :: vals(:,:),ders(:,:),adiff(:,:)
double precision, allocatable :: abinv(:,:),valsinv(:,:),dersinv(:,:)



!
!  Calculate alpha(b,dlambda) and its derivative with respect to
!  lambda
!

k  = chebdata%k

allocate(vals(k,nints),ders(k,nints))
vals = 0
ders = 0

call chebdiff(k,chebdata%xs,adiff)

do int=1,nints
a = ab(1,int)
b = ab(2,int)
do i=1,k
vals(i,int) = phases(i,int)%alpha(k,phases(i,int)%nints) 
end do
ders(:,int) = matmul(adiff*2/(b-a),vals(:,int))
end do

!
!  Calculate the inverse function (hijack the kummer phase function
!  inverse routine to do it)
!

call kummer_phase_inverse(nints,ab,k,chebdata%xs,chebdata%aintl,chebdata%u, &
   vals,ders,nintsinv,abinv,valsinv,dersinv)


val0 = vals(1,1)
val1 = vals(k,nints)

iroot1 = ceiling((val0 - pi/2)/pi)
iroot2 = ceiling((val1 - pi/2)/pi)
ndlambda = iroot2-iroot1+1
allocate(dlambdas(ndlambda))


do i=1,ndlambda
xx = pi/2 + (i-1+iroot1)*pi
call chebpw_eval(nintsinv,abinv,k,chebdata%xs,valsinv,xx,root)
dlambdas(i) = root
end do

end subroutine






subroutine  kummer_sturm_error(chebdata,nints,ab,phases,int,derr)
implicit double precision (a-h,o-z)

type(sturm_phase), allocatable, intent(in)      :: phases(:,:)
integer, intent(in)                             :: nints
double precision, allocatable , intent(in)      :: ab(:,:)
type(chebexps_data), intent(in)                 :: chebdata

!
!
!

double precision, allocatable  :: ts(:),dnus(:),avals(:,:),apvals(:,:),avals0(:),apvals0(:),derrs(:)
type(sturm_phase)              :: phase

eps0  = epsilon(0.0d0)
eps   = eps0*100

nts   = 100
ndnus = 100

allocate(ts(nts),dnus(ndnus),avals(nts,ndnus),apvals(nts,ndnus))
allocate(avals0(nts),apvals0(nts),derrs(nts))



do i=1,nts
call random_number(dd)
ts(i) = dd
end do

a = ab(1,int)
b = ab(2,int)

do i=1,ndnus
!call random_number(dd)
dd = (i-1.0d0)/(ndnus-1.0d0)
dnus(i) = a + dd * (b-a)
end do

call quicksort(nts,ts)
call quicksort(ndnus,dnus)


call kummer_sturm_interpolate(chebdata,nints,ab,phases,nts,ts,ndnus,dnus,avals,apvals)



errmax1 = 0
errmax2 = 0

do i = 1, ndnus

call kummer_sturm_phase(eps,dnus(i),chebdata,phase)
call kummer_sturm_phase_eval(phase,nts,ts,avals0,apvals0)

!derrs  = abs( (apvals0-apvals(:,i)) / apvals0 )
derrs  = abs( (avals0-avals(:,i)) / avals0 )

! call prin2("dnu = ",dnus(i))
! call prin2("derrs = ",derrs)

derr2   = maxval(derrs)
errmax2 = max(errmax2,derr2)

end do

derr = errmax2

end subroutine






subroutine chebeval_two(a,b,n,xs,vals1,vals2,x,val1,val2)
implicit double precision (a-h,o-z)

integer :: n
double precision ::  xs(n),vals1(n),vals2(n),x,val1,val2


eps0 = epsilon(0.0d0)

xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0
sum3=0

dd1 = 1.0d0

do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0

diff = xx-xs(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0) then
val1 = vals1(i)
val2 = vals2(i)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals1(i)
sum2 = sum2+dd
sum3 = sum3+dd*vals2(i)
dd   = - dd
end do

val1 = sum1/sum2
val2 = sum3/sum2

end subroutine


subroutine kummer_sturm_interpolate(chebdata,nints,ab,phases,nts,ts,ndnus,dnus,avals,apvals)
implicit double precision (a-h,o-z)

type(chebexps_data), intent(in)                 :: chebdata
type(sturm_phase), allocatable, intent(in)      :: phases(:,:)
integer, intent(in)                             :: nints
double precision, allocatable , intent(in)      :: ab(:,:)

integer, intent(in)                             :: nts,ndnus
double precision, intent(in)                    :: ts(nts), dnus(ndnus)
double precision, intent(out)                   :: avals(nts,ndnus),apvals(nts,ndnus)
!
!

double precision, allocatable :: avals0(:,:),apvals0(:,:)

k = chebdata%k

if (nts .eq. 0 .OR. ndnus .eq. 0) return

allocate(avals0(nts,k),apvals0(nts,k))


!
!  Traverse the list of dnu's, find all dnus in a given interval
!  and interpolate.
!

int0 = 1
i1   = 1

0100 continue


do int = int0,nints
b = ab(2,int)
if (dnus(i1) .lt. b) goto 1000
end do

1000 continue
a = ab(1,int)



do i2=i1,ndnus-1
if (dnus(i2+1) .gt. b) exit
end do


!
!  Do the interpolation
!

do j=1,k
dnu = (b-a)/2 * chebdata%xs(j) + (b+a)/2
call kummer_sturm_phase_eval(phases(j,int),nts,ts,avals0(:,j),apvals0(:,j))
end do

do i=i1,i2
do l=1,nts
call chebeval_two(a,b,k,chebdata%xs,avals0(l,:),apvals0(l,:),dnus(i),avals(l,i),apvals(l,i))
end do
end do


if (i2 .lt. ndnus) then
i1 = i2+1
goto 0100
endif



end subroutine


subroutine  kummer_sturm_error0(chebdata,k,a,b,phases,derr)
implicit double precision (a-h,o-z)

type(sturm_phase), intent(in)                   :: phases(k)
type(chebexps_data), intent(in)                 :: chebdata

!
!
!

double precision, allocatable  :: ts(:),dnus(:),avals(:,:),apvals(:,:),avals0(:),apvals0(:),derrs(:)
type(sturm_phase)              :: phase

eps0  = epsilon(0.0d0)
eps   = eps0*100

nts   = 40
ndnus = 6


allocate(ts(nts),dnus(ndnus),avals(nts,ndnus),apvals(nts,ndnus))
allocate(avals0(nts),apvals0(nts),derrs(nts))


do i=1,nts
call random_number(dd)
ts(i) = dd
end do

call elapsed(t1)
n = t1
call random_seed(n)

do i=1,ndnus
call random_number(dd)
!dd = (i-1.0d0)/(ndnus-1.0d0)
dnus(i) = a + dd * (b-a)
end do

call quicksort(nts,ts)
call quicksort(ndnus,dnus)


call kummer_sturm_interpolate0(chebdata,k,a,b,phases,nts,ts,ndnus,dnus,avals,apvals)


errmax1 = 0
errmax2 = 0

do i = 1, ndnus

call kummer_sturm_phase(eps,dnus(i),chebdata,phase)
call kummer_sturm_phase_eval(phase,nts,ts,avals0,apvals0)

derrs   = abs( (avals0-avals(:,i)) / avals0 )
derr2   = maxval(derrs)
errmax2 = max(errmax2,derr2)

end do

derr = errmax2

end subroutine




subroutine kummer_sturm_interpolate0(chebdata,k,a,b,phases,nts,ts,ndnus,dnus,avals,apvals)
implicit double precision (a-h,o-z)

type(chebexps_data), intent(in)                 :: chebdata
type(sturm_phase), intent(in)                   :: phases(k)

integer, intent(in)                             :: nts,ndnus
double precision, intent(in)                    :: ts(nts), dnus(ndnus)
double precision, intent(out)                   :: avals(nts,ndnus),apvals(nts,ndnus)
!
!

double precision, allocatable :: avals0(:,:),apvals0(:,:)

k = chebdata%k

if (nts .eq. 0 .OR. ndnus .eq. 0) return

allocate(avals0(nts,k),apvals0(nts,k))


!
!  Do the interpolation
!

do j=1,k
dnu = (b-a)/2 * chebdata%xs(j) + (b+a)/2
call kummer_sturm_phase_eval(phases(j),nts,ts,avals0(:,j),apvals0(:,j))
end do

do i=1,ndnus
do l=1,nts
call chebeval_two(a,b,k,chebdata%xs,avals0(l,:),apvals0(l,:),dnus(i),avals(l,i),apvals(l,i))
end do
end do



end subroutine


end module
 
