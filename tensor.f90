!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!  This module contains code for representing a function of two variables given on
!  a rectangle [a,b] x [c,d] as a tensor product of Chebyshev polynomials;  that is,
!  in the form
!                    
!                n-1   n-1                  
!     f(x,y) =   sum   sum     a_ij T_{1,i}(x) T_{2,j}(y)                                     (1)
!                j=0   i=0
!
!
!  where T_{1,i} denotes the Chebyshev polynomial of degree i on the interval [a,b]
!  and T_{2,i} denote the Chebyshev polynomial of degree i on the interval [c,d].
!
!  It also contains code for storing a compressed version of (1) which reduces
!  the number of terms in the expansion by neglecting those with coefficients which
!  are small in magnitude.
!
!  More specifically, the compressed form of (1) is a sum of the form
!
!                m-1   n_j                  
!     f(x,y) =   sum   sum     a_ij T_{1,i}(x) T_{2,j}(y)                                     (2)
!                j=0   i=0
!
!
!  with the integers m < n and -1 <= n_j < m_1 chosen such that
!
!      |a_ij| < \epsilon     for j > m
!
!      |a_ij| < \epsilon     for i > n_j.
!
!  The following subroutines should be regarded as public:
!
!    tensor_umatrix - return the matrix which takes the values of a funtion f at
!      nodes of a tensor product of Clenshaw-Curtis quadrature formulas to the 
!      coefficients in the expansion (1)
!
!    tensor_coefs - a utility routine for applying the matrix returned by
!      tensor_umatrix to the vector of values; its raison d'etre is that
!      it provides a mechanism for reshaping the array of values implicitly
!
!    tensor_eval - calculate the value of an expansion of the form (1)  at 
!      a specified point
!
!    tensor_evalder - calculate the value of an expansion of the form (1) and
!      its derivative at a specified point
!
!    tensor_compressed_coefs - return an array which stores the compresses coefficients
!      in an expansion of the form (2) 
!
!    tensor_compressed_eval - evaluate a compressed expansion of the form (2) at
!      a specified point
!
!    tensor_compressed_eval2 - a special version of the routine which evaluates
!      two expansions of the form (2) simultaneously; this is somewhat faster
!      than evaluating them separately because it obviates the need to 
!      evaluate the Chebyshev polynomials twice
!
!    tensor_compressed_evalder - evaluate a compressed expansion of the form (2)
!      and its derivative at a specified point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tensor

use chebyshev

contains

subroutine tensor_umatrix(n,u)
implicit double precision (a-h,o-z)

integer, intent(in)                         :: n
double precision, allocatable, intent(out)  :: u(:,:)

!
!  Return the matrix which takes the vector of values
!
!      f( xs(1)    , xs(1) )
!      f( xs(2)    , xs(1) )
!               .
!               .
!      f( xs(n)    , xs(1) )
!      f( xs(1)    , xs(2) )
!      f( xs(2)    , xs(2) )                                                      (2)
!               .
!               .
!
!      f( xs(1)    , xs(n) )
!      f( xs(2)    , xs(n) )
!               .
!               .
!      f( xs(n )   , xs(n) ),
!
!  where xs(1), ..., xs(n) are the nodes of the n-point Clenshaw-Curtis quadrature
!  formula to the vector of coefficients
!
!      a_00
!      a_01
!      a_10
!      a_02                                                                       (3)
!      a_11
!      a_20
!       .
!       .
!       .
!    
!
!  in the expansion (1).  Note that there are (n^2 + n)/2 such coefficients, so that
!  the matrix u is dimensioned ( (n**2 + n)/2 , n*n ).  Also note tht the coefficients
!  in (3) are sorted by the degree of the corresponding term in the expansion (1).
!
!  Input parameters:
!    n - the number of points in the Chebyshev grid on [-1,1] 
!    
!  Output parameters:
!    u - the ((n+1)*n/2, n*n) matrix which takes the vector of values (2) 
!     to the vector of coefficients (3)
!
!

double precision, allocatable :: polsx(:),polsy(:),xs(:)


allocate(u(n*n,n*n))

!
!  Fetch the Chebyshev nodes
!
allocate(polsx(n+1), polsy(n+1),xs(n))

call chebpts(n,xs)


!
!  Use the identity
!
!
!              n  ''                        { \delta_ij * (n)   if i = 0 or i = n
!             \sum   T_i(x_k) T_j(x_k)  =   {
!              k=0                          { \delta_ij * (n/2) if 0 < i < n
!
!  to form the matrix u0 which takes the values of an expansion of the form
!
!               n-1  n-1
!      f(x,y) = sum  sum   a_{ij} T_i(x) T_j(x) 
!               i=0  j=0
!
!  to the coefficients in (1).
!

nx = n
ny = n

do i1=1,nx
x    = xs(i1)
call chebs(x,nx,polsx)
if (i1 .eq. 1 .OR. i1 .eq. nx) polsx = polsx/2.0d0

do i2=1,ny
y    = xs(i2)
call chebs(y,ny,polsy)
if (i2 .eq. 1 .OR. i2 .eq. ny) polsy = polsy/2.0d0

do j1=1,nx
do j2=1,ny
val    = polsy(j2)*polsx(j1)

dscaley = 1.0d0/(ny-1)
if (j2 .gt. 1 .AND. j2 .lt. ny) dscaley=dscaley*2

dscalex = 1.0d0/(nx-1)
if (j1 .gt. 1 .AND. j1 .lt. nx) dscalex=dscalex*2

u(j1 + nx*(j2-1),i1+(i2-1)*nx) = val * dscaley * dscalex

end do
end do

end do
end do


end subroutine


subroutine tensor_coefs(n,u,vals,coefs)
implicit double precision (a-h,o-z)

integer          ::  n
double precision ::  vals(n*n),coefs(n*n),u(n*n,n*n)

!
!  This is a simple utility code which applies the matrix u returned by
!  tensor_umatrix to a vector.  Its raison d'etre is that it allows
!  the user to shape the input vector vals as an (n,n) matix.  That is,
!  this routine is used to reshape the input vector vals.
!
!  Input parameters:
!    n - the number of points in the Chebyshev grid
!    u - the ((n+1)*n/2,n*n) matrix returned by tensor_umatrix  
!    vals - a (n,n) matrix containing the values of the expansion
!
!  Output parameters:
!    coefs- the coefficients in the expansion
!

coefs = matmul(u,vals)

end subroutine




subroutine tensor_eval(n,a,b,c,d,coefs,x,y,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n
double precision, intent(in)  :: a,b,c,d,x,y,coefs(n,n)
double precision, intent(out) :: val

!
!  Evaluate an expansion of the form (1) at a specified point.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid on [-1,1]
!    (a,b) - the interval over which x varies
!    (c,d) - the interval over which y varies
!    coefs - the array of coefficients in the expansion
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!
!   val - the value of the expansion at the desired point
!   derx - the value of the derivative of the expansion w.r.t. x
!   dery - the value of the derivative of the expansion w.r.t. y
!

double precision :: polsx(n),polsy(n)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebs(xx,n-1,polsx)
call chebs(yy,n-1,polsy)

nn = 0
val  = 0
do j=1,n
do i=1,n
if ( (i-1)**2 + (j-1)**2 .le. n**2) then
nn = nn + 1
val  = val  + coefs(i,j) * polsx(i)*polsy(j)
endif

end do
end do

end subroutine


subroutine tensor_evalder(n,a,b,c,d,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)

integer, intent(in)           :: n
double precision, intent(in)  :: a,b,c,d,x,y,coefs(n,n)
double precision, intent(out) :: val

!
!  Evaluate an expansion of the form (1) at a specified point as well
!  as its derivative.
!
!  Input parameters:
!
!    n - the number of points in the Chebyshev grid on [-1,1]
!    (a,b) - the interval over which x varies
!    (c,d) - the interval over which y varies
!    coefs - the array of coefficients in the expansion
!    (x,y) - the point at which to evaluate the expansion
!
!  Output parameters:
!
!   val - the value of the expansion at the desired point
!   derx - the value of the derivative of the expansion w.r.t. x
!   dery - the value of the derivative of the expansion w.r.t. y
!

double precision :: polsx(n+10),polsy(n+10),dersx(n+10),dersy(n+10)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebders(xx,n-1,polsx,dersx)
call chebders(yy,n-1,polsy,dersy)

val   = 0
derx  = 0 
dery  = 0
nn    = 0

do j=1,n
do i=1,n
dd    = coefs(i,j)


if ( (i-1)**2+(j-1)**2 .le. (n-1)**2) then

nn    = nn + 1
val   = val   + dd * polsx(i)*polsy(j)
derx  = derx  + dd * dersx(i)*polsy(j)
dery  = dery  + dd * polsx(i)*dersy(j)

endif
end do
end do

derx = derx * 2/(b-a)
dery = dery * 2/(d-c)

print *,nn,n**2
end subroutine


subroutine tensor_compressed_coefs(eps,n,coefs,ncompcoefs,compcoefs)
implicit double precision (a-h,o-z)

double precision, intent(in)               :: eps,coefs(n,n)
double precision, allocatable,intent(out)  :: compcoefs(:)

!
!
!  Input parameters:
!
!  Output parameters:
!
!

integer, allocatable :: iptrs(:)

allocate(iptrs(n))

!
!  For each j, find the smallest n_j such that |a_ij| < \epsilon for all
!  i > n_j and record the index in the iptrs array.
!

ncompcoefs = 0
nx         = 0
ny         = 0

dmax = maxval(abs(coefs))

do j=0,n-1
nj = -1

do i=0,n-1
if (abs(coefs(i+1,j+1)/dmax)  .gt. eps)  nj = i
end do

iptrs(j+1)   = nj
ncompcoefs   = ncompcoefs+nj+1
if (nj .ge. 0) ny = max(ny,j)
nx = max(nj,nx)

end do


ncompcoefs = ncompcoefs + 2 + ny+1
allocate(compcoefs(ncompcoefs))


compcoefs(1) = nx
compcoefs(2) = ny


!
!  Copy out the nonzero coefficients.
!

iptr = 3

do j=0,ny
nj              = iptrs(j+1)
compcoefs(iptr) = nj
iptr            = iptr+1

do i=0,nj
compcoefs(iptr) = coefs(i+1,j+1)
iptr = iptr+1
end do
end do


end subroutine





subroutine tensor_compressed_eval(compcoefs,a,b,c,d,x,y,val)
implicit double precision (a-h,o-z)

double precision :: compcoefs(:)

double precision :: polsx(100),polsy(100)

nx     = compcoefs(1)
ny     = compcoefs(2)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebs(xx,nx,polsx)
call chebs(yy,ny,polsy)


val  = 0
iptr = 3

do j=0,ny
nj    =  compcoefs(iptr)
iptr  = iptr+1
do i=0,nj
val    = val + compcoefs(iptr)*polsx(i+1)*polsy(j+1)
iptr   = iptr+1
end do
end do

end subroutine


subroutine tensor_compressed_eval2(compcoefs,iptr0,a,b,c,d,x,y,val)
implicit double precision (a-h,o-z)

double precision :: compcoefs(:)

double precision :: polsx(100),polsy(100)

nx     = compcoefs(iptr0)
ny     = compcoefs(iptr0+1)


xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebs(xx,nx,polsx)
call chebs(yy,ny,polsy)


val  = 0
iptr = iptr0+2

do j=0,ny
nj    =  compcoefs(iptr)
iptr  = iptr+1
do i=0,nj
val    = val + compcoefs(iptr)*polsx(i+1)*polsy(j+1)
iptr   = iptr+1
end do
end do

end subroutine



subroutine tensor_compressed_eval22(compcoefs1,compcoefs2,iptr1,iptr2,a,b,c,d,x,y,val1,val2)
implicit double precision (a-h,o-z)

double precision :: compcoefs1(:),compcoefs2(:)
double precision :: polsx(100),polsy(100)

nx1     = compcoefs1(iptr1)
ny1     = compcoefs1(iptr1+1)

nx2     = compcoefs2(iptr2)
ny2     = compcoefs2(iptr2+1)

nx      = max(nx1,nx2)
ny      = max(ny1,ny2)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebs(xx,nx,polsx)
call chebs(yy,ny,polsy)

val1   = 0
iptr   = iptr1+2

do j=0,ny1
nj    =  compcoefs1(iptr)
iptr  = iptr+1
do i=0,nj
val1    = val1 + compcoefs1(iptr)*polsx(i+1)*polsy(j+1)
iptr   = iptr+1
end do
end do

val2   = 0
iptr   = iptr2+2
do j=0,ny2
nj    = compcoefs2(iptr)
iptr  = iptr+1
do i=0,nj
val2   = val2 + compcoefs2(iptr)*polsx(i+1)*polsy(j+1)
iptr   = iptr+1
end do
end do

end subroutine



subroutine tensor_compressed_evalder(compcoefs,a,b,c,d,x,y,val,derx,dery)
implicit double precision (a-h,o-z)

double precision :: compcoefs(:)
double precision :: polsx(100),polsy(100),dersx(100),dersy(100)

nx     = compcoefs(1)
ny     = compcoefs(2)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebders(xx,nx,polsx,dersx)
call chebders(yy,ny,polsy,dersy)


val  = 0
derx = 0
dery = 0
iptr = 3

do j=0,ny
nj    =  compcoefs(iptr)
iptr  = iptr+1
do i=0,nj
val    = val  + compcoefs(iptr)*polsx(i+1)*polsy(j+1)
derx   = derx + compcoefs(iptr)*dersx(i+1)*polsy(j+1)
dery   = dery + compcoefs(iptr)*polsx(i+1)*dersy(j+1)
iptr   = iptr+1
end do
end do

derx = derx * 2/(b-a)
dery = dery * 2/(d-c)

end subroutine


subroutine tensor_compressed_evalder2(compcoefs,iptr0,a,b,c,d,x,y,val,derx,dery)
implicit double precision (a-h,o-z)

double precision :: compcoefs(:)
double precision :: polsx(100),polsy(100),dersx(100),dersy(100)

nx     = compcoefs(iptr0)
ny     = compcoefs(iptr0+1)

xx = (x-(b+a)/2)*2/(b-a)
yy = (y-(d+c)/2)*2/(d-c)

call chebders(xx,nx,polsx,dersx)
call chebders(yy,ny,polsy,dersy)


val  = 0
derx = 0
dery = 0
iptr = iptr0+2

do j=0,ny
nj    =  compcoefs(iptr)
iptr  = iptr+1
do i=0,nj
val    = val  + compcoefs(iptr)*polsx(i+1)*polsy(j+1)
derx   = derx + compcoefs(iptr)*dersx(i+1)*polsy(j+1)
dery   = dery + compcoefs(iptr)*polsx(i+1)*dersy(j+1)
iptr   = iptr+1
end do
end do

derx = derx * 2/(b-a)
dery = dery * 2/(d-c)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine tensor_compressed2_coefs(eps,n,coefs,ncompcoefs,compcoefs,compis,compjs)
implicit double precision (a-h,o-z)

double precision, intent(in)               :: eps,coefs(n,n)
double precision, allocatable,intent(out)  :: compcoefs(:)
integer*1, allocatable, intent(out)        :: compis(:),compjs(:)


!
!  Count the number of coeffcients whose magniutdes are above the threshold.
!

dmax = maxval(abs(coefs))

ncompcoefs = 0

do i=1,n
do j=1,n

coef = coefs(i,j)
dd   = abs(coef)/dmax
if (i+j .le. n+1 .AND. dd .gt. eps) then
ncompcoefs = ncompcoefs + 1
endif

end do 
end do

allocate(compis(ncompcoefs+3),compjs(ncompcoefs),compcoefs(ncompcoefs))


!
!   Copy out the coefficients
! 

ncompcoefs = 0
nx         = 0 
ny         = 0

do i=1,n
do j=1,n

coef = coefs(i,j)
dd   = abs(coef)/dmax

if (dd .gt. eps) then
ncompcoefs         = ncompcoefs + 1
compis(ncompcoefs+3)  = i-1
compjs(ncompcoefs)    = j-1
nx                    = max(i-1,nx)
ny                    = max(j-1,ny)
compcoefs(ncompcoefs) = coef
endif

end do 
end do


compis(1) = ncompcoefs
compis(2) = nx
compis(3) = ny

end subroutine


subroutine tensor_compressed2_eval(eps,n,ncompcoefs,compcoefs,compis,compjs,a,b,c,d,x,y,val)
implicit double precision (a-h,o-z)

double precision, intent(in)               :: eps
double precision, intent(in)               :: compcoefs(ncompcoefs)
integer*1, intent(in)                      :: compis(ncompcoefs+3),compjs(ncompcoefs)

double precision :: polsx(100),polsy(100)

!
!  Count the number of coeffcients whose magniutdes are above the threshold.
!

nx   = compis(2)
ny   = compis(3)

xx   = (x-(b+a)/2)*2/(b-a)
yy   = (y-(d+c)/2)*2/(d-c)

call chebs(xx,nx,polsx)
call chebs(yy,ny,polsy)

idx        = 1
val        = 0 
do idx=1,ncompcoefs
i          = compis(idx+3)
j          = compjs(idx)
coef       = compcoefs(idx)
val        = val + coef * polsx(i+1)*polsy(j+1)
end do

end subroutine



end module
