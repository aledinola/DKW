#include "fintrf.h"
!----------------------------------------------------------------!
! Gateway routine
! MAC
! Compile with mex -O -v -R2018a sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex
! mex -v FOPTIMFLAGS="-O3 -qopenmp "  sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex
! WINDOWS 
! mex -v OPTIMFLAGS="/O3 /Qopenmp /Qprec-div- /QxHost /DNDEBUG" -R2018a sub_vfi_onestep_mex.f90 -output sub_vfi_onestep_mex
! USAGE in Matlab
! [val_c_new,pol_kp] = sub_vfi_onestep_mex(val_c,val0_u,kp_bar,B_hat,profit_mat,k_grid,...
!	b_grid,pi_x,theta,q,delta,psi,int32(do_howard),int32(n_howard))
! where val_c must be a 3-dim (nk,nb,nx) array and val_c_new as well
! @Alessandro Di Nola, July 29, 2022
!----------------------------------------------------------------!
subroutine mexFunction(nlhs,plhs,nrhs,prhs)
    implicit none 
    ! Arguments in mexFunction
    mwPointer :: plhs(*), prhs(*)
    integer   :: nlhs, nrhs
    
    ! Function declarations
    mwPointer, external :: mxGetPr, mxGetNumberOfDimensions, mxGetDimensions
    mwPointer, external :: mxCreateNumericArray
    integer, external   :: mxIsNumeric
    integer(4), external   :: mxClassIDFromClassName
    
    !  Array size information:
    mwSize :: ndim
    mwSize :: dims(3)
    mwSize :: nk, nb, nx
    
    ! Array type information
    integer(4) :: classid, classid2, ComplexFlag
    
    ! Pointers to input/output arrays
    ! - Inputs
    mwPointer :: val_c,val0_u,kp_bar,B_hat,profit_mat,k_grid,b_grid,pi_x,&
    	theta,q,delta,psi,do_howard,n_howard

    ! - Outputs
    mwPointer :: val_c_new,pol_kp

    ! Others
    character (len=3) :: c
    character (len=50) :: txt
    
    ! USEFUL INFO FROM MATLAB WEBSITE
    ! #include "fintrf.h"
    ! mwPointer mxCreateNumericArray(ndim, dims, classid, ComplexFlag)
    ! mwSize ndim
    ! mwSize dims(ndim)
    ! integer*4 classid, ComplexFlag
!----------------------------------------------------------------!
    
    ! Check for proper number of arguments. 
    if(nrhs /= 14) then
        call mexErrMsgIdAndTxt('MATLAB:sub_vfi_onestep_mex:nInput','14 inputs required.')
    elseif(nlhs > 2) then
        call mexErrMsgIdAndTxt('MATLAB:sub_vfi_onestep_mex:nOutput','Too many output arguments.')
    endif
    
    ! Check that the input is a number.
    if(mxIsNumeric(prhs(1)) == 0) then
        call mexErrMsgIdAndTxt('MATLAB:sub_vfi_onestep_mex:NonNumeric','First input must be a numeric variable.')
    endif
    
    !Get no. of dimensions of input array, dim:(n1,n2,n3)
    ndim = mxGetNumberOfDimensions(prhs(1))
    
    ! Get size of input argument (also size of output argument)
    call mxCopyPtrToPtrArray(mxGetDimensions(prhs(1)),dims,mxGetNumberOfDimensions(prhs(1)))
    
    nk = dims(1)
    nb = dims(2)
    nx = dims(3)
    
    classid     = mxClassIDFromClassName('double')
    classid2    = mxClassIDFromClassName('int32')
    ComplexFlag = 0
    
    ! DEBUG Display on screen 
    ! Specify alert to be printed
    !write(txt,'(A)') 'n3 = ' 
    ! Convert nrhs to char
    !write(c,'(i3)') n3
    !! Call the API
    !call mexPrintf(txt)
    !call mexPrintf(c)
    
    ! Retrieve input arguments pointers
    ! val_c,val0_u,kp_bar,B_hat,profit_mat,k_grid,b_grid,pi_x,theta,q,delta
    val_c      = mxGetPr(prhs(1))
    val0_u     = mxGetPr(prhs(2))
    kp_bar     = mxGetPr(prhs(3))
    B_hat      = mxGetPr(prhs(4))
    profit_mat = mxGetPr(prhs(5))
    k_grid     = mxGetPr(prhs(6))
    b_grid     = mxGetPr(prhs(7))
    pi_x       = mxGetPr(prhs(8))
    theta      = mxGetPr(prhs(9))
    q          = mxGetPr(prhs(10))
    delta      = mxGetPr(prhs(11))
    psi        = mxGetPr(prhs(12))
    do_howard  = mxGetPr(prhs(13))
    n_howard   = mxGetPr(prhs(14)) 

    ! Create arrays for output arguments
    plhs(1) = mxCreateNumericArray(ndim,dims,classid,ComplexFlag)
    plhs(2) = mxCreateNumericArray(ndim,dims,classid2,ComplexFlag)
    
    ! Get pointer for output argument
    val_c_new = mxGetPr(plhs(1))
    ! TODO: pol_kp is real or integer??
    pol_kp    = mxGetPr(plhs(2))
    
    ! Call the computational subroutine (LHS: output, RHS: inputs)
    ! [val_c_new,pol_kp] = sub_vfi_onestep_mex(val_c,val0_u,kp_bar,B_hat,profit_mat,k_grid,...
    !	b_grid,pi_x,theta,q,delta,int32(do_howard),int32(n_howard))
    call sub_vfi_onestep(%val(val_c_new),%val(pol_kp), %val(val_c),%val(val0_u),%val(kp_bar),&
        %val(B_hat),%val(profit_mat),%val(k_grid),%val(b_grid),%val(pi_x),%val(theta),&
        %val(q),%val(delta),%val(psi),%val(do_howard),%val(n_howard),nk,nb,nx)
    
    !call vfi_one_step_mex_v1(%val(val_new), %val(val),%val(nb),%val(nx))
    
end subroutine mexFunction
!============================================================================!

module mymod
implicit none
    contains
    
PURE FUNCTION locate(xx,x)
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: xx
    REAL(8), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,il,im,iu
    n=size(xx)
    il=0
    iu=n+1
    do
	    if (iu-il <= 1) exit
	    im=(iu+il)/2
	    if (x >= xx(im)) then
		    il=im
	    else
		    iu=im
	    end if
    end do
    if (x == xx(1)) then
	    locate=1
    else if (x == xx(n)) then
	    locate=n-1
    else
	    locate=il
    end if
END FUNCTION locate
!===============================================================================!

PURE FUNCTION locate_equi(x_grid,xi) result(jl)
    IMPLICIT NONE
	! Declare inputs
    REAL(8), INTENT(IN) :: x_grid(:)
    REAL(8), INTENT(IN) :: xi
	! Declare function result 
    INTEGER :: jl
    INTEGER :: n
	REAL(8) :: step, xi_min 
    
	n      = size(x_grid)
    step   = x_grid(2)-x_grid(1)
	xi_min = xi-x_grid(1)
	
	jl     = min(n-1,max(1,floor(xi_min/step)+1))
	
END FUNCTION locate_equi
!===============================================================================!

pure subroutine find_loc(jstar,omega, a_grid,aopt) 
     implicit none
     ! Declare Inputs:
     real(8), intent(in) :: a_grid(:)
     real(8), intent(in) :: aopt
     ! Declare outputs:
     integer, intent(out) :: jstar
     real(8), intent(out) :: omega
     ! Declare local variables:
     integer :: na
     
     !Find j s.t. a_grid(j)<=aopt<a_grid(j+1)
     !for j=1,..,N-1
     ! omega is the weight on a_grid(j) and lies in [0,1]
 
     na = size(a_grid)
 
     jstar = max(min(locate_equi(a_grid,aopt),na-1),1)
     !Weight on a_grid(j)
     omega = (a_grid(jstar+1)-aopt)/(a_grid(jstar+1)-a_grid(jstar))
     !omega = max(min(omega,1.0d0),0.0d0)
     
end subroutine find_loc

FUNCTION linint(x,y,xi)
    ! Purpose: linear interpolation of function y on grid x at interpolation point xi
    !          To make it pure, cannot use PRINT or STOP statements
    IMPLICIT NONE
    REAL(8), DIMENSION(:), INTENT(IN) :: x,y
    REAL(8), INTENT(IN) :: xi
    REAL(8) :: linint
    REAL(8) :: a,b,d
    INTEGER :: n,i
    n=size(x)
    IF (size(y)/=n) THEN
	    PRINT *, 'linint: x and y must be of the same size'
        PAUSE
	    STOP 'program terminated by linint'
    END IF
	! If the grid x is NOT equally spaced, use "locate"
	! otherwise use "locate_equi", faster.
    !i=max(min(locate(x,xi),n-1),1)
    i = locate_equi(x,xi)
	d = x(i+1)-x(i)
    !IF (d == 0.0) STOP 'bad x input in splint'
    a      = (x(i+1)-xi)/d
    b      = (xi-x(i))/d
    linint = a*y(i)+b*y(i+1)
END FUNCTION linint

function fun_adjcost(kprime,k,theta,delta) result(adj)
    implicit none
    ! Variables
    ! INPUTS
    ! kprime: Next-period capital, scalar
    ! k: Current period capital, scalar
    ! theta: resale value, in (0,1)
    real(8), intent(in) :: kprime, k 
    real(8), intent(in) :: theta, delta 
    real(8) :: adj
    
    ! Body of fun_adjcost
    if (kprime>=(1.0d0-delta)*k) then
        adj = kprime-(1.0d0-delta)*k
    else
        adj = theta*(kprime-(1.0d0-delta)*k)
    endif
    
end function fun_adjcost
    
end module mymod
!===============================================================================!

subroutine sub_vfi_onestep(val_c_new,pol_kp,&
    val_c,val0_u,kp_bar,B_hat,profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,do_howard,n_howard,nk,nb,nx)
    ! TODO: as inputs, need to have both this and next period's b_grid and b_hat 
    ! For the steady state, b_grid and b_hat do not depend on time. But for the transition, 
    !	they change over time!
    use omp_lib 
    use mymod, only: linint, fun_adjcost
    implicit none
    ! Declare input variables
    integer(4), intent(in) :: nk,nb,nx
    integer(4), intent(in) :: do_howard,n_howard
    real(8), intent(in)    :: val_c(nk,nb,nx)
    real(8), intent(in)    :: val0_u(nk,nb,nx), kp_bar(nk,nb,nx)
    real(8), intent(in)    :: B_hat(nk,nx), profit_mat(nk,nx)
    real(8), intent(in)    :: k_grid(nk), b_grid(nk,nb), pi_x(nx,nx)
    real(8), intent(in)    :: theta, q, delta, psi
    ! Declare output variables
    real(8), intent(out)    :: val_c_new(nk,nb,nx)
    integer(4), intent(out) :: pol_kp(nk,nb,nx)
    ! Declare local variables
    integer(4) :: k_c, b_c, x_c, xp_c, kp_c, max_ind,h_c
    logical :: liq1, liq2
    real(8) :: k_val, b_val, profit_val, kp_ub, kp_val, bprime, v0_int
    real(8), allocatable :: EVx(:,:), b_grid_kp(:), rhs_vec(:)
    real(8), allocatable :: val0_c(:,:,:), val0(:,:,:)
    integer, allocatable :: is_c(:,:,:)
    
    ! The number of threads cannot be greater than the number of cores on the computer.
    ! For haomin's imac, this is 10. 
    call omp_set_num_threads(8)
    
    ! Initialize outputs
    val_c_new = 0.0d0 !(nk,nb,nx)
    pol_kp    = 1     !(nk,nb,nx)

    ! STEP 2 - Impose liquidation, eq. (24)

    ! Value of constrained firm before exit
    allocate(val0_c(nk,nb,nx))
    val0_c = 0.0d0
    do x_c = 1,nx
        do b_c = 1,nb
            do k_c = 1,nk
                k_val = k_grid(k_c)
                ! TODO: this b_grid should be b_grid_next
                b_val = b_grid(k_c,b_c)
                profit_val = profit_mat(k_c,x_c)
                liq1 = profit_val-b_val+theta*(1.0d0-delta)*k_val<0.0d0
                liq2 = val_c(k_c,b_c,x_c)<theta*(1.0d0-delta)*k_val-b_val
                if (liq1 .or. liq2) then
                    val0_c(k_c,b_c,x_c) = theta*(1.0d0-delta)*k_val-b_val
                else
                    val0_c(k_c,b_c,x_c) = val_c(k_c,b_c,x_c)
                endif
            enddo
        enddo
    enddo
    
    ! STEP 3 - Do equation (23)
    allocate(val0(nk,nb,nx),is_c(nk,nb,nx))
    val0 = 0.0d0
    is_c = 1
    do x_c = 1,nx
        do b_c = 1,nb
            do k_c = 1,nk
            	! TODO: b_grid and b_hat should both be next period's
                b_val = b_grid(k_c,b_c)
                if (b_val<=B_hat(k_c,x_c)) then
                    val0(k_c,b_c,x_c) = val0_u(k_c,b_c,x_c)
                    ! TODO: is_c should be based on this period's b_grid and this period's b_hat
                    is_c(k_c,b_c,x_c) = 0
                else
                    val0(k_c,b_c,x_c) = val0_c(k_c,b_c,x_c)
                endif
            enddo
        enddo
    enddo
    
    ! STEP 4 - Solve eq. (25)
    
    allocate(EVx(nk,nb),rhs_vec(nk))
    
    ! Maximize over k' on the grid
	!$omp parallel default(shared) private(x_c,xp_c,b_c,k_c,EVx,k_val,b_val,profit_val,&
	!$ kp_ub,rhs_vec,kp_val,b_grid_kp,bprime,v0_int,max_ind) 
    !$omp do collapse(1)
    do x_c = 1,nx
	    ! Compute expected value in eq. (25)
        EVx = 0.0d0 !(k',b')
        do xp_c = 1,nx
            EVx = EVx+val0(:,:,xp_c)*pi_x(x_c,xp_c)
        enddo
        do b_c = 1,nb
            do k_c = 1,nk                
                k_val = k_grid(k_c)
                ! TODO: this b_grid is this period's b_grid
                b_val = b_grid(k_c,b_c)
                profit_val = profit_mat(k_c,x_c)
                kp_ub = kp_bar(k_c,b_c,x_c)

                if (is_c(k_c,b_c,x_c)==1 .and. profit_val-b_val+theta*(1.0d0-delta)*k_val>=0.0d0) then
                    ! If no forced liquidation
                    ! Set lower and upper bound for maximization over k' in (25)
                    ! b'>=b_grid(1)
                    !kp_lb = q*b_grid_k(1)+profit_val-b_val+k_val
                
                    ! Create the RHS of Bellman eq. 25 for each k' \in k_grid
                    rhs_vec = -100000.0d0
                    do kp_c = 1,nk
                        kp_val    = k_grid(kp_c)
                        ! TODO: this b_grid should be next period's b_grid
                        b_grid_kp = b_grid(kp_c,:)
                        bprime = max(b_grid_kp(1), (1.0d0/q)*(b_val-profit_val+fun_adjcost(kp_val,k_val,theta,delta)))
                        if (kp_val<=kp_ub) then 
                            v0_int = linint(b_grid_kp,EVx(kp_c,:),bprime)
                            rhs_vec(kp_c) = q*(psi*(theta*(1.0d0-delta)*kp_val-bprime)+(1.0d0-psi)*v0_int)
                        else
                            exit    
                        endif
                    enddo !kp_c
                    max_ind = maxloc(rhs_vec,dim=1)
                    pol_kp(k_c,b_c,x_c)    = max_ind
                    val_c_new(k_c,b_c,x_c) = rhs_vec(max_ind)
                else
                    ! forced liquidation: val_c is irrelevant
                    val_c_new(k_c,b_c,x_c) = theta*(1.0d0-delta)*k_val-b_val
                endif
            enddo !k
        enddo !b
    enddo !x
	!$omp enddo
  	!$omp end parallel
  	
	! Now do Howard

	! Howard acceleration
	if (do_howard==1) then
		call sub_howard(val_c_new,pol_kp,val0_u,kp_bar,B_hat,profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,n_howard,nk,nb,nx)
	endif

end subroutine sub_vfi_onestep   
!============================================================================! 

subroutine sub_howard(val_c_new,pol_kp_ind,&
    val0_u,kp_bar,B_hat,profit_mat,k_grid,b_grid,pi_x,theta,q,delta,psi,n_howard,nk,nb,nx)
    
    use mymod, only: linint, fun_adjcost, find_loc
    implicit none
    ! Declare input variables
    integer(4), intent(in) :: nk,nb,nx
    real(8), intent(in)    :: val0_u(nk,nb,nx), kp_bar(nk,nb,nx)
    real(8), intent(in)    :: B_hat(nk,nx), profit_mat(nk,nx)
    real(8), intent(in)    :: k_grid(nk), b_grid(nk,nb), pi_x(nx,nx)
    real(8), intent(in)    :: theta, q, delta, psi
    integer(4), intent(in) :: n_howard
    integer(4), intent(in)   :: pol_kp_ind(nk,nb,nx)
    ! Declare output variables
    real(8), intent(inout)   :: val_c_new(nk,nb,nx)
    ! Declare local variables
    integer(4) :: k_c, b_c, x_c, xp_c, kp_c, max_ind, h_c
    logical :: liq1, liq2
    real(8) :: k_val, b_val, profit_val, kp_ub, kp_val, bprime, v0_int
    integer :: jstar
    real(8) :: omega
    real(8), allocatable :: EVx(:,:), b_grid_kp(:)
    real(8), allocatable :: val0_c(:,:,:), val0(:,:,:)
    integer, allocatable :: is_c(:,:,:), left_arr(:,:,:)
    real(8), allocatable :: omega_arr(:,:,:), bprime_arr(:,:,:)
    
    allocate(left_arr(nk,nb,nx),omega_arr(nk,nb,nx),bprime_arr(nk,nb,nx))
    do x_c = 1,nx
    	do b_c = 1,nb
    		do k_c = 1,nk
                k_val      = k_grid(k_c)
                b_val      = b_grid(k_c,b_c)
                profit_val = profit_mat(k_c,x_c)
                kp_c       = pol_kp_ind(k_c,b_c,x_c)
                kp_val     = k_grid(kp_c)
                b_grid_kp  = b_grid(kp_c,:)
    			bprime     = max(b_grid_kp(1), (1.0d0/q)*(b_val-profit_val+fun_adjcost(kp_val,k_val,theta,delta)))
    			! It seems that we don't need to store bprime in bprime_arr
                bprime_arr(k_c,b_c,x_c) = bprime
                call find_loc(jstar,omega, b_grid_kp,bprime)
                left_arr(k_c,b_c,x_c)  = jstar
                omega_arr(k_c,b_c,x_c) = omega
    		enddo
    	enddo
    enddo
    
    allocate(val0_c(nk,nb,nx))
    allocate(val0(nk,nb,nx),is_c(nk,nb,nx))
    allocate(EVx(nk,nb))
    
    do h_c = 1,n_howard
    
    ! STEP 2 - Impose liquidation, eq. (24)

    ! Value of constrained firm before exit
    val0_c = 0.0d0
    do x_c = 1,nx
        do b_c = 1,nb
            do k_c = 1,nk
                k_val = k_grid(k_c)
                b_val = b_grid(k_c,b_c)
                profit_val = profit_mat(k_c,x_c)
                liq1 = profit_val-b_val+theta*(1.0d0-delta)*k_val<0.0d0
                liq2 = val_c_new(k_c,b_c,x_c)<theta*(1.0d0-delta)*k_val-b_val
                if (liq1 .or. liq2) then
                    val0_c(k_c,b_c,x_c) = theta*(1.0d0-delta)*k_val-b_val
                else
                    val0_c(k_c,b_c,x_c) = val_c_new(k_c,b_c,x_c)
                endif
            enddo
        enddo
    enddo
    
    ! STEP 3 - Do equation (23)
    val0 = 0.0d0
    is_c = 1
    do x_c = 1,nx
        do b_c = 1,nb
            do k_c = 1,nk
                b_val = b_grid(k_c,b_c)
                if (b_val<=B_hat(k_c,x_c)) then
                    val0(k_c,b_c,x_c) = val0_u(k_c,b_c,x_c)
                    is_c(k_c,b_c,x_c) = 0
                else
                    val0(k_c,b_c,x_c) = val0_c(k_c,b_c,x_c)
                endif
            enddo
        enddo
    enddo
    
    ! STEP 4 - Solve eq. (25)
    
    do x_c = 1,nx
        !disp(x_c)
        ! Compute expected value in eq. (25)
        EVx = 0.0d0 !(k',b')
        do xp_c = 1,nx
            EVx = EVx+val0(:,:,xp_c)*pi_x(x_c,xp_c)
        enddo
        do b_c = 1,nb
            do k_c = 1,nk
                k_val = k_grid(k_c)
                b_val = b_grid(k_c,b_c)
                profit_val = profit_mat(k_c,x_c)
                !kp_ub = kp_bar(k_c,b_c,x_c)

                if (is_c(k_c,b_c,x_c)==1 .and. profit_val-b_val+theta*(1.0d0-delta)*k_val>=0.0d0) then
                    jstar = left_arr(k_c,b_c,x_c)
                    omega = omega_arr(k_c,b_c,x_c)
					kp_c       = pol_kp_ind(k_c,b_c,x_c)
					kp_val     = k_grid(kp_c)
					bprime     = bprime_arr(k_c,b_c,x_c)
                    !v0_int    = linint(b_grid_kp,EVx(kp_c,:),bprime)  
                    ! I precomputed the interpolation weights for speed
                    v0_int = omega*EVx(kp_c,jstar)+(1.0d0-omega)*EVx(kp_c,jstar+1)
                    !val_c_new(k_c,b_c,x_c) = q*v0_int
                    val_c_new(k_c,b_c,x_c) = q*(psi*(theta*(1.0d0-delta)*kp_val-bprime)+(1.0d0-psi)*v0_int)
                else
                    ! forced liquidation: val_c is irrelevant
                    val_c_new(k_c,b_c,x_c) = theta*(1.0d0-delta)*k_val-b_val
                endif
            enddo !k
        enddo !b
    enddo !x
    	
    ! Update	
    	
    enddo ! close howard iterations h_c
    	
end subroutine sub_howard   
!============================================================================! 


   
!program main
!implicit none
!
!write(*,*) "sub_vfi_onestep: Hello, world!"
!pause
!
!end program main
