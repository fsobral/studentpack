program studentpack

  use packmod
  
  implicit none

  ! PARAMETERS
  character(len=15) :: LTEXSOL = 'solution.tex', &
                       JSONSOL = 'solution.json'
  integer           :: MAXMEM = 5
  
  ! LOCAL SCALARS
  logical      :: checkder
  integer      :: allocerr,hnnzmax,inform,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm

  ! LOCAL ARRAYS
  character(len=80)     :: specfnm,outputfnm,vparam(10)
  logical               :: coded(11)
  logical,      pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:),xb(:,:),maxmindist(:)
  real(kind=8), allocatable :: tmpx(:)
  integer     , pointer :: seed(:),btrial(:)

  interface
     subroutine findgnite(tmpx,W,H,ntrials,ssize,seed,MINDIST,ptype)
       implicit none
       ! SCALAR ARGUMENTS
       integer, intent(in) :: ntrials,ptype,ssize
       real(kind=8), intent(in) :: H,MINDIST,W
       !
       ! ARRAY ARGUMENTS
       real(kind=8), allocatable, intent(inout) :: tmpx(:)
       integer, intent(inout) :: seed(ssize)
       !
       ! LOCAL SCALARS
       logical      :: checkder
       integer      :: allocerr,hnnzmax,inform,jcnnzmax,m,n,nvparam
       real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
            f,nlpsupn,snorm
     end subroutine findgnite
  end interface
  
  ! LOCAL SCALARS
  integer :: i,j,ntrial,ntrials,ptype,status,ssize,nmem
  real(kind=8) :: H,vover,valoc,W

  ! FUNCTIONS
  real(kind=8) :: maxaloc,minover,overlap
  
  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalfu,myevalgu,myevalhu, &
       myevalc,myevaljac,myevalhc,myevalfc,myevalgjac,myevalgjacp,&
       myevalhl,myevalhlp

  ! WRITE PRESENTATION

  write(*,7000)

  ! READ PROBLEM DATA

  ndim = 2

  write(*,*) 'Enter the minimum distance between students:'
  read(*,*) MINDIST

  write(*,*) 'Enter the approximate size of the chair (W H):'
  read(*,*) cW,cH

  write(*,*) 'Enter the dimensions of the room (W H): '
  read(*,*) W,H

  write(*,*) 'Enter the number of trials:'
  read(*,*) ntrials

  write(*,*) 'Enter the number of fixed points'
  read(*,*) nfix

  allocate(fcoord(2,nfix),frad(nfix),stat=status)
  if ( status .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  write(*,*) 'Enter the coordinates and radius of the fixed points:'
  do j = 1,nfix
     read(*,*) fcoord(1,j),fcoord(2,j),frad(j)
  end do

6000  write(*,*) 'Problem type: 1- Max radius 2- Max chairs 3- Max radius with rows 4- Max chairs with rows '
  read(*,*) ptype

  if ( ptype .ne. 1 .and. ptype .ne. 2 .and. &
       ptype .ne. 3 .and. ptype .ne. 4 ) goto 6000     

  if ( ptype .eq. 1 .or. ptype .eq. 3 ) then
     write(*,*) 'Enter the number of chairs: '
     read(*,*) nite
  end if
  
  ! INITIALIZE REGIONS

  nregw = CEILING(W / MINDIST)
  nregh = CEILING(H / MINDIST)

  call RANDOM_SEED(SIZE=ssize)

  allocate(start(0:nregw+1,0:nregh+1),br(2,nregw*nregh), &
           seed(ssize),stat=status)

  if ( status .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  start = 0
  br    = 0
  nbr   = 0

  ! If problem type is 2, try to find the maximum number of chairs,
  ! respecting the minimum distance.
  if ( ptype .eq. 2 .or. ptype .eq. 4 ) then
     call findgnite(tmpx,W,H,ntrials,ssize,seed,MINDIST,ptype)
  end if

  ! Finish allocating the structure of regions
  
  allocate(diagb(ndim,ndim,nite),next(nite))
  
  if ( status .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  ! Number of variables

  n = 2 * nite + 1

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),xb(n,MAXMEM),btrial(MAXMEM), &
           maxmindist(MAXMEM),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  ! Rows

  if ( ptype .eq. 4 ) then
     xb(1:n, 1) = tmpx
     nmem = 1

     ! We need to deallocate tmpx that was used for findgnite with
     ! ptype = 4.
     deallocate(tmpx,stat=allocerr)
     if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error in main program.'
        stop
     end if
     
     goto 6001
  end if
  
  if ( ptype .eq. 3 ) then
     call generate_x(xb(1:n,1), n, W, H, cH, cW, nite)
     nmem = 1
     goto 6001
  end if
  
  ! We assume that the right part of the room  is NOT a wall
  do i = 1,nite
     l(2 * i - 1) = cW / 2.0D0
     u(2 * i - 1) = W
     l(2 * i)     = cH / 2.0D0
     u(2 * i)     = H - cH / 2.0D0
  end do

  l(n) = max(MINDIST, sqrt(cW **2 + cH ** 2) / 2.0D0)
  u(n) = 1.5D0 * sqrt(W **2 + H ** 2)

  ! Constraints

  m = 0

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  ! Coded subroutines

  coded(1:3)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
  coded(4:11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = n
  hnnzmax  = ((nite + nfix) * ((nite + nfix) + 1) / 2) * &
       ndim * (ndim + 3) + &
       (nite + nfix) * (ndim * (ndim + 1) / 2) + 1

  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting

  epsfeas   = 1.0d-02
  epsopt    = 1.0d-02

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  outputfnm = ''
  specfnm   = ''

  vparam(1) = 'LARGEST-PENALTY-PARAMETER-ALLOWED 1.0D+50'
  vparam(2) = 'PENALTY-PARAMETER-INITIAL-VALUE 1.0D+1'
  vparam(3) = 'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED'
  
  nvparam = 3

  ! OPTIMIZE THE BOX-CONSTRAINED PENALIZED PROBLEM

  btrial     = 0
  maxmindist = 0.0D0
  nmem       = 0

  do ntrial = 1,ntrials
  
     seed = 123456789.0D0 + ntrial
     call RANDOM_SEED(PUT=seed)
  
     call RANDOM_NUMBER(x)
     do i = 1, n
        x(i) = l(i) + x(i) * (u(i) - l(i))
     end do
  
     call algencan(myevalfu,myevalgu,myevalhu,myevalc,myevaljac,myevalhc, &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
       checkder,f,cnorm,snorm,nlpsupn,inform)

     ! TEST SOLUTION
     vover = minover(n,x)
     valoc = maxaloc(n,x,l,u)

     ! Sugestao do Felipe
     x(n) = vover
     
     if ( vover .ge. MINDIST .and. valoc .le. ERR) then
        call packsort(n,x,nite,ndim)
        do i = 1,nmem
           if ( ( vover .gt. maxmindist(i) + ERR ) .or. &
                ( ABS(vover - maxmindist(i)) .le. ERR .and. &
                  MAXVAL(ABS(x(1:n) - xb(1:n,i))) .gt. ERR ) ) then

              if ( nmem .lt. MAXMEM ) nmem = nmem + 1
              
              do j = nmem,i + 1,-1
                 xb(:,j) = xb(:,j - 1)
                 maxmindist(j) = maxmindist(j - 1)
                 btrial(j) = btrial(j - 1)
              end do
              xb(:,i) = x
              maxmindist(i) = vover
              btrial(i) = ntrial

              exit

           end if
        end do
        if ( nmem .eq. 0 ) then
           xb(:,1)       = x
           maxmindist(1) = vover
           btrial(1)     = ntrial
           nmem          = nmem + 1
        end if
     end if

     write(*,8020) ntrial,maxmindist(1),vover,x(n),valoc
  
  end do

  if ( nmem .eq. 0 ) then
     call tojson(n,nmem,xb,nite,W,H,JSONSOL,.false.)
     stop
  end if
  
  ! RUN CONSTRAINED PROBLEM
  
  deallocate(equatn,linear,lambda)

  m = 1

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  coded(1:6)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
  coded(7:11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub

  ! Parameters setting

  epsfeas   = 1.0d-08
  epsopt    = 1.0d-08

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  do j = 1,nmem
  
     equatn(1) = .true.
     lambda(1) = PEN * overlap(n,xb(1:n,j))
     linear(1) = .false.

     ! Coded subroutines

     ! seed = 123456789.0D0 + btrial(j)
     ! call RANDOM_SEED(PUT=seed)
  
     ! call RANDOM_NUMBER(x)
     ! do i = 1, n
     !    x(i) = xb(i,j)
     ! end do

     call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,xb(1:n,j),l,u,m,lambda,equatn,linear,coded, &
       checkder,f,cnorm,snorm,nlpsupn,inform)

     ! TEST SOLUTION
     vover = minover(n,xb(1:n,j))
     valoc = maxaloc(n,xb(1:n,j),l,u)

     ! Sugestao do Felipe
     xb(n,j) = vover
     
     write(*,8020) btrial(j),maxmindist(j),vover,xb(n,j),valoc

  end do

  ! Free Algencan structure
  deallocate(lambda,equatn,linear,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

6001 continue
  
  call drawsol(nite,W,H,n,nmem,xb(1:n,1:nmem),LTEXSOL)
  
  if (xb(n,1) .lt. mindist) then
  	  call tojson(n,nmem,xb(1:n,1:nmem),nite,W,H,JSONSOL,.false.)
  else
  	  call tojson(n,nmem,xb(1:n,1:nmem),nite,W,H,JSONSOL,.true.)
  end if


  ! Free structures
  
  deallocate(x,l,u,xb,maxmindist,btrial, stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

  deallocate(start,next,br,diagb,seed,fcoord,frad,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

  stop

  ! NON-EXECUTABLE STATEMENTS

7000 format(/,' Welcome to STUDENTPACK !!!',/                              &
       /,' This program uses ALGENCAN and the strategy',          &
       ' described in:',/                                       &
       /,' to appear.',/)
8020 format(' Trial = ',I4,' Min Dist (best = ',F12.8,') = ',F12.8, &
       ' Fobj = ',1P,D9.1,' Feasibility = ',1P,D9.1)

end program studentpack

!     ******************************************************
!     ******************************************************

subroutine packsort(n,x,nite,ndim)

  use packmod, only: ERR
  
  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n,nite,ndim

  ! ARRAY ARGUMENTS
  real(kind=8), intent(inout) :: x(n)

  ! LOCAL SCALARS
  integer      :: i,ind1,ind2,j,k,t
  real(kind=8) :: tmp

  do i = 1,nite
     t = i
     ! Look for the 'smallest' item
     do j = i + 1, nite
        do k = 1,ndim
           ind1 = (t - 1) * ndim + k
           ind2 = (j - 1) * ndim + k

           tmp = x(ind2) - x(ind1)

           if ( tmp .lt. - ERR ) then
              t = j
              exit
           end if
           if ( abs(tmp) .le. -ERR ) continue
           if ( tmp .gt. ERR ) exit
           
        end do
     end do
     ! Swap items
     if ( t .ne. i ) then
        do k = 1,ndim
           ind1 = (i - 1) * ndim + k
           ind2 = (t - 1) * ndim + k
           !
           tmp = x(ind1)
           x(ind1) = x(ind2)
           x(ind2) = tmp
        end do
     end if
  end do
  
end subroutine packsort

!     ******************************************************
!     ******************************************************

subroutine findgnite(tmpx,W,H,ntrials,ssize,seed,MINDIST,ptype)

  ! This subroutine find the greatest value of the number of items
  ! that can be packed inside the given region. The value 'nite' is
  ! updated in the global scalar 'nite' in module 'packmod'.
  
  use packmod, only: ndim,nfix,nite,cH,cW,ERR,next,frad,diagb

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: ntrials,ptype,ssize
  real(kind=8), intent(in) :: H,MINDIST,W

  ! ARRAY ARGUMENTS
  real(kind=8), allocatable, intent(inout) :: tmpx(:)
  integer, intent(inout) :: seed(ssize)
  
  ! LOCAL SCALARS
  logical      :: checkder
  integer      :: allocerr,hnnzmax,inform,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm

  ! LOCAL ARRAYS
  character(len=80)     :: specfnm,outputfnm,vparam(10)
  logical               :: coded(11)
  logical,      pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! LOCAL SCALARS
  integer      :: i,lnite,unite,ntrial
  real(kind=8) :: farea,minover,vover,mxaloc,maxaloc

  ! EXTERNAL SUBROUTINES
  external :: myevalfu,myevalgu,myevalhu,myevalc,myevaljac,myevalhc,&
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! We assume that the lower bound is feasible
  farea = W * H
  do i = 1,nfix
     farea = farea - ACOS(-1.0D0) * frad(i) ** 2
  end do
  
  unite = CEILING(W * H / (ACOS(-1.0D0) * (MINDIST / 3.0D0) ** 2))
  lnite = INT(farea / (ACOS(-1.0D0) * (MINDIST) ** 2))
  
8000 nite = INT((lnite + unite) / 2)

  allocate(diagb(ndim,ndim,nite),next(nite))
  
  if ( unite - lnite .le. 1 ) goto 9000

  ! Number of variables
  
  n = 2 * nite + 1

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in subroutine findgnite.'
     stop
  end if

  ! Rows
  if ( ptype .eq. 4 ) then

     call generate_x(x, n, W, H, cH, cW, nite)

     if ( x(n) .lt. MINDIST ) then

        write(*,9050) nite, 0
        unite = nite
        
     else

        write(*,9051) nite,0,0.0D0,0.0D0
        lnite = nite

        ! Save successful solution
        if ( ALLOCATED(tmpx) ) then
           deallocate(tmpx,stat=allocerr)
           if ( allocerr .ne. 0 ) then
              write(*,*) 'Deallocation error in subroutine findgnite.'
              stop
           end if
        end if
        allocate(tmpx(n),stat=allocerr)
        if ( allocerr .ne. 0 ) then
           write(*,*) 'Allocation error in subroutine findgnite.'
           stop
        end if
        tmpx = x
        
     end if

     goto 8100
  end if
  
  ! Constraints

  m = 0

  allocate(l(n),u(n),equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in subroutine findgnite.'
     stop
  end if

  ! We assume that the right part of the room IS NOT a wall
  do i = 1,nite
     l(2 * i - 1) = cW / 2.0D0
     u(2 * i - 1) = W
     l(2 * i)     = cH / 2.0D0
     u(2 * i)     = H - cH / 2.0D0
  end do

  l(n) = max(MINDIST, sqrt(cW **2 + cH ** 2) / 2.0D0)
  u(n) = 1.5D0 * sqrt(W **2 + H ** 2)

  ! Coded subroutines

  coded(1:3)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
  coded(4:11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 0
  hnnzmax  = ((nite + nfix) * ((nite + nfix) + 1) / 2) * &
       ndim * (ndim + 3) + &
       (nite + nfix) * (ndim * (ndim + 1) / 2) + 1

  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting

  epsfeas   = 1.0d-02
  epsopt    = 1.0d-02

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  outputfnm = ''
  specfnm   = ''

  vparam(1) = 'LARGEST-PENALTY-PARAMETER-ALLOWED 1.0D+50'
  vparam(2) = 'PENALTY-PARAMETER-INITIAL-VALUE 1.0D+1'
  vparam(3) = 'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED'
  
  nvparam   = 3

  ! OPTIMIZE THE BOX-CONSTRAINED PENALIZED PROBLEM

  do ntrial = 1,ntrials

     seed = 123456789.0D0 + ntrial
     call RANDOM_SEED(PUT=seed)

     call RANDOM_NUMBER(x)
     do i = 1, n
        x(i) = l(i) + x(i) * (u(i) - l(i))
     end do

     call algencan(myevalfu,myevalgu,myevalhu,myevalc,myevaljac, &
     myevalhc,myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,&
     jcnnzmax,hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,&
     outputfnm,specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,   &
     linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)
     
     ! TEST SOLUTION
     vover = minover(n,x)
     mxaloc = maxaloc(n,x,l,u)

     if ( vover .ge. MINDIST .and. mxaloc .le. ERR ) exit

  end do

  deallocate(l,u,lambda,equatn,linear,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in subroutine findgnite.'
     stop
  end if

  if ( ntrial .gt. ntrials ) then

     write(*,9050) nite,ntrials
     unite = nite

  else

     write(*,9051) nite,ntrial,vover,mxaloc
     lnite = nite

  end if

8100 deallocate(x,diagb,next,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in subroutine findgnite.'
     stop
  end if
  
  goto 8000

9000 return

9050 format('No packing for ',I4,' chairs. Trials: ',I5)
9051 format('Found packing for ',I4,' chairs. Trials: ',I5,&
            ' Dist: ', 1P,D9.1,' Feas: ',1P,D9.1)

end subroutine findgnite
  
!     ******************************************************************
!     ******************************************************************

real(kind=8) function minover(n,x)

  use packmod, only: nbr,start,next,br

  implicit none
  
  ! SCALAR ARGUMENTS
  integer :: n

  !     ARRAY ARGUMENTS
  real(kind=8) :: x(n)

  !     LOCAL SCALARS
  integer :: i,j,k,l
  real(kind=8) :: mnrover

  call classify(n,x)

  minover = 1.0D+20
  do k = 1,nbr
     i = br(1,k)
     j = br(2,k)
     l = start(i,j)
10   if ( l .ne. 0 ) then
        minover = min( minover, mnrover(n,x,l,next(l)) )
        minover = min( minover, mnrover(n,x,l,start(i-1,j+1)) )
        minover = min( minover, mnrover(n,x,l,start(i,  j+1)) )
        minover = min( minover, mnrover(n,x,l,start(i+1,j+1)) )
        minover = min( minover, mnrover(n,x,l,start(i+1,j))   )
        l = next(l)
        go to 10
     end if
  end do

end function minover

!     ******************************************************************
!     ******************************************************************

real(kind=8) function mnrover(n,x,i,lstart)

  use packmod, only: ndim,next

  implicit none
  
  !     SCALAR ARGUMENTS
  integer :: i,lstart,n

  !     ARRAY ARGUMENTS
  real(kind=8) :: x(n)

  !     LOCAL SCALARS
  integer :: ind1,ind2,j,k
  real(kind=8) :: dist

  mnrover = 1.0D+20

  j = lstart

10 if ( j .ne. 0 ) then

     dist = 0.0d0
     do k = 1,ndim
        ind1 = ( i - 1 ) * ndim + k
        ind2 = ( j - 1 ) * ndim + k
        dist = dist + ( x(ind1) - x(ind2) ) ** 2
     end do
     dist = sqrt( dist )

     mnrover = min( mnrover,  dist )

     j = next(j)
     go to 10

  end if

end function mnrover

!     ******************************************************
!     ******************************************************

real(kind=8) function maxaloc(n,x,l,u)

  use packmod, only: ndim,nfix,nite,cH,cW,fcoord,frad
  
  implicit none
  
  !     SCALAR ARGUMENTS
  integer :: n

  !     ARRAY ARGUMENTS
  real(kind=8) :: l(n),u(n),x(n)

  !     LOCAL SCALARS
  integer :: i,j,k,ind
  real(kind=8) :: alloc,tmp

  maxaloc = 0.0d0

  do i = 1,n - 1

     alloc   = max( l(i) - x(i), x(i) - u(i) )
     maxaloc = max( maxaloc, alloc )

  end do

  ! Check feasibility with respect the forbidden circles

  do i = 1,nfix
     tmp = frad(i) + sqrt(cH ** 2 + cW ** 2) / 2.0D0
     do j = 1,nite
        alloc = 0.0D0
        do k = 1,ndim
           ind = (j - 1) * ndim + k
           alloc = alloc + (x(ind) - fcoord(k,i)) ** 2
        end do
        maxaloc = max( maxaloc, tmp - sqrt(alloc) ) 
     end do
  end do

end function maxaloc

!     ******************************************************
!     ******************************************************

subroutine drawsol(nite,W,H,n,nmem,x,solfile)

  use packmod, only: nfix,fcoord,frad

  implicit none
  
  !     This subroutine generate a metapost file with the
  !     graphical representation of the problem.

  !     SCALAR ARGUMENTS
  integer, intent(in)      :: n,nite,nmem
  real(kind=8), intent(in) :: W,H

  !     ARRAY ARGUMENTS
  real(kind=8), intent(in)      :: x(n,nmem)
  character(len=12), intent(in) :: solfile

  !     LOCAL SCALARS
  integer :: i,j
  real(kind=8) :: scale,radius

  open(unit=10,file=solfile)

  write(10,01)

  do j = 1,nmem

     if ( j .gt. 1 ) write(10,12)

     radius = x(n,j) / 2.0D0
	 
     ! SCALING
     scale = min(20.0D0 / (H + 2.0D0) * radius, &
                 10.0D0 / (W + 2.0D0 * radius))
     write(10,10) 2.0D0*radius,nite,scale,scale

     ! CLASSROOM
     write(10,25) W,H
     ! PROFESSOR
     write(10,26) W+1.0D0,W+5.0,H


     ! CIRCULAR ITEMS
     do i = 1,nite
        write(10,20) x(2*i-1,j),x(2*i,j),i,x(2*i-1,j),x(2*i,j),radius
     end do

     ! FIXED ITEMS
     do i = 1,nfix
        write(10,21) fcoord(1,i),fcoord(2,i),fcoord(1,i),&
                     fcoord(2,i),frad(i)
     end do

     write(10,30)

  end do

  write(10,31)

  close(10)

  ! NON-EXECUTABLE STATEMENTS

01 format('\documentclass{article}',/,'\usepackage{tikz}',/, &
          '\usepackage[top=1cm,bottom=1cm,left=1cm,right=1cm]{geometry}',/, &
          '\begin{document}',/)
10 format('\begin{flushleft} \LARGE Minimum distance: ',&
        F20.10,'\\ Number of chairs: ',I5,'\end{flushleft}',/,&
        '\begin{tikzpicture}[scale=',F6.3, &
        ', every node/.style={scale=',F6.3,'}]')
12 format(/,'\newpage',/)
20 format('\filldraw (',F20.10,',',F20.10,') circle (1pt) node[above right] {',I4,'};',/, &
        '\filldraw[fill=blue!40!white, fill opacity=0.2, draw=black] (', &
        F20.10,',',F20.10,') circle (',F20.10,'cm);')
21 format('\filldraw[red] (',F20.10,',',F20.10,') circle (1pt);',/, &
        '\filldraw[fill=red!40!white, fill opacity=0.2, draw=black] (', &
        F20.10,',',F20.10,') circle (',F20.10,'cm);')
25 format('\draw (0,0) rectangle (',F20.10,',',F20.10,');')
26 format('\draw[draw=red] (',F20.10,',0.0) rectangle (',F20.10,',',F20.10,');')
30 format('\end{tikzpicture}')
31 format('\end{document}') 

end subroutine drawsol

!     ******************************************************
!     ******************************************************

subroutine tojson(n,nmem,x,nite,W,H,solfile,foundsol)

  implicit none
  
  ! This subroutine generate a JSON file with the solution

  ! SCALAR ARGUMENTS
  integer, intent(in)      :: n,nite,nmem
  logical, intent(in)      :: foundsol
  real(kind=8), intent(in) :: W,H

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)      :: x(n,nmem)
  character(len=15), intent(in) :: solfile

  ! LOCAL SCALARS
  integer :: i,j

  open(unit=10,file=solfile)

  if ( .not. foundsol ) then
     write(10,22)
     goto 99
  end if

  ! CLASSROOM
  write(10,10) nite,x(n,1),W,H

  do j = 1,nmem
     
     ! CIRCULAR ITEMS
     write(10,11) x(n,j)
     do i = 1,nite - 1
        write(10,20) x(2*i-1,j),x(2*i,j)
     end do
     write(10,21) x(2*nite-1,j),x(2*nite,j)

     if ( j .lt. nmem ) write(10,23)
     if ( j .eq. nmem ) write(10,24)

  end do

99 close(10)

  !     NON-EXECUTABLE STATEMENTS

10 format('{',/,'  "found_solution": true,',/, &
        '  "number_items": ',I5,',',/, &
        '  "min_distance": ',F20.10,',',/, &
        '  "width": ',F20.10,',',/, &
        '  "height": ',F20.10,',',/, &
        '  "solutions": [')
11 format('    {"min_distance": ',F20.10,',',/, &
          '     "positions"   : [ ')
20 format('                       [',F20.10,',',F20.10,'],')
21 format('                       [',F20.10,',',F20.10,']')
23 format('                     ]},')
24 format('                     ]}',/,']}')
22 format('{',/,'   "found_solution": false',/,'}')

end subroutine tojson

! ******************************************************************
! ******************************************************************

subroutine myevalfu(n,x,f,flag)

  use packmod, only: PEN
  
  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! LOCAL SCALARS
  real(kind=8) :: overlap

  flag = 0

  f = - x(n) + PEN * overlap(n,x)

end subroutine myevalfu

! ******************************************************************
! ******************************************************************

subroutine myevalgu(n,x,g,flag)

  use packmod, only: PEN
  
  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! LOCAL SCALARS
  integer :: i

  flag = 0

  call goverlap(n,x,g)

  do i = 1,n - 1
     g(i) = g(i) * PEN
  end do
  
  g(n) = - 1.0D0 + PEN * g(n)

end subroutine myevalgu

! ******************************************************************
! ******************************************************************

subroutine myevalhu(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  use packmod, only: PEN
  
  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,n
  integer,      intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  ! LOCAL SCALARS
  integer :: i

  flag = 0
  lmem = .false.

  call hoverlap(n,x,hrow,hcol,hval,hnnz,lim,lmem)

  do i = 1,hnnz
     hval(i) = hval(i) * PEN
  end do

end subroutine myevalhu

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  flag = 0

  f = - x(n)

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  flag = 0

  g    =   0.0D0
  g(n) = - 1.0D0

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,n
  integer,      intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = 0
  lmem = .false.

  hnnz = 0

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: ind,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)

  ! LOCAL SCALARS
  real(kind=8) :: overlap

  flag = 0

  if ( ind .eq. 1 ) then
     c = overlap(n,x)
  else
     flag = - 1
  end if

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  ! LOCAL SCALARS
  integer :: i

  flag = 0
  lmem = .false.
  
  if ( ind .eq. 1 ) then

     jcnnz = n

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     call goverlap(n,x,jcval)
     do i = 1,jcnnz
        jcvar(i) = i
     end do

  else

     flag = - 1

  end if

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  flag = 0
  lmem = .false.

  if ( ind .eq. 1 ) then
     call hoverlap(n,x,hcrow,hccol,hcval,hcnnz,lim,lmem)
  else
     flag = - 1
  end if

end subroutine myevalhc

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: m,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: c(m)

  flag = - 1

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  flag = - 1

end subroutine myevalgjac

! ******************************************************************
! ******************************************************************

subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,   intent(inout) :: gotj
  integer,   intent(in)    :: m,n
  integer,   intent(out)   :: flag
  character, intent(in)    :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)    :: x(n)
  real(kind=8), intent(inout) :: p(m),q(n)
  real(kind=8), intent(out)   :: g(n)

  flag = - 1

end subroutine myevalgjacp

! ******************************************************************
! ******************************************************************

subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,hlnnz
  real(kind=8), intent(in)  :: sf

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hlcol(lim),hlrow(lim)
  real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
  real(kind=8), intent(out) :: hlval(lim)

  flag = - 1

end subroutine myevalhl

! ******************************************************************
! ******************************************************************

subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(inout) :: goth
  integer,      intent(in)    :: m,n
  integer,      intent(out)   :: flag
  real(kind=8), intent(in)    :: sf

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
  real(kind=8), intent(out) :: hp(n)

  flag = - 1

end subroutine myevalhlp

!     ******************************************************************
!     ******************************************************************

real(kind=8) function overlap(n,x)

  use packmod, only: nbr,nite,ndim,nfix,br,next,start,fcoord,frad,cH,cW

  implicit none
  
  !     SCALAR ARGUMENTS
  integer :: n

  !     ARRAY ARGUMENTS
  real(kind=8) :: x(n)

  !     LOCAL SCALARS
  integer :: i,j,k,l,ind
  real(8) :: dist,nrdist,tmp

  !     CLASSIFY ITEMS BY REGIONS
  call classify(n,x)

  !     COMPUTE SPARSE OVERLAPPING
  overlap = 0.0d0
  do k = 1,nbr
     i = br(1,k)
     j = br(2,k)
     l = start(i,j)
10   if ( l .ne. 0 ) then
        overlap =                                     &
             overlap + nrdist(n,x,l,next(l))        + &
                       nrdist(n,x,l,start(i-1,j+1)) + &
                       nrdist(n,x,l,start(i,  j+1)) + &
                       nrdist(n,x,l,start(i+1,j+1)) + &
                       nrdist(n,x,l,start(i+1,j))
        l = next(l)
        go to 10
     end if
  end do

  ! FIXED ITEMS
  do j = 1,nfix
     tmp = frad(j) + sqrt(cH ** 2 + cW ** 2) / 2.0D0
     do i = 1,nite
        dist = 0.0d0
        do k = 1,ndim
           ind  = ( i - 1 ) * ndim + k
           dist = dist + ( x(ind) - fcoord(k,j) ) ** 2
        end do
        dist = max( 0.0d0, tmp ** 2 - dist )
        
        overlap = overlap + dist ** 2
     end do
  end do

end function overlap

!     ******************************************************************
!     ******************************************************************

subroutine goverlap(n,x,g)

  use packmod, only: nbr,nite,ndim,nfix,br,next,start,fcoord,frad,cH,cW

  implicit none

  !     SCALAR ARGUMENTS
  integer, intent(in) :: n

  !     ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  !     LOCAL SCALARS
  integer      :: i,j,k,l,ind
  real(kind=8) :: dist,tmp

  !     CLASSIFY ITEMS BY REGIONS
  call classify(n,x)

  !     COMPUTE SPARSE OVERLAPPING

  g = 0.0D0

  do k = 1,nbr
     i = br(1,k)
     j = br(2,k)
     l = start(i,j)
10   if ( l .ne. 0 ) then
        call gnrdist(n,x,g,l,next(l)) 
        call gnrdist(n,x,g,l,start(i-1,j+1))
        call gnrdist(n,x,g,l,start(i,  j+1))  
        call gnrdist(n,x,g,l,start(i+1,j+1))  
        call gnrdist(n,x,g,l,start(i+1,j  ))
        l = next(l)
        go to 10
     end if
  end do

  ! FIXED ITEMS
  do j = 1,nfix
     tmp = frad(j) + sqrt(cH ** 2 + cW ** 2) / 2.0D0
     do i = 1,nite
        
        dist = 0.0d0
        do k = 1,ndim
           ind  = ( i - 1 ) * ndim + k
           dist = dist + ( x(ind) - fcoord(k,j) ) ** 2
        end do

        dist = max( 0.0d0, tmp ** 2 - dist )

        if ( dist .ne. 0.0D0 ) then
           do k = 1,ndim
              ind = (i - 1) * ndim + k
              g(ind) = g(ind) - 4.0D0 * dist * (x(ind) - fcoord(k,j))
           end do
        end if
        
     end do
  end do

end subroutine goverlap

!     ******************************************************************
!     ******************************************************************

subroutine hoverlap(n,x,hlin,hcol,hval,nnzh,lim,lmem)

  use packmod, only: nbr,ndim,nite,nfix,dd,br,diagb,next,start,&
                     fcoord,frad,cH,cW

  implicit none
  
  !     SCALAR ARGUMENTS
  integer, intent(in)  :: n,lim
  integer, intent(out) :: nnzh
  logical, intent(out) :: lmem

  !     ARRAY ARGUMENTS
  integer, intent(out)      :: hcol(lim),hlin(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  !     LOCAL SCALARS
  integer      :: i,j,k,l,ind1,ind2
  real(kind=8) :: dist,tmp,tmpr

  !     CLASSIFY ITEMS BY REGIONS
  call classify(n,x)

  !     INITALIZE DIAGONAL BLOCKS
  diagb = 0.0D0
  dd    = 0.0D0

  !     COMPUTE SPARSE OVERLAPPING SECOND DERIVATIVES

  nnzh = 0

  do k = 1,nbr
     i = br(1,k)
     j = br(2,k)
     l = start(i,j)
10   if ( l .ne. 0 ) then
        call hnrdist(n,x,hlin,hcol,hval,nnzh,l,next(l),lim) 
        call hnrdist(n,x,hlin,hcol,hval,nnzh,l,start(i-1,j+1),lim)
        call hnrdist(n,x,hlin,hcol,hval,nnzh,l,start(i,  j+1),lim)  
        call hnrdist(n,x,hlin,hcol,hval,nnzh,l,start(i+1,j+1),lim)  
        call hnrdist(n,x,hlin,hcol,hval,nnzh,l,start(i+1,j  ),lim)
        l = next(l)
        go to 10
     end if
  end do

  ! FIXED ITEMS
  do j = 1,nfix
     tmpr = frad(j) + sqrt(cH ** 2 + cW ** 2) / 2.0D0
     do i = 1,nite
        dist = 0.0d0
        do k = 1,ndim
           ind1  = ( i - 1 ) * ndim + k
           dist = dist + ( x(ind1) - fcoord(k,j) ) ** 2
        end do

        if ( dist .le. tmpr ** 2 ) then
           do k = 1,ndim
              ind1 = (i - 1) * ndim + k

              tmp = 8.0d0 * ( x(ind1) - fcoord(k,j) ) ** 2 &
                   - 4.0d0 * ( tmpr ** 2 - dist )
              diagb(k,k,i) = diagb(k,k,i) + tmp

              do l = k + 1,ndim
                 ind2 = (i - 1) * ndim + l
                 tmp = 8.0D0 * (x(ind1) - fcoord(k,j)) * &
                      (x(ind2) - fcoord(l,j))
                 nnzh = nnzh + 1
                 hlin(nnzh) = ind2
                 hcol(nnzh) = ind1
                 hval(nnzh) = tmp
              end do
           end do
        end if
        
     end do
  end do

  !     ADD DIAGONAL BLOCKS TO THE SPARSE STRUCTURE
  do k = 1,nite
     do j = 1,ndim
        do i = j,ndim
           if ( diagb(i,j,k) .ne. 0.0d0 ) then
              nnzh = nnzh + 1
              hlin(nnzh) = ( k - 1 ) * ndim + i
              hcol(nnzh) = ( k - 1 ) * ndim + j
              hval(nnzh) = diagb(i,j,k)
           end if
        end do
     end do
  end do

  ! DERIVATIVES W.R.T 'D'
  nnzh = nnzh + 1
  hlin(nnzh) = n
  hcol(nnzh) = n
  hval(nnzh) = dd

end subroutine hoverlap

!     ******************************************************************
!     ******************************************************************

real(kind=8) function nrdist(n,x,i,lstart)

  use packmod, only: ndim,next

  implicit none
  
  !     SCALAR ARGUMENTS
  integer :: i,lstart,n

  !     ARRAY ARGUMENTS
  real(kind=8) :: x(n)

  !     LOCAL SCALARS
  integer :: ind1,ind2,j,k
  real(kind=8) :: dist

  nrdist = 0.0d0

  j = lstart

10 if ( j .ne. 0 ) then

     dist = 0.0d0
     do k = 1,ndim
        ind1 = ( i - 1 ) * ndim + k
        ind2 = ( j - 1 ) * ndim + k
        dist = dist + ( x(ind1) - x(ind2) ) ** 2
     end do

     dist = max( 0.0d0, x(n) ** 2 - dist )
     nrdist = nrdist + dist ** 2

     j = next(j)
     go to 10

  end if

end function nrdist

!     ******************************************************************
!     ******************************************************************

subroutine gnrdist(n,x,g,i,lstart)

  use packmod, only: ndim,next

  implicit none

  !     SCALAR ARGUMENTS
  integer, intent(in) :: i,lstart,n

  !     ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  !     LOCAL SCALARS
  integer :: ind1,ind2,j,k
  real(kind=8) :: dist,tmp

  j = lstart

10 if ( j .ne. 0 ) then

     dist = 0.0d0
     do k = 1,ndim
        ind1 = ( i - 1 ) * ndim + k
        ind2 = ( j - 1 ) * ndim + k
        dist = dist + ( x(ind1) - x(ind2) ) ** 2
     end do

     dist = max( 0.0d0, x(n) ** 2 - dist )

     if ( dist .ne. 0.d0 ) then
        do k = 1,ndim
           ind1 = ( i - 1 ) * ndim + k
           ind2 = ( j - 1 ) * ndim + k
           tmp = 4.0d0 * dist * ( x(ind1) - x(ind2) )
           g(ind1) = g(ind1) - tmp 
           g(ind2) = g(ind2) + tmp
        end do
        g(n) = g(n) + 4.0D0 * x(n) * dist
     end if

     j = next(j)
     go to 10

  end if

end subroutine gnrdist

!     ******************************************************************
!     ******************************************************************

subroutine hnrdist(n,x,hlin,hcol,hval,nnzh,i,lstart,lim)

  use packmod, only: ndim,dd,diagb,next

  implicit none
  
  !     SCALAR ARGUMENTS
  integer, intent(in)  :: i,lstart,lim,n
  integer, intent(out) :: nnzh

  !     ARRAY ARGUMENTS
  integer, intent(out)      :: hcol(lim),hlin(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  !     LOCAL SCALARS
  integer :: col,ind1,ind2,ind3,ind4,j,k,l,lin
  real(kind=8) :: dist,tmp

  j = lstart

10 if ( j .ne. 0 ) then

     if ( j .gt. i ) then
        col = i
        lin = j
     else
        col = j
        lin = i
     end if

     dist = 0.0d0
     do k = 1,ndim
        ind1 = ( col - 1 ) * ndim + k
        ind2 = ( lin - 1 ) * ndim + k
        dist = dist + ( x(ind1) - x(ind2) ) ** 2
     end do

     if ( dist .le. x(n) ** 2 ) then
        do k = 1,ndim
           ind1 = ( col - 1 ) * ndim + k
           ind2 = ( lin - 1 ) * ndim + k
           tmp = 8.0d0 * ( x(ind1) - x(ind2) ) ** 2 &
                - 4.0d0 * ( x(n) ** 2 - dist )
           if ( tmp .ne. 0.0d0 ) then
              ! H(ind1,ind1) = H(ind1,ind1) + tmp
              diagb(k,k,col) = diagb(k,k,col) + tmp
              ! H(ind2,ind2) = H(ind2,ind2) + tmp
              diagb(k,k,lin) = diagb(k,k,lin) + tmp
              ! H(ind2,ind1) = H(ind2,ind1) - tmp
              nnzh = nnzh + 1
              hlin(nnzh) = ind2
              hcol(nnzh) = ind1
              hval(nnzh) = - tmp
           end if
           do l = 1,k - 1
              ind3 = ( col - 1 ) * ndim + l
              ind4 = ( lin - 1 ) * ndim + l
              tmp = 8.0d0 * ( x(ind3) - x(ind4) ) &
                   * ( x(ind1) - x(ind2) )
              if ( tmp .ne. 0.0d0 ) then
                 ! H(ind1,ind3) = H(ind1,ind3) + tmp
                 diagb(k,l,col) = diagb(k,l,col) + tmp
                 ! H(ind2,ind4) = H(ind2,ind4) + tmp
                 diagb(k,l,lin) = diagb(k,l,lin) + tmp
                 ! H(ind2,ind3) = H(ind2,ind3) - tmp
                 nnzh = nnzh + 1
                 hlin(nnzh) = ind2
                 hcol(nnzh) = ind3
                 hval(nnzh) = - tmp
              end if
           end do
           do l = k + 1,ndim
              ind3 = ( col - 1 ) * ndim + l
              ind4 = ( lin - 1 ) * ndim + l
              tmp = 8.0d0 * ( x(ind3) - x(ind4) ) &
                   * ( x(ind1) - x(ind2) )
              if ( tmp .ne. 0.0d0 ) then
                 ! H(ind2,ind3) = H(ind2,ind3) - tmp
                 nnzh = nnzh + 1
                 hlin(nnzh) = ind2
                 hcol(nnzh) = ind3
                 hval(nnzh) = - tmp
              end if
           end do
           ! Partial derivative w.r.t. D (x(n))
           tmp = 8.0D0 * x(n) * (x(ind1) - x(ind2))
           if ( tmp .ne. 0.0d0 ) then
              ! H(n,ind1) -= tmp
              nnzh = nnzh + 1
              hlin(nnzh) = n
              hcol(nnzh) = ind1
              hval(nnzh) = - tmp
              ! H(n,ind2) += tmp
              nnzh = nnzh + 1
              hlin(nnzh) = n
              hcol(nnzh) = ind2
              hval(nnzh) = tmp
           end if
        end do
        dd = dd + 4.0D0 * (x(n) ** 2 - dist) + 8.0D0 * x(n) ** 2
     end if

     j = next(j)
     go to 10

  end if

end subroutine hnrdist

!     ******************************************************************
!     ******************************************************************

subroutine classify(n,x)

  use packmod, only: nbr,ndim,nite,nregh,nregw,br,next,start

  implicit none

  !     SCALAR ARGUMENTS
  integer, intent(in) :: n

  !     ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  !     LOCAL SCALARS
  integer :: i,j,k

  !     CLEAN-UP THE START STRUCTURE

  do k = 1,nbr
     i = br(1,k)
     j = br(2,k)
     start(i,j) = 0
  end do

  !     FILL-IN THE START STRUCTURE AGAIN
  nbr = 0
  do k = 1,nite
     call region(nregw,nregh,x(n)/2.0D0,x((k-1)*ndim+1),&
                 x((k-1)*ndim+2),i,j)
     if ( start(i,j) .eq. 0 ) then
        nbr = nbr + 1
        br(1,nbr) = i
        br(2,nbr) = j
     end if
     next(k) = start(i,j)
     start(i,j) = k
  end do

end subroutine classify

!     ******************************************************************
!     ******************************************************************

subroutine region(nregw,nregh,iterad,x,y,i,j)

  implicit none
  
  !     SCALAR ARGUMENTS
  integer, intent(in)      :: nregh,nregw
  integer, intent(out)     :: i,j
  real(kind=8), intent(in) :: iterad,x,y

  i = 1 + int( ( x + nregw * iterad ) / ( 2.0d0 * iterad ) )
  i = min( max( 1, i ), nregw )

  j = 1 + int( ( y + nregh * iterad ) / ( 2.0d0 * iterad ) )
  j = min( max( 1, j ), nregh )

end subroutine region

subroutine generate_x(x, n, W, H, ch, cw, A)
  implicit none
  integer :: i, j, caso  , lin(6), cox, coy, cb
  ! n is the number of variables
  ! A is the number of chairs
  integer, intent(in) :: n, A
  !x is the array with the points we will generate and give back
  real(kind=8), intent(inout) :: x(n)
  ! P and L are the dimensions of the room
  real(kind=8), intent(in):: W, H, ch, cw

  real(kind=8) :: Tz(A,2), Ts(A,2), Td(A,2), Pdist_t(A,6)
  real(kind=8) :: ajuste_t(A,6), provaL, P, L, limbase(A,6)
  real(kind=8) :: dist_t(1,6), dist, aux, Z(A,2), S(A,2), D(A,2)
  real(kind=8) :: camada, camadadupla, v(6), so_t(2)
  real(kind=8) :: pit2(4), pit1(4), coo(A,2), cod(A,2)

  !***adicionei***
  real(kind=8) :: raz
  integer :: m1, m2, alfa, omega,  linha1
  !***adicionei***

  !***adicionei***
  real(kind=8) :: d_v
  integer :: caso_v
  !***adicionei***

  P = W - cw / 2.0D0
  L = H - ch

  !***adicionei***
  raz = P/L
  m1 = NINT(SQRT(A*raz))
  m2 = NINT(SQRT(A/raz))
  alfa = MIN (m1, m2) - 2
  omega = MAX (m1, m2) + 2

  if (alfa .lt. 1) then
     alfa = 1
  end if

  if (omega .gt. A) then
     omega = A
  end if

  !=========  Zig-zag pattern   =========
  ! O_O_O
  ! _O_O_O
  ! O_O_O

  Z(1,1) = A
  Z(1,2) = 1

  do i=alfa,omega
     Z(i,1) = CEILING(A/real(i,8))
     Z(i,2) = i
  end do



  ! Zig-zag sizes (Tz):


  ! With more than 6 students, they will be distibuted in vertices
  !of equilateral triangles. The array "Tz" stores the bases
  !and the heights of them in function of the side "l".

  !The first column of Tz has the sum of the bases of the fist line.
  !The second one has the sum of the heights.

  do i=alfa,omega
     ! With n triangles, we have n - 1 bases and the sum of half a
     ! base of the last one of the second line.
     ! We also have n-1 heights
     Tz(i,1) = Z(i,1) - 0.5
     Tz(i,2) = (Z(i,2)-1.0)*SQRT(3.0)/2.0
  end do

  Tz(1,2) = 0.000001; !considering the width of a single line



  !read *, tt

  ! Bases in P
  do i=alfa,omega
     Pdist_t(i,1) = P / Tz(i, 1)
     provaL = Pdist_t(i,1) * Tz(i,2)
     if (provaL .gt. L) then
        ajuste_t(i,1) = provaL
        Pdist_t(i,1) = L / Tz(i,2)
     end if
     limbase(i,1)= Pdist_t(i,1)*Tz(i,1) !***adicionei***
  end do

  ! Heights in P
  do i=alfa,omega
     Pdist_t(i,2) = P / Tz(i,2)
     provaL = Pdist_t(i,2) * Tz(i,1)
     if (provaL .gt. L) then
        ajuste_t(i,2) = provaL
        Pdist_t(i,2) = L / Tz(i,1)
     end if
     limbase(i,2)= Pdist_t(i,2)*Tz(i,1) !***adicionei***
  end do



  !=========  Sandwich pattern   =========
  ! O_O_O
  ! _O_O_
  ! O_O_O



  S(1,1) = A
  S(1,2) = 1

  do i=alfa,omega
     S(i,1) = CEILING((A - CEILING(real(i)/2.0))/ real(i)) + 1.0
     S(i, 2) = i
  end do



  ! Sandwich sizes:
  ! First column has the sum of the bases of the first line
  ! Second column has the sum of the heights

  do i=alfa,omega
     ! With n triangles, we have n-1 bases and the sum of half a base
     !of the last one of the second line.
     ! We also have n-1 heights
     Ts(i,1) = S(i,1) - 1.0
     Ts(i,2) = (S(i,2) - 1.0)*SQRT(3.0)/2.0
  end do

  Ts(1,2) = 0.000001 ! Considering the width of a single row





  !Bases in P
  do i=alfa,omega
     Pdist_t(i,3) = P/Ts(i,1)
     provaL = Pdist_t(i,3) * Ts(i,2)
     if (provaL .gt. L) then
        ajuste_t(i,3) = provaL
        Pdist_t(i,3) = L / Ts(i,2)
     end if
     limbase(i,3)= Pdist_t(i,3)*Ts(i,1) !***adicionei***
  end do

  !Heights in P
  do i=alfa,omega
     Pdist_t(i,4) = P/Ts(i,2)
     provaL = Pdist_t(i,4)*Ts(i,1)
     if (provaL .gt. L) then
        ajuste_t(i,4) = provaL
        Pdist_t(i,4) = L/Ts(i,1)
     end if
     limbase(i,4)= Pdist_t(i,4)*Ts(i,1) !***adicionei***
  end do



  !=========  Double Layer pattern   =========
  ! O_O_O
  ! O_O_O
  ! _O_O_
  ! O_O_O

  camadadupla = 0

  D(1,1)=A
  D(1,2)=1

  do i=alfa,omega
     D(i,1)= CEILING(real((A - 1 + CEILING(real(i)/2.0))/real(i)))
     D(i,2)= i
  end do


  ! Double layer sizes:
  ! First row has the sum of the bases of the first line
  ! Second row has the sum of the heights

  do i=alfa,omega
     ! With n triangles, we have n-1 bases and the sum of half a base
     !of the last triangle of the second line.
     Td(i,1)= D(i,1)-1
     ! Here we have to consider the double layer before we think about the
     !triangles.
     Td(i,2)= (S(i,2) - 2.0)*SQRT(3.0)/2.0 + 1.0
  end do

  Td(1,2) = 0.000001 ! Considerando a altura de uma única fileira



  !Bases in P
  do i=alfa,omega
     Pdist_t(i,5)= P/Td(i,1)
     provaL = Pdist_t(i,5) * Td(i,2)
     if (provaL .gt. L) then
        ajuste_t(i,5) = provaL
        Pdist_t(i,5) = L / Td(i,2)
     end if
     limbase(i,5)= Pdist_t(i,5)*Td(i,1) !***adicionei***
  end do

  !Heights in P
  do i=alfa,omega
     Pdist_t(i,6)= P / Td(i,2)
     provaL = Pdist_t(i,6) * Td(i,1)
     if (provaL .gt. L) then
        ajuste_t(i,6) = provaL
        Pdist_t(i,6) = L / Td(i,1)
     end if
     limbase(i,6)= Pdist_t(i,6)*Td(i,1) !***adicionei***
  end do



  ! ==================================

  !indentifies the maximum distance and if the case demands to change P with L

  !dist_t = MAXVAL(Pdist_t)
  do i = 1,6
     aux = Pdist_t(alfa,i)
     lin(i) = 1
     do j = alfa,omega                                      !***Pode ser de alfa até omega?
        if(Pdist_t(j,i) > aux)then
           aux = Pdist_t(j,i)
           lin(i) = j                         !***adicionei
        end if
     end do
     dist_t(1,i) = aux
  end do

  dist = MAXVAL(dist_t)

  !caso = MAXLOC(dist_t)
  do i = 1,6
     if(dist_t(1,i) == dist) then
        caso = i
        EXIT
     end if
  end do



  !***adicionei***

  pit1(1) = ((dist_t(1,1)*P/limbase(lin(1),1))/2)**2
  pit1(2) = ((dist_t(1,2)*L/limbase(lin(2),2))/2)**2
  pit1(3) = ((dist_t(1,3)*P/limbase(lin(3),3))/2)**2
  pit1(4) = ((dist_t(1,4)*L/limbase(lin(4),4))/2)**2


  pit2(1) = dist_t(1,1)**2  * 0.75
  pit2(2) = dist_t(1,2)**2  * 0.75
  pit2(3) = dist_t(1,3)**2  * 0.75
  pit2(4) = dist_t(1,4)**2  * 0.75



  v(1)=sqrt(pit1(1)+pit2(1))
  v(2)=sqrt(pit1(2)+pit2(2))
  v(3)=sqrt(pit1(3)+pit2(3))
  v(4)=sqrt(pit1(4)+pit2(4))
  v(5)=dist_t(1,5)
  v(6)=dist_t(1,6)

  !caso_v = MAXLOC(v)
  !d_v = MAXVAL(v)
  d_v = v(1)
  caso_v = 1
  do i = 2,6
     if(v(i) > d_v)then
        d_v = v(i)
        caso_v = i
     end if
  end do

  if (d_v .gt. dist) then
     dist=dist_t(1,caso_v)
     caso=caso_v
  end if

  !***adicionei***

  !rounding the milimiters     !***Isso eu tirei pra comparar as medidas
  !dist = FLOOR(dist*1000.0)/1000.0

  !verifies if the bases were constructed over L or over P and change them if needed
  cb=0
  if (MOD(real(caso),2.0) == 0) then
     aux = P
     P = L
     L = aux
     cb=1
  end if

  !Verifies the double layer
  if (caso == 5 .OR. caso == 6) then
     camadadupla = 1
  end if



  !==========  Exporting the coordinates  ==========
  camada = 1.0

  coo(1,1) = 0.0
  coo(1,2) = 0.0

  do i=2,A
     coo(i,1)=coo(i-1,1)+dist;
     coo(i,2)=coo(i-1,2);
     if (coo(i,1) > (limbase(lin(caso), caso))*1.0001) then
        camada = camada + 1
        coo(i,2) = coo(i,2) + (dist * SQRT(3.0)/2.0)
        coo(i,1) = dist / 2.0 * (1.0 - MOD(camada,2.0))
        if (camadadupla == 1) then
           coo(i,2) = dist
           coo(i,1) = 0
           camadadupla = 0
           camada = 1
        end if
     end if
  end do

  cod=coo

  so_t=MAXVAL(cod, dim=1)


  if ((P-so_t(1)) > (L-so_t(2))) then
     if  (so_t(1) .ne. 0) then
        do i=1,A
           cod(i,1)= cod(i,1)*P/(so_t(1))
        end do
     end if
  else
     if  (so_t(2) .ne. 0) then
        do i=1,A
           cod(i,2)= cod(i,2)*L/(so_t(2))
        end do
     end if
  end if

  dist=v(caso)

  cox=1
  coy=2

  if (cb==1) then
     cox=2
     coy=1
  end if


  do i=1,n-1
     linha1=(i+1)/2
     if (MOD(real(i),2.0) == 1) then
        x(i)=cod(linha1,cox)
        !   xp(i)=cop(linha1,cox)
     else
        x(i)=cod(linha1,coy)
        !  xp(i)=cop(linha1,coy)
     end if
  end do

  x(n) = dist

  do i = 1, A

     x(2 * i - 1) = x(2 * i - 1) + cw / 2.0D0
     x(2 * i)     = x(2 * i)     + ch / 2.0D0
     
  end do

  return
end subroutine generate_x
