program testes
  use packmod
  implicit none

  character(len=15) :: LTEXSOL = 'solution.tex', &
                       JSONSOL = 'solution.json'
  integer           :: MAXMEM = 1, nmem = 1, n, status

	real(kind=8) :: W, H, t, qnt
	real(kind=8), pointer :: x(:), xp(:), xb(:,:)
	
  nregw = CEILING(W / MINDIST)
  nregh = CEILING(H / MINDIST)

  allocate(start(0:nregw+1,0:nregh+1), stat=status)
  if ( status .ne. 0 ) then
    write(*,*) 'Allocation error in main program'
    stop
  end if
  
  write(*,*) 'Enter the minimum distance between students:'
  read(*,*) MINDIST
  
  write(*,*) 'Enter the dimensions of the room (W H): '
  read(*,*) W,H
  
  write(*,*) 'Enter the number of chairs: '
  read(*,*) nite
  
  n=1+nite*2

	!qnt = 0.1        !Coeficiente de perturbacao
	
	allocate(x(n),xp(n),xb(n,MAXMEM))

	call generate_x(x, n, W, H, nite)

	write(*,*) 'Solution'
	write(*,*)  x
	
	
	xb(:,1) = x(:)
		
	call drawsol(nite,W,H,n,nmem,xb(1:n,1:nmem),LTEXSOL)
	call tojson(n,nmem,xb(1:n,1:nmem),nite,W,H,JSONSOL,.true.)


end program testes


subroutine generate_x(x, n, W, H, A)
	implicit none
	integer :: i, j, caso  , lin(6), cox, coy, cb
	! n is the number of variables
	! A is the number of chairs
	integer, intent(in) :: n, A
	!x is the array with the points we will generate and give back
	real(kind=8), intent(inout) :: x(n)
	! P and L are the dimensions of the room
	real(kind=8), intent(in):: W, H

	real(kind=8) :: Tz(A,2), Ts(A,2), Td(A,2), Pdist_t(A,6)
        real(kind=8) :: ajuste_t(A,6), provaL, P, L, limbase(A,6)
	real(kind=8) :: dist_t(1,6), dist, aux, Z(A,2), S(A,2), D(A,2)
        real(kind=8) :: camada, camadadupla, v(6), so_t(2)
        real(kind=8) :: pit2(4), pit1(4), coo(A,2), cod(A,2)

	!***adicionei***
	real(kind=8) :: raz,t
	integer :: m1, m2, alfa, omega,  linha1, ultimalinha
	!***adicionei***

	!***adicionei***
	real(kind=8) :: q, d_v
	integer :: caso_v
	!***adicionei***

	P = W
	L = H

	PRINT *, 'Entrou na subrotina'

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
		aux = Pdist_t(1,i)
		lin(i) = 1
		do j = 1,A                                      !***Pode ser de alfa até omega?
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
      
       !Falta converter nesse formato
       x(n) = dist
      ! xp(n)= dist


      return
      end subroutine generate_x



subroutine perturbation_x(x, xp, n, A, qnt)
	implicit none
	integer, intent(in) :: n, A
	!x is the array with the points we will generate and give back
	real(kind=8), intent(inout) :: x(n), xp(n)
	real(kind=8), intent(in) :: qnt
	real(kind=8) :: cop(A,2)
	integer :: i, linha1, hv, ultimalinha


        
	if (x(2) == x(4)) then
	 hv=2
	else
	 hv=1
	end if


	do i=1,n-1,2

	   cop((i+1)/2,1)=x(i)
	   cop((i+1)/2,2)=x(i+1)

	end do



	i=1;
	do while (cop(1,hv)==cop(i,hv))

	cop(i,hv)=(1-mod(real(i),2.0))*qnt*x(n)/2
	i=i+1
	if (i .gt. A) then
	    exit
	end if
	end do

	linha1=i-1

	    i=1

	do while (cop(A,hv)==cop(A-i+1,hv))
	cop(A-i+1,hv)=cop(A,hv)-(1-mod(i,2))*0.1*x(n)/2
	i=i+1
	if (i .gt. A) then

	    exit
	end if
	end do


	ultimalinha=A-i

	do i=linha1,ultimalinha

	cop(i+1,hv)=cop(i+1,hv)+(0.1*x(n))*((-1)**(i+1))

	end do

	do i=1,n-1
	 linha1=(i+1)/2
	 if (MOD(real(i),2.0) == 1) then
	     xp(i)=cop(linha1,1)
	 else
	     xp(i)=cop(linha1,2)
	 end if
	end do

	xp(n)= x(n)


	return
end subroutine perturbation_x
       
       
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
	real(kind=8) :: scale, minover, radius
	
	open(unit=10,file=solfile)
	write(10,01)
	do j = 1,nmem
		if ( j .gt. 1 ) write(10,12)
		
		!radius = minover(n,x(1:n,j)) / 2.0D0
		radius = x(n,j) / 2.0D0

		! SCALING
		scale = min(20.0D0 / (H + 2.0D0) * radius, &
					 10.0D0 / (W + 2.0D0 * radius))
		write(*,*) 'haha: ',2.0D0*radius,nite,scale,scale
		write(10,10) 2.0D0*radius,nite,scale,scale

		! CLASSROOM
		write(10,25) W,H
		! PROFESSOR
		write(10,26) W+1.0D0,W+5.0,H


		! CIRCULAR ITEMS
		do i = 1,nite
			write(10,20) x(2*i-1,j),x(2*i,j),i,x(2*i-1,j),x(2*i,j), radius
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
