!Algorithme Simplexe
!M.L. Sahari (c)12/2003-01/2004
implicit none
real xr,att
real,dimension(:),allocatable::C,d,bbar,Mu,x,xB,xE,Cz,z,quotient
real, dimension(:,:),allocatable ::A,B,Binv
integer, dimension(:),allocatable ::Bi,Ei,r,k
integer m,n,i,j,ind,alg
!---------------------------------Inititialisation--------------------------------------------
write(*,'(a32)',advance='no'),'Donner le nombre de variables n='
read *,n
write(*,'(a46)',advance='no'),'Donner le nombres d-equations (contraintes) m='
read *,m
!----------------------------------------------------------------------------------------------
allocate(x(n),xB(m),xE(n-m),C(n),d(m),bbar(m),Mu(m),Bi(m),Ei(n-m),Cz(n-m),z(n-m),r(1),k(1),quotient(m))
allocate(A(m,n),B(m,m),Binv(m,m))
do i=1,m
	write(*,'(a3)',advance='no'),'c('
	write(*,'(i1)',advance='no'),i
	write(*,'(a2)',advance='no'),')='
	read *,c(i)
enddo
do i=1,m
	write(*,'(a3)',advance='no'),'b('
	write(*,'(i1)',advance='no'),i
	write(*,'(a2)',advance='no'),')='
	read *,d(i)
enddo
do i=1,m
	do j=1,n
		write(*,'(a3)',advance='no'),'a('
		write(*,'(i1)',advance='no'),i
		write(*,'(a1)',advance='no'),','
		write(*,'(i1)',advance='no'),j
		write(*,'(a2)',advance='no'),')='
		read *,A(i,j)
	enddo
enddo
!----------------------------------------------------------------------------------------------
print *,'Donner les indices de la prèmiere  Base'
do i=1,m
	write(*,'(a1)',advance='no'),'B'
	write(*,'(i1)',advance='no'),i
	write(*,'(a1)',advance='no'),'='
	read *,Bi(i)
enddo
!trouver les indices hors base
ind=1
do i=1,n
	do j=1,m
		if (i==Bi(j)) then
			goto 10
		end if
	enddo
	ei(ind)=i
	ind=ind+1
10 enddo
!-------------------------------Affichage des données------------------------------------------
print *,'c=',c
print *,'A=',a
print *,'b=',d
!-------------------------------Affichage des résultats de la première iter--------------------
print *,'Iteration Num 01'
print *,'Bi=',Bi
print *,'Ei=',ei
print *,'XB=',d
print *,'XE=',xE
!----------------------------------------------------------------------------------------------
algo: do alg=1,10
do i=1,m
	B(:,i)=a(:,Bi(i))
enddo
!calcul de l'inverse de B
call invmat(B,Binv,m)
print *,'La matrice B'
do j=1,m
	print *,B(j,:)
enddo
print *,'L-invese de la matrice B'
do j=1,m
	print *,Binv(j,:)
enddo
!---------------------------------------------------------------------------------------------
!Calcul de B(bar)-----------------------------------------------------------------------------
Bbar=matmul(Binv,d)
!Calcul de Zj---------------------------------------------------------------------------------
do i=1,n-m
	Z(i)=dot_product(C(Bi),matmul(Binv,A(:,Ei(i))))
enddo
!l'indice de la vaiable entrante--------------------------------------------------------------
do i=1,n-m
	CZ(i)=C(Ei(i))-Z(i)
enddo
r = MAXLOC ((CZ) , MASK = CZ > 0) ! une condition
if (r(1)==0) exit algo
!Calcul de Mu-r--------------------------------------------------------------------------------
Mu=matmul(Binv,A(:,Ei(r(1))))
!Calcul du quotient----------------------------------------------------------------------------
do i=1,m
	if(Mu(i)>0) then
		quotient(i)=bbar(i)/Mu(i)
	else
		quotient(i)=-1
	end if
enddo
!l'indice de la vaiable sortante---------------------------------------------------------------
k = MINLOC ((quotient) , MASK = quotient > 0) ! une condition
if (k(1)==0) then
	goto 20
endif
!la nouvelle base------------------------------------------------------------------------------
ind=Bi(k(1))
Bi(k)=Ei(r)
Ei(r)=ind
!Solution realisable
xr=quotient(k(1)) !variable entrante
do i=1,m
	x(Bi(i))=bbar(i)-(Mu(i)*xr)
enddo
x(Bi(k))=xr
x(Ei(:))=0
!-------------------------------Affichage des résultats----------------------------------------
print *,'Iteration N°',alg+1
print *,'Bi=',Bi
print *,'Ei=',ei
print *,'XB=',x(Bi(:))
print *,'Cr-zr=',CZ(r(1))
print *,'bbar_k/mu_kr=',quotient(k(1))
print *,'XE=',xE
print *,'Z(x)=',dot_product(C,x)
read *,att
enddo algo
print*, 'La solution:',x
print *,'est optimale !!!'
goto 30
20 print*, 'La fonction objectif n-est pas bornée !!!'
contains
!----------------------------------------------------------------------------------------------
!les sous-programme----------------------------------------------------------------------------
subroutine invmat(B,Binv,m)
use imsl 
integer m
real, dimension(m,m)::B,Binv,AINV(m,m)
CALL LINRG (m, B, m, Binv, m) 
end subroutine invmat
!-----------------------------------------------------------------------------------------------
30 end