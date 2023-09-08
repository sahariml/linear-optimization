!Algorithme Simplexe
!M.L. Sahari (c)12/2003-01/2004
implicit none
real xr,att
real,dimension(:),allocatable::C,d,bbar,Mu,x,xB,xE,Cz,z,quotient
real, dimension(:,:),allocatable ::A,B,Binv
integer, dimension(:),allocatable ::Bi,Ei,r,k
integer m,n,i,j,ind,alg,nv
!---------------------------------Inititialisation--------------------------------------------
write(*,'(a34)',advance='no'),'Donner le nombre de variables nv='
read *,nv
write(*,'(a46)',advance='no'),'Donner le nombres d-equations (contraintes) m='
read *,m
n=m+nv
!----------------------------------------------------------------------------------------------
allocate(x(n),xB(m),xE(n-m),C(n),d(m),bbar(m),Mu(m),Bi(m),Ei(n-m),Cz(n-m),z(n-m),r(1),k(1),quotient(m))
allocate(A(m,n),B(m,m),Binv(m,m))
do i=1,nv
	write(*,'(a2)',advance='no'),'c('
	write(*,'(i1)',advance='no'),i
	write(*,'(a2)',advance='no'),')='
	read *,c(i)
enddo
do i=1,m
	write(*,'(a2)',advance='no'),'b('
	write(*,'(i1)',advance='no'),i
	write(*,'(a2)',advance='no'),')='
	read *,d(i)
enddo
do i=1,m
	do j=1,nv
		write(*,'(a3)',advance='no'),'a('
		write(*,'(i1)',advance='no'),i
		write(*,'(a1)',advance='no'),','
		write(*,'(i1)',advance='no'),j
		write(*,'(a2)',advance='no'),')='
		read *,A(i,j)
	enddo
enddo
do i=1,m
	a(i,nv+i)=1
enddo
!----------------------------------------------------------------------------------------------
!---------------------------Les indices de la prèmiere  Base-----------------------------------
do i=1,m
	Bi(i)=nv+i
enddo
!---------------------------Les indices hors base----------------------------------------------
do i=1,nv
	Ei(i)=i
10 enddo
!-------------------------------Affichage des données------------------------------------------
print *,'c=',c
print *,'A=',A
print *,'b=',d
!-------------------------------Affichage des résultats de la première iter--------------------
print *,'Iteration Num 01'
print *,'*********************'
print *,'Bi=',Bi
print *,'Ei=',Ei
print *,'XB=',d
print *,'XE=',xE
read *,att
!----------------------------------------------------------------------------------------------
algo: do alg=1,10
do i=1,m
	B(:,i)=A(:,Bi(i))
enddo
!------------------------------calcul de l'inverse de B---------------------------------------
call invmat(B,Binv,m)
!----------------------------------Calcul de B(bar)--------------------------------------------
Bbar=matmul(Binv,d)
!------------------------------------Calcul de Zj----------------------------------------------
do i=1,n-m
	Z(i)=dot_product(C(Bi),matmul(Binv,A(:,Ei(i))))
enddo
!-----------------------------l'indice de la vaiable entrante----------------------------------
do i=1,n-m
	CZ(i)=C(Ei(i))-Z(i)
enddo
r = MAXLOC ((CZ) , MASK = CZ > 0) ! une condition
if (r(1)==0) exit algo
!-------------------------------------Calcul de Mu-r-------------------------------------------
Mu=matmul(Binv,A(:,Ei(r(1))))
!------------------------------------Calcul du quotient----------------------------------------
do i=1,m
	if(Mu(i)>0) then
		quotient(i)=bbar(i)/Mu(i)
	else
		quotient(i)=-1
	end if
enddo
!---------------------------------l'indice de la vaiable sortante-------------------------------
k = MINLOC ((quotient) , MASK = quotient > 0) ! une condition
if (k(1)==0) then
	goto 20
endif
!--------------------------------------la nouvelle base----------------------------------------
ind=Bi(k(1))
Bi(k)=Ei(r)
Ei(r)=ind
!-------------------------------------Solution realisable--------------------------------------
xr=quotient(k(1)) !variable entrante
do i=1,m
	x(Bi(i))=bbar(i)-(Mu(i)*xr)
enddo
x(Bi(k))=xr
x(Ei(:))=0
!-------------------------------Affichage des résultats----------------------------------------
print *,'Iteration Num ',alg+1
print *,'*********************'
!----------------------------------------------------------------------------------------------
print *,'La matrice B'
do j=1,m
	print *,B(j,:)
enddo
print *,'L-invese de la matrice B'
do j=1,m
	print *,Binv(j,:)
enddo
!----------------------------------------------------------------------------------------------
print *,'Bi=',Bi
print *,'Ei=',Ei
print *, 'La solution de base:'
print *, 'XB=',x(Bi(:))
print *,'La valeur de la fonction objectif'
print *,'Z=',dot_product(C,x)
!----------------------------------------------------------------------------------------------
read *,att
enddo algo
print*, 'La solution X = ',x
print *,'est optimale !!!'
read *,att
goto 30
20 print*, 'La fonction objectif n-est pas bornée !!!'
read *,att
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