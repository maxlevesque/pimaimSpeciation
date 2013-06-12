program polym

! analyze polymer species in mixtures

implicit none

integer :: nmax
parameter (nmax=2000)

integer :: i,j,k,l
integer :: imark,jmark
integer :: nconfigs,nequil,nanion,ncation1,ncation2,nion
integer :: nneighb,nshar,nchain,nat,ichain,jchain,ntri,nmol

integer, dimension(10)        :: shar
integer, dimension(nmax)      :: chain,nF,nBe,ncharge
integer, dimension(nmax,nmax) :: atoms

double precision :: boxlength,halfbox,halfboxrec
double precision :: du,dx,dy,dz,r2,totfrac,nchargetot
double precision :: fracijBe,fracijF

double precision, dimension(3) :: rdfmin, rdfmin2

double precision, dimension(nmax)   :: x,y,z
double precision, dimension(0:nmax,0:nmax) :: spec 

character*80 :: filein1

!Reading the input file
open(10,file='speciation.inpt')

read(10,*) nconfigs
read(10,*) nanion
read(10,*) ncation1
read(10,*) ncation2
read(10,*) boxlength

nion=nanion+ncation1+ncation2

do i=1,3
   read(10,*) rdfmin(i)
   rdfmin2(i)=rdfmin(i)**2
enddo

read(10,'(a)') filein1

open(11,file=filein1,status='old')
close(10)

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

nchargetot=0
nmol=0
totfrac=0

do i=0,nmax
   do j=0,nmax
      spec(i,j)=0
   enddo
enddo

!loop over nconfigs
do l=1,nconfigs

!Reading positions.out
do i=1,nion
   read(11,*) x(i),y(i),z(i)
enddo

!Initialize variables
do i=1,nmax
   chain(i)=0
   do j=1,nmax
      atoms(i,j)=0
   enddo
enddo

nchain=0
ntri=0

!loop 1 over cations

do i=1,ncation2

   imark=i+nanion+ncation1
   nneighb=0
   

   !loop 2 over cations
   do j=1,i-1

      jmark=j+nanion+ncation1

      dx=x(imark)-x(jmark)
      dy=y(imark)-y(jmark)
      dz=z(imark)-z(jmark)

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      !tests  if they are bonded
      !2 tests: d(Be-Be)<rdfmin and presence of 1 or more bridging F

       if (r2.le.rdfmin2(3)) then

          call linkage(imark,jmark,nshar,nanion,rdfmin2(2),boxlength,x,y,z,shar)

          if (nshar.ne.0) then 
             nneighb=nneighb+1
             ichain=chain(imark)
             jchain=chain(jmark)


   !1st case: i is not in a chain then it goes in the same as j
             if (ichain.eq.0) then
                chain(imark)=jchain
   !1st element of atoms contains the number of atoms in the chain
   ! then the others are the atom labels of atoms contained 
   !in the chain
                atoms(jchain,1)=atoms(jchain,1)+1
                atoms(jchain,atoms(jchain,1)+1)=imark


   !2nd case: i and j are in the same chain
             elseif (ichain.eq.jchain) then
                ntri=ntri+1


    ! 3rd case: i is in a chain of greater label than j
    ! then all the atoms in the same chain as i go into
    ! the chain containing j
    ! and the former chain containing i is filled with 0
             elseif (ichain.gt.jchain) then
                nat=atoms(ichain,1)
    ! must separate i and the other atoms of his old chain
                do k=1,nat
                   if (atoms(ichain,1+k).ne.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(ichain,1+k).eq.imark) then
                      chain(atoms(ichain,1+k))=jchain
                      atoms(jchain,1)=atoms(jchain,1)+1
                      atoms(jchain,atoms(jchain,1)+1)=atoms(ichain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(ichain,k)=0
                enddo



   ! 4th case: i is in a chain of smaller label than j
   !then all the atoms in the same chain as j go into
   !the chain containing i
   !and the former chain containing j is filled with 0
             else
                nat=atoms(jchain,1)
   ! must separate j and the other atoms of his old chain
                do k=1,nat
                   if (atoms(jchain,1+k).ne.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat
                   if (atoms(jchain,1+k).eq.jmark) then
                      chain(atoms(jchain,1+k))=ichain
                      atoms(ichain,1)=atoms(ichain,1)+1
                      atoms(ichain,atoms(ichain,1)+1)=atoms(jchain,1+k)
                   endif
                enddo
                do k=1,nat+1
                   atoms(jchain,k)=0
                enddo
             endif

   ! Add the shared F (given by subroutine linkage) in the good chain 
             do k=1,nshar
                chain(shar(k))=chain(imark)
                atoms(chain(imark),1)=atoms(chain(imark),1)+1
               atoms(chain(imark),atoms(chain(imark),1)+1)=shar(k)
             enddo
          endif
       endif
    enddo
  
  ! if i has no neighbour then creation of a new chain
    if (nneighb.eq.0) then
       nchain=nchain+1
       chain(imark)=nchain
       atoms(nchain,1)=1
       atoms(nchain,2)=imark
    endif
 enddo


 ! Add the unshared F in the good chains
 do i=1,nanion

    j=0

    do while (chain(i).eq.0 .and. j.lt.ncation2)
       
      j=j+1
      jmark=nanion+ncation1+j
      dx=x(jmark)-x(i) 
      dy=y(jmark)-y(i) 
      dz=z(jmark)-z(i) 

      dx=dx-boxlength*int(dx*halfboxrec)
      dy=dy-boxlength*int(dy*halfboxrec)
      dz=dz-boxlength*int(dz*halfboxrec)

      r2=dx**2+dy**2+dz**2

      if (r2.le.rdfmin2(2)) then
         chain(i)=chain(jmark)
         atoms(chain(i),1)=atoms(chain(i),1)+1
         atoms(chain(i),atoms(chain(i),1)+1)=i
      endif

    enddo
  ! if F is alone then creates a new chain
     if (chain(i).eq.0) then
        nchain=nchain+1
        chain(i)=nchain
        atoms(nchain,1)=1
        atoms(nchain,2)=i
     endif
 enddo






 do i=1,nmax
    nF(i)=0
    nBe(i)=0
    ncharge(i)=0
 enddo

 do i=1,nchain
    
       do j=1,atoms(i,1)
          if (atoms(i,j+1).le.nanion) then
             nF(i)=nF(i)+1
             ncharge(i)=ncharge(i)-1
          else
             nBe(i)=nBe(i)+1
             ncharge(i)=ncharge(i)+2
          endif
       enddo
       nchargetot=nchargetot+ncharge(i)
       
       if (nF(i).ne.0) then
          spec(nBe(i),nF(i))=spec(nBe(i),nF(i))+1
          nmol=nmol+1
       endif
       

 enddo




 enddo

 open(21,file='speciationBe.dat')
 open(22,file='speciationF.dat')

 do i=0,nmax
    do j=0,nmax
       fracijBe=100.0*(spec(i,j))*float(i)/(float(nconfigs)*float(ncation2))
       totfrac=totfrac+fracijBe
       fracijF=100.0*(spec(i,j))*float(j)/(float(nconfigs)*float(nanion))
       if(fracijBe.gt. 0.1) then
         write(21,*) fracijBe,'% de Be ',i,'F ',j
       endif
       if(fracijF.gt. 0.1) then
         write(22,*) fracijF,'% de Be ',i,'F ',j
       endif
    enddo
 enddo
 nchargetot=nchargetot/float(nconfigs)
!~  write(6,*) nchargetot,totfrac ! unverbose because usually lots and lots of serial executions are done and this completely mess the standard output.

 close(21)
 close(22)
 close(11)


 end



subroutine linkage(imark,jmark,nshar,nanion,rdfcut2,boxlength,x,y,z,shar)

implicit none

integer :: nmax

parameter (nmax=2000)

integer :: i,j,k,nshar,nanion
integer :: imark,jmark
integer, dimension(10) :: shar

double precision, dimension(nmax) :: x,y,z
double precision dx1,dx2,dy1,dy2,dz1,dz2,dr1,dr2,rdfcut2
double precision boxlength,halfbox,halfboxrec

nshar=0

halfbox=boxlength/2.0d0
halfboxrec=1.0d0/halfbox

do i=1,nanion

   dx1=x(i)-x(imark)
   dy1=y(i)-y(imark)
   dz1=z(i)-z(imark)

   dx1=dx1-boxlength*int(dx1*halfboxrec)
   dy1=dy1-boxlength*int(dy1*halfboxrec)
   dz1=dz1-boxlength*int(dz1*halfboxrec)

   dx2=x(i)-x(jmark)
   dy2=y(i)-y(jmark)
   dz2=z(i)-z(jmark)

   dx2=dx2-boxlength*int(dx2*halfboxrec)
   dy2=dy2-boxlength*int(dy2*halfboxrec)
   dz2=dz2-boxlength*int(dz2*halfboxrec)

   dr1=dx1**2+dy1**2+dz1**2
   dr2=dx2**2+dy2**2+dz2**2

   if (dr1.le.rdfcut2 .and. dr2.le.rdfcut2) then
      nshar=nshar+1
      shar(nshar)=i
   endif
enddo

return

end
