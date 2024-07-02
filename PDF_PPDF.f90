!  PDF_PPDF.f90 
!****************************************************************************
!
!****************************************************************************

program PDF_PPDF
implicit none

! Variables
INTEGER*4::i,Natoms,Natoms1,K,l,m,lay,number_of_histogram_bins
character(len=40)::Str(5) 
character(len=50)::filename,out_filename
REAL*8::VL(2,1:3)
REAL*8::VX,VY,VZ
REAL*8,allocatable::gr(:),gr11(:),gr22(:),gr12(:)
REAL*8::V,rmax,drr,dr, particle_density,rho,rh11,rh22,rh12,rh21,v_1,r,r1,r2
integer,allocatable::ID(:),iType(:)
REAL*8,allocatable::SR4(:),SR5(:),SR6(:)  

open(1,file='PDF_PPDF.INI',action='read')
write(*,*)"**************************************************"
read(1,*)
read(1,*)filename
write(*,*)'filename:',filename
read(1,*)
read(1,*)rmax
write(*,*)'cutoff radius:',rmax
read(1,*)
read(1,*)number_of_histogram_bins
write(*,*)'number of histogram bins:',number_of_histogram_bins
read(1,*)
read(1,*)out_filename
write(*,*)'out_filename:',out_filename
write(*,*)"**************************************************"
close(1)


! Body of PDF_PPDF
open(1,file=Trim(adjustl(filename)),status='old',action='read')
do i=1, 3 
    read(1,*)
end do ! 空过前三行 空行的时候直接进行读入，不用填写变量

read(1,*)Natoms
read(1,*)str(4) ! 继续空一行fil

do k=1,3
    read(1,*)VL(1,k),VL(2,k)
end do

VX=VL(2,1)-VL(1,1)
VY=VL(2,2)-VL(1,2)
VZ=VL(2,3)-VL(1,3)
V = VX*VY*VZ

read(1,*)str(5) !继续空一行
!=====read Particle coordinates===========
allocate(ID(Natoms),iType(Natoms))          ! NTotal is total atoms 给定义的数组分配大小
allocate(SR4(Natoms),SR5(Natoms),SR6(Natoms))
Natoms1 = 0
do k=1, Natoms
    read(1,*)ID(k),iType(k),SR4(k),SR5(k),SR6(k)
    if(iType(k)==1)then
        Natoms1 = Natoms1 + 1
    end if
end do 
close(1)


!======compute g(r)==================
allocate(gr(number_of_histogram_bins),gr11(number_of_histogram_bins))
allocate(gr12(number_of_histogram_bins),gr22(number_of_histogram_bins))

gr = 0
gr11 = 0
gr12 = 0
gr22 = 0

dr = rmax/number_of_histogram_bins
drr = 0
do l=1, Natoms
    do m=1, Natoms
        if(l==m)then
            cycle
        end if   
        call drx(SR4(l),SR5(l),SR6(l),SR4(m),SR5(m),SR6(m),drr,VX,VY,VZ)
        if(drr<=rmax)then
            lay = drr / dr + 1  ! 壳层数
            gr(lay) = gr(lay)+1  ! 壳层内的原子数
            
            if(iType(l)==1.and.iType(m)==1)then
                gr11(lay) = gr11(lay)+1
            else if(iType(l)==2.and.iType(m)==2)then
                gr22(lay) = gr22(lay)+1
            else if(iType(l)==1.and.iType(m)==2)then
                gr12(lay) = gr12(lay)+1
!             else if(iType(l)==2.and.iType(m)==1)then
!                 gr21(lay) = gr21(lay)+1
            end if    
        end if      
    end do
end do    

r1 = 0
r2 = dr
rho = 0
write(*,*)'-------------------------------------------'
do i=1,100
    v_1 = (4*3.14159265358979323846*r2*r2*r2-4*3.14159265358979323846*r1*r1*r1)/3   ! 壳层体积
    rho = (gr(i)*V)/(v_1*Natoms*Natoms)
    rh11 = (gr11(i)*V)/(v_1*Natoms1*Natoms1)
    rh22 = (gr22(i)*V)/(v_1*(Natoms-Natoms1)*(Natoms-Natoms1))
    rh12 = (gr12(i)*V)/(v_1*Natoms1*(Natoms-Natoms1))
!    rh21 = (gr21(i)*V)/(v_1*Natoms1*(Natoms-Natoms1))
    open(unit=101,file=Trim(adjustl(out_filename)),action = 'write')
    if(i==1)then
        write(101,*)'r tot 1-1 1-2 2-2'    
    end if
    write(101,'(F5.3,1X,5(F10.6))')r2-(dr/2),rho,rh11,rh12,rh22
    r1 = r1 + dr
    r2 = r2 + dr
end do

close(101)
deallocate(ID,iType) 
deallocate(SR4,SR5,SR6)
deallocate(gr,gr11,gr12,gr22)
write(*,*)"complete! Bai bai la, friend!"
PAUSE 
    
end program PDF_PPDF

subroutine drx(x1,y1,z1,x2,y2,z2,drr,VX1,VY1,VZ1)        ! 动态数组不能直接作为子程序的形参传递
    implicit none
    real*8::x1,y1,z1,x2,y2,z2,drr,edge,VX1,VY1,VZ1
  
   
    if (abs(x2-x1)>VX1/2)then
        if(x2>=x1)then
            x2 = x2 - VX1
        else
            x2 = x2 + VX1
        end if
    end if 

      
    if (abs(y2-y1)>VY1/2)then
        if(y2>=y1)then
            y2 = y2 - VY1
        else
            y2 = y2 + VY1
        end if
    end if  
   
    if (abs(z2-z1)>VZ1/2)then
        if(z2>=z1)then
            z2 = z2 - VZ1
        else
            z2 = z2 + VZ1
        end if
    end if 

    drr = sqrt((abs(x2-x1))**2+abs((y2-y1))**2+(abs(z2-z1))**2)
    return 
    
end subroutine drx