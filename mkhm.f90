      !
      !>import hamiltonian from file
      !
      subroutine import_hamiltonian(hop,rvec,nr,no,inflag,file_name,sb2c)
        implicit none
        integer(4),intent(in):: nr, no, inflag
        logical(4),intent(in):: sb2c
        real(8),intent(out):: rvec(nr,3)
        complex(8),intent(out):: hop(nr,no,no)
        character(len=*),intent(in):: file_name

        integer(4) i, j, k, l, m, sw
        real(8),allocatable::ndegen(:)
        real(8):: tmp(3), t1, t2

        sw=0
        select case(inflag)
        case(1)
           Open(100,file=trim(adjustl(file_name)),status='old')
           do i=1,no
              do j=1,no
                 do  k=1,nr
                    !read(100,'(2(E15.8,1X),2(E22.15,1X))')& !2D input
                    read(100,'(3(E10.3,1X),2(E22.15,1X))')& !3D input
                         rvec(k,1),rvec(k,2),rvec(k,3),hop(k,i,j)
                 end do
              end do
           end do
           close(100)
        case(2)
           Open(100,file=trim(adjustl(file_name))//'/ham_r.txt',status='old')
           Open(101,file=trim(adjustl(file_name))//'/irvec.txt',status='old')
           do i=1,nr
              read(101,*)rvec(i,1:3)
              do j=1,no
                 do  k=1,no
                    read(100,*)hop(i,j,k)
                 end do
              end do
           end do
           close(100)
           close(101)
        case(3)
           Open(100,file=trim(adjustl(file_name))//'_hr.dat',status='old')
           read(100,*)
           read(100,*)l !no
           read(100,*)m !nr
           if(.not. no==l)then
              print *,'orbital number conflict between input and hamiltonian'
              stop
           end if
           if(.not. nr==m)then
              print *,'number of hoppings conflict between input and hamiltonian'
              stop
           end if
           allocate(ndegen(nr))
           read(100,'(15F5.0)') (ndegen(i),i=1,nr)
           do  i=1,nr
              do j=1,no
                 do k=1,no
                    read(100,'(5F5.0,2F12.6)')rvec(i,:),l,m,hop(i,k,j)
                    hop(i,k,j)=hop(i,k,j)/ndegen(i)
                 end do
              end do
           end do
           deallocate(ndegen)
           close(100)
        case(4)
           sw=1
           Open(100,file=trim(adjustl(file_name))//'/Hopping.dat',status='old')
           read(100,*)
           do i=1,3
              read(100,*)
           end do
           read(100,*) l,m,i
           if(.not. no==l)then
              print *,'orbital number conflict between input and hamiltonian'
              stop
           end if
           if(.not. nr==m)then
              print *,'number of hoppings conflict between input and hamiltonian'
              stop
           end if
           read(100,*)
           read(100,*)
           do i=1,no
              read(100,*)
           end do
           do  k=1,nr
              do i=1,no
                 do j=1,no
                    read(100,*)rvec(k,:),tmp(:),l,m,hop(k,i,j)
                 end do
              end do
           end do
           close(100)
        case(0) !square lattice
           t1=-1.0d0

           hop=0.0d0
           hop(1,1,1)=t1
           hop(2,1,1)=t1
           hop(3,1,1)=t1
           hop(4,1,1)=t1

           rvec=0.0d0
           rvec(1,1)=1.0d0
           rvec(2,1)=-1.0d0
           rvec(3,2)=1.0d0
           rvec(4,2)=-1.0d0
        case(5) !honeycomb lattice
           t1=-1.0d0

           hop=0.0d0
           hop(1,2,1)=t1  !( 1, 1)
           hop(2,2,1)=t1  !( 0, 1)
           hop(3,1,2)=t1  !(-1,-1)
           hop(4,1,2)=t1  !( 0,-1)
           hop(5,1,2)=t1  !( 0, 0)
           hop(5,2,1)=t1

           rvec=0.0d0
           rvec(1,1)=1.0d0  !( 1, 1)
           rvec(1,2)=1.0d0  !( 1, 1)
           rvec(2,2)=1.0d0  !( 0, 1)
           rvec(3,1)=-1.0d0 !(-1,-1)
           rvec(3,2)=-1.0d0 !(-1,-1)
           rvec(4,2)=-1.0d0 !( 0,-1)
        case(6) !triangular lattice
           t1=-1.0d0

           hop=0.0d0
           hop(1,1,1)=t1
           hop(2,1,1)=t1
           hop(3,1,1)=t1
           hop(4,1,1)=t1
           hop(5,1,1)=t1
           hop(6,1,1)=t1

           rvec=0.0d0
           rvec(1,1)=1.0d0
           rvec(2,1)=-1.0d0
           rvec(3,2)=1.0d0
           rvec(4,2)=-1.0d0
           rvec(5,1)=1.0d0
           rvec(5,2)=1.0d0
           rvec(6,1)=-1.0d0
           rvec(6,2)=-1.0d0
        case(7) !E_g square lattice
           t1=-1.0d0
           t2=t1*0.1d0

           hop=0.0d0
           !nn
           hop(1,1,1)=t1
           hop(2,1,1)=t1
           hop(3,2,2)=t1
           hop(4,2,2)=t1

           hop(1,2,2)=t2
           hop(2,2,2)=t2
           hop(3,1,1)=t2
           hop(4,1,1)=t2
           !nnn
           hop(5,1,2)=t2
           hop(5,2,1)=t2
           hop(6,1,2)=-t2
           hop(6,2,1)=-t2
           hop(7,1,2)=t2
           hop(7,2,1)=t2
           hop(8,1,2)=-t2
           hop(8,2,1)=-t2

           rvec=0.0d0
           !nn
           rvec(1,1)=1.0d0
           rvec(1,2)=0.0d0
           rvec(2,1)=-1.0d0
           rvec(2,2)=0.0d0
           rvec(3,1)=0.0d0
           rvec(3,2)=1.0d0
           rvec(4,1)=0.0d0
           rvec(4,2)=-1.0d0
           !nnn
           rvec(5,1)=1.0d0
           rvec(5,2)=1.0d0
           rvec(6,1)=-1.0d0
           rvec(6,2)=1.0d0
           rvec(7,1)=-1.0d0
           rvec(7,2)=-1.0d0
           rvec(8,1)=1.0d0
           rvec(8,2)=-1.0d0
        case default
           write(*,*) "input flag is wrong"
           stop
        end select

        if(inflag/=0)then
           do  k=1,no
              do j=1,no
                 do i=1,nr
                    l=int(dble(hop(i,j,k))*1.0d6)
                    m=int(aimag(hop(i,j,k))*1.0d6)
                    hop(i,j,k)=cmplx(dble(l)*1.0d-6,dble(m)*1.0d-6)
                 end do
              end do
           end do
        end if

        if(sb2c)call bcc2desc(rvec,nr,sw)

        return
      contains
        !
        !> this subroutine only 122 iron-pnictide 
        !! change axis bcc to sc
        !
        subroutine bcc2desc(r,M,sw)
          implicit none
          integer(4),intent(in):: M,sw
          real(8),intent(inout):: r(M,3)

          integer(4) i
          real(8),dimension(3):: re

          do i=1,M
             select case(sw)
             case(0) !PWSCF
                re(1)=(r(i,1)+r(i,2)-r(i,3))
                re(2)=(-r(i,1)+r(i,2)-r(i,3))
                re(3)=(r(i,1)+r(i,2)+r(i,3))
             case(1) !ecalj
                re(1)=2.0d0*r(i,1)-r(i,3)
                re(2)=2.0d0*r(i,2)-r(i,3)
                re(3)=r(i,3)
             case default !wien BCC
                re(1)=(-r(i,1)+r(i,2)+r(i,3))
                re(2)=(r(i,1)-r(i,2)+r(i,3))
                re(3)=(r(i,1)+r(i,2)-r(i,3))
             end select
             r(i,1)=(re(1)-re(2))*0.5d0
             r(i,2)=(re(1)+re(2))*0.5d0
             !r(i,1)=re(1)
             !r(i,2)=re(2)
             r(i,3)=re(3)
          end do

          return
        end subroutine bcc2desc
      end subroutine import_hamiltonian
      !
      !> import hamiltonian and calculate eigenvalue
      !
      subroutine mkhm(nqx,nqy,nqz,no,nr,eq,en,uni,sb2c,hamname,inflag,mass)
        implicit none
        integer(4),intent(in):: nr,no,nqx,nqy,nqz,inflag
        logical(4),intent(in):: sb2c
        real(8),intent(in):: mass
        real(8),intent(out):: eq(0:nqx-1,0:nqy-1,0:nqz-1,no)
        complex(8),intent(out):: uni(0:nqx-1,0:nqy-1,0:nqz-1,no,no)
        complex(8),intent(inout):: en(0:nqx-1,0:nqy-1,0:nqz-1,no,no)
        character(len=*),intent(in):: hamname

        real(8),parameter:: pi=4.0d0*atan(1.0d0)
        integer(4) i,j,k,l,m,n
        real(8) eq1(no),rvec(nr,3),phase
        complex(8) hop(nr,no,no),en2(no,no)

        ! Lapack blas
        integer(4) info
        real(8) rwork(3*no-2),tmp
        complex(8) work(2*no-1)

        call import_hamiltonian(hop,rvec,nr,no,inflag,hamname,sb2c)

        !$omp parallel
        !$omp workshare
        en=0.0d0
        !$omp end workshare
        do l=1,no
           do m=l,no
              do n=0,nqz-1
                 !$omp do private(i,k,phase)
                 do j=0,nqy-1
                    do i=0,nqx-1
                       do k=1,nr
                          !if (rvec(k,3)==0)then
                             phase=2*i*pi/nqx*rvec(k,1)+2*j*pi/nqy*rvec(k,2)&
                                  +2*n*pi/nqz*rvec(k,3) !k@
                             en(i,j,n,l,m)=en(i,j,n,l,m)&
                                  +hop(k,l,m)*cmplx(cos(phase),-sin(phase))
                          !end if
                       end do
                       en(i,j,n,m,l)=conjg(en(i,j,n,l,m))
                    end do
                 end do
                 !$omp end do
              end do
           end do
        end do

        !$omp workshare
        en=en/mass
        !$omp end workshare

        do n=0,nqz-1
           !$omp do private(i,l,m,en2,eq1,work,rwork,info)
           do j=0,nqy-1
              do i=0,nqx-1
                 en2=en(i,j,n,:,:)
! Lapack
                 call zheev('V','U',no,en2,no,eq1,work,2*no-1,rwork,info)
                 uni(i,j,n,:,:)=en2
                 eq(i,j,n,:)=eq1
              end do
           end do
           !$omp end do
        end do
        !$omp end parallel

        return
      end subroutine mkhm
