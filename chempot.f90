      !
      !> this subroutine calcute chemical potential mu as
      !!
      !! sum_k gr(k,mu)=rfill 
      !!
      !! we calculate mu using Brent's method or bisection method
      !
      subroutine mkrmu(nqx,nqy,nqz,no,temp,rfill,eq,rmu,rmuOLD,esterr,RPA)
        implicit None
        integer(4),intent(in):: nqx,nqy,nqz,no
        logical(4),intent(in):: RPA
        real(8),intent(in):: temp, rfill
        real(8),intent(in),dimension(0:nqx-1,0:nqy-1,0:nqz-1,no):: eq
        real(8),intent(inout):: rmu,rmuOLD,esterr

        integer(4) i
        integer(4),parameter:: itemax=1000
        logical(4),parameter:: T=.true.,F=.false.
        logical(4) flag
        real(8) eps,drmu,rmuL,rmuS,rnL,rnS,rmuM,rn,rnM,yac,rmuc,rnc,rmud,emax,emin

        if(RPA)then
           emax=maxval(eq)
           emin=minval(eq)
        end if
        yac=1.0d0/dble(nqx*nqy*nqz)
        if(esterr>1.0d-2)then
           eps= 1.0d-8
        else
           eps= esterr*1.0d-1
        end if
        drmu= abs(rmu-rmuOLD)*2.0d0
        if (drmu<eps*4.0d0) drmu= eps*4.0d0
        rmuOLD= rmu
        rmuL= rmu+drmu
        rmuS= rmu-drmu

        print '(1x,a,2f22.16)','rmuOLD,drmu=',rmuOLD,drmu
        do i=1,itemax
           if(RPA .and. rmuL < emin)then
              rmuL=emin
              rmu=rmuL
              rnL=0.0d0
           else if(RPA .and. rmuL>emax)then
              rmuL=emax
              rmu=rmuL
              rnL=dble(no)
           else
              rmu= rmuL
              rnL=0.0d0
              call makern(nqx,nqy,nqz,no,temp,rnL,eq,yac,rmu,RPA)
           end if
           print '(1x,a,2f22.16)','rmuL,rnL=   ',rmu,rnL
           if (rnL<rfill) then
              if(abs(rfill-rnL)>2.0d0)then
                 rmuL= rmuL +1.0d0
              else
                 rmuL= rmuL +drmu
              end if
           else
              exit
           end if
           if(i==itemax)then
              print*,'too many!'
              stop
           end if
        end do

        do i=1,itemax
           if(RPA .and. rmuS > emax )then
              rnS=dble(no)
              rmuS=emax
              rmu= rmuS
           else if(RPA .and. rmuS < emin)then
              rnS=0.0d0
              rmuS=emin
              rmu= rmuS
           else
              rmu= rmuS
              rnS=0.0d0
              call makern(nqx,nqy,nqz,no,temp,rnS,eq,yac,rmu,RPA)
           end if
           print '(1x,a,2f22.16)','rmuS,rnS=   ',rmuS,rnS
           if (rnS>rfill) then
              if(abs(rnS-rfill)>2.0d0)then
                 rmuS= rmuS -1.0d0
              else
                 rmuS= rmuS -drmu
              end if
           else
              exit
           end if
           if(i==itemax)then
              print*,'too many!'
              stop
           end if
        end do

        if(.false.)then !secant method
           do i=1,itemax
              rmuM= (rmuL+rmuS)*0.5d0
              if (abs(rmuL-rmuM)<eps .and. abs(rnM-rfill)<eps)exit
              rmu= rmuM
              rnM=0
              call makern(nqx,nqy,nqz,no,temp,rnM,eq,yac,rmu,RPA)
              print '(1x,a,2f22.16)','rmuM,rnM=   ',rmuM,rnM
              if (rnM>rfill) then
                 rmuL= rmuM
                 rnL= rnM
              else if (rnM<rfill) then
                 rmuS= rmuM
                 rnS= rnM
              else
                 rnS=-1
                 exit
              end if
              if(i==itemax)then
                 print*,'too many!'
                 stop
              end if
           end do
           if(rnS==-1)then
              rmu=rmuM
           else
              rmu= (rmuL*(rfill-rnS)+rmuS*(rnL-rfill))/(rnL-rnS)
           end if
        else !Brent's method
           rnS=rnS-rfill
           rnL=rnL-rfill
           rnc=rnS
           rmuc=rmuS
           rmud=0.0d0
           flag=F
           do i=1,itemax
              if(rnc/=rnS .and. rnc/=rnL)then
                 rmuM=(rmuL*rnS*rnc*(rnS-rnc)+rmuS*rnL*rnc*(rnc-rnL)+rmuc*rnL*rnS*(rnL-rnS))&
                      /((rnL-rnS)*(rnc-rnL)*(rnc-rnS))
              else
                 rmuM=rmuL-rnL*(rmuL-rmuS)/(rnL-rnS)
              end if
              if((0.25d0*(3.0d0*rmuS+rmuL) > rmuM .or. rmuM > rmuL) .or. &
                   (flag .and. abs(rmuM-rmuL) >= abs(rmuL-rmuc)*0.5d0) .or. &
                   (flag.eqv.F .and. abs(rmuM-rmuL) >= abs(rmuc-rmud)*0.5d0) .or. &
                   (flag .and. abs(rmuL-rmuc)<1.0d-8) .or. &
                   (flag.eqv.F .and. abs(rmuc-rmud)<1.0d-8))then
                 rmuM= (rmuL+rmuS)*0.5d0
                 flag=T
              else
                 flag=F
              end if
              if(abs(rmuL-rmuM)<eps)exit
              rmu=rmuM
              rnM=0
              call makern(nqx,nqy,nqz,no,temp,rnM,eq,yac,rmu,RPA)
              print '(1x,a,2f22.16,l2)','rmuM,rnM=   ',rmuM,rnM,flag
              rmud=rmuc
              rmuc=rmuL
              rnc=rnL

              if(rnS*(rnM-rfill)<0)then
                 rmuL=rmuM
                 rnL=rnM-rfill
              else
                 rmuS=rmuM
                 rnS=rnM-rfill
              end if

              if(abs(rnS)<abs(rnL))then
                 rmuM=rmuL
                 rmuL=rmuS
                 rmuS=rmuM
                 rnM=rnL
                 rnL=rnS
                 rnS=rnM
              end if
           end do
           if(rnL==rnS)then
              rmu=(rmuS+rmuL)*0.5d0
           else
              rmu= (rmuS*rnL-rmuL*rnS)/(rnL-rnS)
           end if
        end if

        rn=0.0d0
        call makern(nqx,nqy,nqz,no,temp,rn,eq,yac,rmu,RPA)
        print *,'rmu  =',rmu
        print *,'rn   =',rn
        print *,'temp = ',temp
      
        return
      contains
        !
        !> this subroutine calculate the number of electron in FLEX calculation
        ! 
        subroutine makern(nqx,nqy,nqz,no,temp,rn,eq,yac,rmu,RPA)
          implicit none
          integer(4),intent(in):: nqx,nqy,nqz,no
          logical(4),intent(in):: RPA
          real(8),intent(in):: temp,yac,rmu
          real(8),intent(in),dimension(0:nqx-1,0:nqy-1,0:nqz-1,no):: eq
          real(8),intent(out):: rn

          integer(4) i,j,k,n,l
          real(8) tmp, dosk

          tmp=0.0d0
          !$omp parallel private(l,n) reduction(+: tmp)
          do l=1,no
             do n=0,nqz-1
                !$omp do private(i)
                do j=0,nqy-1
                   do i=0,nqx-1
                      if((eq(i,j,n,l)-rmu)>5.0d1*temp)then
                         continue
                      else if((eq(i,j,n,l)-rmu)<-5.0d1*temp)then
                         tmp=tmp+1.0d0
                      else
                         tmp=tmp+1.0d0/(exp((eq(i,j,n,l)-rmu)/temp)+1.0d0)
                      end if
                   end do
                end do
                !$omp end do
             end do
          end do
          !$omp end parallel
          rn=tmp*yac
        end subroutine makern
      end subroutine mkrmu
