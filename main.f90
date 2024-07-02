      program main
        implicit none
        integer(4),parameter:: nqx=32,nqy=32,nqz=2,nr=581*2,no=5
        integer(4),parameter:: Enmax=150,site_num=1, const=1, dq=1
        integer(4),parameter:: nos(site_num)=(/5/), sw_band=1, inflag=1
        real(8),parameter::U=0.8d0,J=0.1d0,temp=2.59d-2,idelta=3.0d-2,mass=1.0d0
        real(8),parameter:: U_weight(site_num)=(/1.0d0/),rfill=2.9375d0,Emax=3.0d0
        real(8),parameter:: cntw=1.0d-2
        character(len=60),parameter:: hamname='000AsP.input'
        logical(4),parameter:: sw_im=.false.

        real(8),parameter:: pi=4.0d0*atan(1.0d0)
        real(8),dimension(0:nqx-1,0:nqy-1,0:nqz-1,no):: eq,ffermi
        complex(8),dimension(0:nqx-1,0:nqy-1,0:nqz-1,no,no):: uni,en
        complex(8),dimension(0:Enmax):: trchi,trchis
        complex(8),dimension(0:Enmax,no+2):: chisorb
        real(8),allocatable,dimension(:,:):: S
        integer(4) l,nco,qx,qy,qz,Nk,iw
        real(8) dE,rmu,rmuOld,esterr,lsp,t1,t2
        complex(8) w

        call cpu_time(t1)
        nco=sum(nos(:)**2)
        rmu=5.8d0
        rmuOLD=rmu-0.1d0
        esterr=1.0d0
        dE=Emax/Enmax
        Nk=nqx*nqy*nqz
        print *,'U=',U,'J=',J
        print *,'k-mesh = ',nqx,nqy,nqz
        print *,'no =',no,'delta = ',idelta
        allocate(S(nco,nco))
        call get_S(nco,S,U,J)
        call mkhm(nqx,nqy,nqz,no,nr,eq,en,uni,.false.,hamname,inflag,mass)
        call mkrmu(nqx,nqy,nqz,no,temp,rfill,eq,rmu,rmuOLD,esterr,.false.)
        ffermi=0.5d0-tanh(0.5d0*(eq-rmu)/temp)*0.5d0
        select case(sw_band)
        case(1,2) !1 G>X>M, 2G>X>M>G
           lsp=-1
           do qx=0,nqx/2,dq
              lsp=lsp+1
              do l=1,no
                 write(40,*)lsp,eq(qx,0,0,l)
              end do
           end do
           close(40)

           lsp=-1
           do qx=0,nqx/2,dq
              qy=0
              qz=0
              lsp=lsp+1
              call calc_chi_q(qx,qy,qz,trchi,trchis,chisorb)
              do iw=0,Enmax
                 if(sw_im)then
                    write(50,*)lsp,iw*de,dble(trchi(iw))
                    write(60,*)lsp,iw*de,dble(trchis(iw))
                 else
                    write(50,*)lsp,iw*de,aimag(trchi(iw))
                    write(60,*)lsp,iw*de,aimag(trchis(iw))
                    do l=1,no+2
                       write(70+l,*)lsp,iw*de,aimag(chisorb(iw,l))
                    end do
                 end if
              end do
              if(Enmax>1)then
                 write(50,*)''
                 write(60,*)''
                 do l=1,no+2
                    write(70+l,*)''
                 end do
              end if
           end do

           do qy=1,nqy/2,dq
              qx=nqx/2
              qz=0
              lsp=lsp+1
              call calc_chi_q(qx,qy,qz,trchi,trchis,chisorb)
              do iw=0,Enmax
                 if(sw_im)then
                    write(50,*)lsp,iw*de,dble(trchi(iw))
                    write(60,*)lsp,iw*de,dble(trchis(iw))
                 else
                    write(50,*)lsp,iw*de,imag(trchi(iw))
                    write(60,*)lsp,iw*de,imag(trchis(iw))
                    do l=1,no+2
                       write(70+l,*)lsp,iw*de,aimag(chisorb(iw,l))
                    end do
                 end if
              end do
              if(Enmax>1)then
                 write(50,*)''
                 write(60,*)''
                 do l=1,no+2
                    write(70+l,*)''
                 end do
              end if
           end do

           if(sw_band==2)then
              do qx=nqx/2-1,0,-dq
                 qy=qx
                 qz=0
                 lsp=lsp+sqrt(2.0d0)
                 call calc_chi_q(qx,qy,qz,trchi,trchis,chisorb)
                 do iw=0,Enmax
                    if(sw_im)then
                       write(50,*)lsp,iw*de,dble(trchi(iw))
                       write(60,*)lsp,iw*de,dble(trchis(iw))
                    else
                       write(50,*)lsp,iw*de,imag(trchi(iw))
                       write(60,*)lsp,iw*de,imag(trchis(iw))
                       do l=1,no+2
                          write(70+l,*)lsp,iw*de,aimag(chisorb(iw,l))
                       end do
                    end if
                 end do
                 if(Enmax>1)then
                    write(50,*)''
                    write(60,*)''
                    do l=1,no+2
                       write(70+l,*)''
                    end do
                 end if
              end do
           end if
           close(50)
           close(60)
           do l=1,no
              close(70+l)
           end do
        case(3) !X>X (z-axis)
           lsp=-1
           do qz=0,nqz-1
              qx=nqx/2
              qy=0
              lsp=lsp+1
              call calc_chi_q(qx,qy,qz,trchi,trchis,chisorb)
              do iw=0,Enmax
                 if(sw_im)then
                    write(50,*)lsp,iw*de,dble(trchi(iw))
                    write(60,*)lsp,iw*de,dble(trchis(iw))
                 else
                    write(50,*)lsp,iw*de,aimag(trchi(iw))
                    write(60,*)lsp,iw*de,aimag(trchis(iw))
                    do l=1,no+2
                       write(70+l,*)lsp,iw*de,aimag(chisorb(iw,l))
                    end do
                 end if
              end do
              if(Enmax>1)then
                 write(50,*)''
                 write(60,*)''
                 do l=1,no+2
                    write(70+l,*)''
                 end do
              end if
           end do
        case default
           w=cmplx(cntw,idelta)
           call calc_qmap(w)
        end select

        deallocate(S)
        call cpu_time(t2)
        print *,'Time = ',t2-t1,'sec.'
      contains
        subroutine calc_chi_q(qx,qy,qz,trchi,trchis,chisorb)
          integer(4),intent(in):: qx,qy,qz
          complex(8),intent(out),dimension(0:Enmax):: trchi,trchis
          complex(8),intent(out),dimension(0:Enmax,no+2):: chisorb

          integer(4) l1,l2,m1,m2,iw,i,j,ini1,ini2,l,m
          integer(4) iq(0:nqx-1),jq(0:nqy-1),kq(0:nqz-1)
          complex(8),dimension(nco,nco)::chi,chis,dumm
          complex(8) w,tmp

          integer(4) ipiv(nco),info
          complex(8),dimension(nco):: work

          !$omp parallel
          !$omp do
          do i=0,nqx-1
             if(i+qx<nqx)then
                iq(i)=i+qx
             else
                iq(i)=i+qx-nqx
             end if
          end do
          !$omp end do
          !$omp do
          do i=0,nqy-1
             if(i+qy<nqy)then
                jq(i)=i+qy
             else
                jq(i)=i+qy-nqy
             end if
          end do
          !$omp end do
          !$omp do
          do i=0,nqz-1
             if(i+qz<nqz)then
                kq(i)=i+qz
             else
                kq(i)=i+qz-nqz
             end if
          end do
          !$omp end do
          !$omp workshare
          trchis(:)=0.0d0
          trchi(:)=0.0d0
          chisorb(:,:)=0.0d0
          !$omp end workshare

          !$omp do private(iw,w,l,l1,l2,m,m1,m2,i,j,ini1,ini2,chi,chis,dumm,ipiv,info,work,tmp)
          do iw=0,Enmax
             if(sw_im)then
                w=cmplx(0.0d0,dble(2*iw)*pi*temp)
             else
                w=cmplx(iw*dE,idelta) !+iw*dE*5.0d-2)
             end if

             ini1=0
             do i=1,site_num
                ini2=0
                do j=1,site_num
                   do l1=1,nos(i)
                      do l2=1,nos(i)
                         l=(l1-1)*nos(i)+l2+ini1
                         do m1=1,nos(j)
                            do m2=1,nos(j)
                               m=(m1-1)*nos(j)+m2+ini2
                               chi(m,l)=irr_chi(iq,jq,kq,w,l1,l2,m1,m2)
                            end do
                         end do
                      end do
                   end do
                   ini2=ini2+nos(j)*nos(j)
                end do
                ini1=ini1+nos(i)*nos(i)
             end do
             chis(:,:)=0.0d0
             dumm(:,:)=0.0d0
             do l=1,nco
                do l1=1,nco
                   !$omp simd
                   do m=1,nco
                      dumm(m,l)=dumm(m,l)-chi(m,l1)*S(l1,l)
                   end do
                   !$omp end simd
                end do
                dumm(l,l)=dumm(l,l)+1.0d0
             end do
             call zgetrf(nco,nco,dumm,nco,ipiv,info)
             call zgetri(nco,dumm,nco,ipiv,work,nco,info)
             do l=1,nco
                do l1=1,nco
                   !$omp simd
                   do m=1,nco
                      chis(m,l)=chis(m,l)+dumm(m,l1)*chi(l1,l)
                   end do
                   !$omp end simd
                end do
             end do
             ini1=0
             do i=1,site_num
                ini2=0
                do j=1,site_num
                   if(i==j)then
                      do l1=1,nos(i)
                         do l2=1,nos(i)
                            if(l1==l2)then
                               l=(l1-1)*nos(i)+l2+ini1
                               do m1=1,nos(j)
                                  do m2=1,nos(j)
                                     m=(m1-1)*nos(j)+m2+ini2
                                     if(m1==m2)then ! .and. (l1==m1))then
                                        trchis(iw)=trchis(iw)+chis(m,l)
                                        trchi(iw)=trchi(iw)+chi(m,l)
                                        if(l1==2 .and. m1==3)chisorb(iw,no+1)=chis(m,l)
                                        if(l1==2 .and. m1==4)chisorb(iw,no+2)=chis(m,l)
                                        if(l1==m1)chisorb(iw,l1)=chis(m,l)
                                     end if
                                  end do
                               end do
                            end if
                         end do
                      end do
                   end if
                end do
             end do
             tmp=trchis(iw)
             if(abs(aimag(trchis(iw)))<1.0d-6)trchis(iw)=cmplx(dble(tmp),0.0d0)
             if(abs(dble(trchis(iw)))<1.0d-6)trchis(iw)=cmplx(0.0d0,aimag(tmp))
             tmp=trchi(iw)
             if(abs(aimag(trchi(iw)))<1.0d-6)trchi(iw)=cmplx(dble(tmp),0.0d0)
             if(abs(dble(trchi(iw)))<1.0d-6)trchi(iw)=cmplx(0.0d0,aimag(tmp))
          end do
          !$omp end do
          !$omp end parallel
        end subroutine calc_chi_q

        subroutine calc_qmap(w)
          complex(8),intent(in):: w

          integer(4) i,j,l,m,l1,l2,m1,m2,qx,qy,qz,ini1,ini2
          integer(4) iq(0:nqx-1),jq(0:nqy-1),kq(0:nqz-1)
          complex(8),dimension(nco,nco)::chi,chis,dumm
          complex(8),dimension(0:nqx-1):: trchi,trchis
          complex(8) tmp

          integer(4) ipiv(nco),info
          complex(8),dimension(nco):: work

          qy=0
          do i=0,nqy-1
             if(i+qy<nqy)then
                jq(i)=i+qy
             else
                jq(i)=i+qy-nqy
             end if
          end do

          do qz=0,nqz-1,dq
             trchis(:)=0.0d0
             trchi(:)=0.0d0
             do i=0,nqz-1
                if(i+qz<nqz)then
                   kq(i)=i+qz
                else
                   kq(i)=i+qz-nqz
                end if
             end do
             !$omp parallel do private(qx,iq,l,l1,l2,m,m1,m2,i,j,ini1,ini2,chi,chis,dumm,ipiv,info,work,tmp)
             do qx=0,nqx-1,dq
                do i=0,nqx-1
                   if(i+qx<nqx)then
                      iq(i)=i+qx
                   else
                      iq(i)=i+qx-nqx
                   end if
                end do
                ini1=0
                do i=1,site_num
                   ini2=0
                   do j=1,site_num
                      do l1=1,nos(i)
                         do l2=1,nos(i)
                            l=(l1-1)*nos(i)+l2+ini1
                            do m1=1,nos(j)
                               do m2=1,nos(j)
                                  m=(m1-1)*nos(j)+m2+ini2
                                  chi(m,l)=irr_chi(iq,jq,kq,w,l1,l2,m1,m2)
                               end do
                            end do
                         end do
                      end do
                      ini2=ini2+nos(j)*nos(j)
                   end do
                   ini1=ini1+nos(i)*nos(i)
                end do

                chis(:,:)=0.0d0
                dumm(:,:)=0.0d0
                do l=1,nco
                   do l1=1,nco
                      !$omp simd
                      do m=1,nco
                         dumm(m,l)=dumm(m,l)-chi(m,l1)*S(l1,l)
                      end do
                      !$omp end simd
                   end do
                   dumm(l,l)=dumm(l,l)+1.0d0
                end do
                call zgetrf(nco,nco,dumm,nco,ipiv,info)
                call zgetri(nco,dumm,nco,ipiv,work,nco,info)
                do l=1,nco
                   do l1=1,nco
                      do m=1,nco
                         chis(m,l)=chis(m,l)+dumm(m,l1)*chi(l1,l)
                      end do
                   end do
                end do

                ini1=0
                do i=1,site_num
                   ini2=0
                   do j=1,site_num
                      if(i==j)then
                         do l1=1,nos(i)
                            do l2=1,nos(i)
                               l=(l1-1)*nos(i)+l2+ini1
                               do m1=1,nos(j)
                                  do m2=1,nos(j)
                                     m=(m1-1)*nos(j)+m2+ini2
                                     if((l1==l2) .and. (m1==m2))then ! .and. (l1==m1))then
                                        trchis(qx)=trchis(qx)+chis(m,l)
                                        trchi(qx)=trchi(qx)+chi(m,l)
                                     end if
                                  end do
                               end do
                            end do
                         end do
                      end if
                      ini2=ini2+nos(j)*nos(j)
                   end do
                   ini1=ini1+nos(i)*nos(i)
                end do
                tmp=trchis(qx)
                if(abs(aimag(trchis(qx)))<1.0d-6)trchis(qx)=cmplx(dble(tmp),0.0d0)
                if(abs(dble(trchis(qx)))<1.0d-6)trchis(qx)=cmplx(0.0d0,aimag(tmp))
                tmp=trchi(qx)
                if(abs(aimag(trchi(qx)))<1.0d-6)trchi(qx)=cmplx(dble(tmp),0.0d0)
                if(abs(dble(trchi(qx)))<1.0d-6)trchi(qx)=cmplx(0.0d0,aimag(tmp))
             end do
             !$omp end parallel do
             if(qy==0)then
                do qx=0,nqx-1,dq
                   write(200,*)qx/dq,qz/dq,dble(trchi(qx)),aimag(trchi(qx))
                   write(210,*)qx/dq,qz/dq,dble(trchis(qx)),aimag(trchis(qx))
                end do
                write(200,*)''
                write(210,*)''
             end if
          end do
          close(200)
          close(210)

        end subroutine calc_qmap

        complex function irr_chi(iq,jq,kq,w,o1,o2,o3,o4)
          integer(4),intent(in):: o1,o2,o3,o4
          integer(4),intent(in):: iq(0:nqx-1),jq(0:nqy-1),kq(0:nqz-1)
          complex(8),intent(in):: w

          integer(4) i,j,k,l,m
          real(8) eps
          complex(8) chi,unitmp

          chi=0.0d0
          eps=idelta*1.0d-3
          do l=1,no
             do m=1,no
                do k=0,nqz-1
                   do j=0,nqy-1
                      do i=0,nqx-1
                         unitmp=uni(iq(i),jq(j),kq(k),o1,l)*conjg(uni(iq(i),jq(j),kq(k),o3,l))*uni(i,j,k,o4,m)*conjg(uni(i,j,k,o2,m))
                         if(abs(w)==0.0d0 .and. abs(eq(i,j,k,m)-eq(iq(i),jq(j),kq(k),l))<1.0d-9)then
                            chi=chi+unitmp*ffermi(i,j,k,m)*(1.0d0-ffermi(i,j,k,m))/temp
                         else if(abs(ffermi(iq(i),jq(j),kq(k),l)-ffermi(i,j,k,m))>eps)then
                            chi=chi+unitmp*(ffermi(iq(i),jq(j),kq(k),l)-ffermi(i,j,k,m))/(w+eq(i,j,k,m)-eq(iq(i),jq(j),kq(k),l))
                         end if
                      end do
                   end do
                end do
             end do
          end do
          irr_chi=chi/Nk
        end function irr_chi

        subroutine get_S(nco,S,U,J)
          integer(4),intent(in):: nco
          real(8),intent(in):: U,J
          real(8),intent(out):: S(nco,nco)
          integer(4) l1,l2,m1,m2,l,m,i,nsq,init

          init=0
          S=0.0d0
          do i=1,site_num
             nsq=nos(i)*nos(i)
             do l=1,nsq
                l1=(l-1)/nos(i)+1
                l2=mod(l-1,nos(i))+1
                do m=1,nsq
                   m1=(m-1)/nos(i)+1
                   m2=mod(m-1,nos(i))+1
                   if((l1==l2) .and. (m1==m2))then
                      if(l1==m1)then
                         S(m+init,l+init)=U
                      else
                         S(m+init,l+init)=J
                      end if
                   else if((l1==m1) .and. (l2==m2))then
                      S(m+init,l+init)=U-2*J
                   else if((l1==m2) .and. (l2==m1))then
                      S(m+init,l+init)=J
                   end if
                end do
             end do
             init=init+nos(i)*nos(i)
          end do
        end subroutine get_S
      end program main
