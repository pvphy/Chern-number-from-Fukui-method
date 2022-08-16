module array

  implicit none
  double precision, allocatable::mi(:),th(:),ph(:),rwork(:),evl_s(:)
  integer,allocatable::x1(:),y1(:),x(:),y(:)	
  complex*16,allocatable:: h(:,:),bb(:,:),work(:),mat(:,:),h_(:,:)
  complex*16,allocatable::mat_1(:,:),mat_2(:,:),mat_3(:,:),vec(:) 
  complex*16,allocatable::mat_4(:,:),mat_5(:,:),mat_6(:,:) 

  !complex*16,allocatable,dimension(:,:):: h3,h1,h4,h2
  complex*16,allocatable::mat_1_up(:,:),mat_2_up(:,:),mat_3_up(:,:)
  complex*16,allocatable::mat_4_up(:,:),mat_5_up(:,:),mat_6_up(:,:)
  complex*16,allocatable::mat_4_down(:,:),mat_5_down(:,:),mat_6_down(:,:)


  complex*16,allocatable::mat_1_down(:,:),mat_2_down(:,:),mat_3_down(:,:),mat_temp(:,:),mat_temp_down(:,:)
end module array

module global

  implicit none
  double precision::t,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,spec_f,spec,kx,ky
  double precision::pi,lambda1,lambda,mag_field,eta
  integer::i_1,i_2,i_3,m,i,j,k_x,k_y,l,b1,d,d1,d2,ie,sigma,a1,lx1,ly1,jx,jy,lx,ly,nx,ny,dkx,dky
  complex*16::ch_no
endmodule global



program chern_no
  use array
  use global	
  integer::grid
  double precision::k2,temp,get_mu_s

  print*, 'give dimension'
  read*,d
  print*, 'hopping parameter t'
  read*,t
  print*,'lambda='
  read*,lambda
  print*,'B='
  read*,mag_field
  !print*, 'Ua='
  !read*,ua
  !print*, 'give Ub'
  !read*,ub
  !print*, 'give Uc'
  !read*,uc
  ua=0
  ub=0
  uc=0
  d2=(3*d**2)/4  ! no of sites
  d1=2*d2  
  ! print*,'unit cell dimension',d/2,d/2              
  !print*,'unit cell dimension of Lieb lattice',d/2,d/2              
  allocate(x1((d)**2),y1((d)**2),x(d2),y(d2),mat(d1,d1),h(d1,d1),mi(d2),th(d2),ph(d2))
  allocate(rwork(3*d1-2),evl_s(d1),work(2*d1-1))
  allocate(h_(d1,d1))
  allocate(mat_temp((d/2)**2,d1),mat_temp_down((d/2)**2,d1))
  call matgen
  call dos
  call bloch_states
  !call check_states
  
  do i=1,18
    call chern_input
    call chern_number
    write(77,*)i,ch_no
  enddo

      print*,'chern number for bottom_band'
      mat_temp=mat_1
      call chern_number
      print*,ch_no
      print*,
      mat_temp=mat_2
      call chern_number
      print*,ch_no
      print*,

      print*,'chern number for flat_band'
      mat_temp=mat_3
      call chern_number
      print*,ch_no
      print*,
      mat_temp=mat_4
      call chern_number
      print*,ch_no
      print*,


      print*,'chern number for top_band'
      mat_temp=mat_5
      call chern_number
      print*,ch_no
      print*,

      mat_temp=mat_6    
      call chern_number
      print*,ch_no
      
      !print*,'1_up_down'
     ! mat_temp=mat_1_up
     ! call chern_number
     ! print*,ch_no
     ! print*,
     ! mat_temp=mat_1_down
     ! call chern_number
     ! print*,ch_no




endprogram chern_no

!********************************************************
!********************Matrix  Generation******************
!********************************************************

subroutine matgen
  !use mtmod
    use global
    use array
    implicit none 
    integer::k,px,py,q,nn,lda,lwmax,info,lwork,a,b,ii,id1,ji,jd,i1,ij,ik,s,p,r1,r2
    integer::site_b(2*(d/2)**2),site_c(2*(d/2)**2),a_1
    double precision::ortho,de,tm1,ua22,xx,yy,zz!,mi,th,ph
    !complex*16::rr,h_up(d1/2,d1/2),h_down(d1/2,d1/2)
    pi=4*atan(1.0)
    evl_s=0.0d0
    !22 format((24f8.4,2x))
    lda=d1
    lwmax=d1
    lwork=(2*d1-1)

    open(unit=11,file='level2.dat')
    open(unit=15,file='liebmatrix.dat')
    !sites of the square lattice 
    x1=0
    y1=0
    p=0			
    do l=1,d
    do k=1,d
    
    p=p+1
    x1(p)=l
    y1(p)=k

    if(((mod((y1(p)),2)==0) .and.(mod((x1(p)),2)==0)))then
    cycle 
    endif
    write(11,*) x1(p),y1(p)
    enddo
    enddo
    rewind(11)

    ! lieb lattice co-ordinate

    open(unit=11,file='level2.dat')
    x=0
    y=0
    do i=1,d2
    

    read(11,*) r1,r2
    x(i)=r2
    y(i)=r1
    !print*,i,x(i),y(i)
    enddo
    
    m_sum=0.0d0
    h=complex(0.0d0,0.0d0)
    h_=complex(0.0d0,0.0d0)
   




    do l=1,d2	   
      i=x(l)
      j=y(l)
      ii=1
      id1=-1
      ji=1                
      jd=-1 
          
      if (i.eq.1) id1=-i+(d)
      if (i.eq.(d)) ii=1-i
      if (j.eq.1) jd=-j+(d)
      if (j.eq.(d)) ji=1-j
            
    do k=1,d2	
    if(((x(k).eq.(i+ii)).and.(y(k).eq.j)).or.((x(k).eq.i).and.(y(k).eq.(j+ji))))then
              a=2*l-1
              b=2*k-1
              h(a,b)=t
              h(a+1,b+1)=t
            !	print*,a,b,a+2,b+2
    endif
    if(((x(k).eq.(i+id1)).and.(y(k).eq.j)).or.((x(k).eq.i).and.(y(k).eq.(j+jd))))then
              a=2*l-1
              b=2*k-1
              h(a,b)=t
              h(a+1,b+1)=t	
        !print*,a,b,a+2,b+2	 					 
      endif	
      
      
      
      
      if((x(k).eq.i).and.(y(k).eq.j))then
      
      if((mod(x(l),2).gt.0 ) .and.(mod(y(l),2)==0 ))then
      u=ub
      endif
      if((mod(x(l),2).gt.0 ) .and.(mod(y(l),2).gt.0 ))then
      u=ua
      endif
      if((mod(x(l),2)==0 ) .and.(mod(y(l),2).gt.0 ))then
      u=uc
      endif
      
      ! xx=mi(l)*sin(th(l))*cos(ph(l))
      ! yy=mi(l)*sin(th(l))*sin(ph(l))
      ! zz=mi(l)*cos(th(l))
      ! write(32,*) xx,yy,zz
      
      a=2*l-1
      b=2*k-1				 		
        
      !mi(l)=m
      !mi(l)=1
      !th(l)=th1
      ! th(l)=pi/4
      !ph(l)=ph1
      !print*,u	  
      h(a,b)=(u/2.0d0)*(-mi(l)*cos(th(l)))
      h(a,b+1)=(-u/2.0d0)*mi(l)*sin(th(l))*complex(cos(ph(l)),-sin(ph(l)))
      h(a+1,b)=conjg(h(a,b+1))
      h(a+1,b+1)=(u/2.0d0)*(-mi(l)*(-cos(th(l))))   
            
      m_sum=m_sum+((mi(l))**2.0d0)*(U/4.0d0)	 
      !print*,m_sum   
    endif
    enddo
    enddo
    
 !************** spin orbit coupling***********************************

  do l=1,d2

    !if((x(l)==2).or.(x(l)==8)) print*,a,b

    do k=1,d2 

      
      if((abs(x(l)-x(k))==1).and.(abs(y(l)-y(k))==1))then
        
       ! print*,l,k
        a=2*l-1
        b=2*k-1
       !if(x(l)==1) print*,l,'here'
       !if(mod(y(l),2)==0) print*,2*l-1
        if(a.lt.b)then
          if((mod(y(l),2)==0).or.(mod(y(l),2)==0))then
            h(a,b)=-(x(l)-x(k))*(lambda)*complex(0,1)
            h(a+1,b+1)=conjg(h(a,b))
            h(b,a)=conjg(h(a,b))
            h(b+1,a+1)=conjg(h(a+1,b+1))
           ! print*,a,b    !,-(x(l)-x(k)),y(l)-y(k)
          else
            h(a,b)= (x(l)-x(k))*(lambda)*complex(0,1)
            h(a+1,b+1)=conjg(h(a,b))
            h(b,a)=conjg(h(a,b))
            h(b+1,a+1)=conjg(h(a+1,b+1))
           ! print*,a,b   !,x(l)-x(k),y(l)-y(k)
          endif
        endif

       ! if(a.gt.b)then
        !  h(a,b)=(lambda)*complex(0,-1)
        !  h(a+1,b+1)=(lambda)*complex(0,1)
         ! print*,a,b,a+1,b+1
        !endif

      endif        
    enddo
  enddo

 !**************************************************************

!boundary condition for soc

  do l=1,d2
    
    if((x(l)==d))then
     a=2*l-1
     b=a+2
     a_1=2*l-1-(d+(d/2)-1)*2
     !print*,a,b
     !print*,a_1
     h(a,b)=-(lambda)*complex(0,1)
     h(b,a)=conjg(h(a,b))
     h(a+1,b+1)=conjg(h(a,b))
     h(b+1,a+1)=conjg(h(a+1,b+1))
     ! print*,a,b,'3'
     if(a_1.lt.0)then
      a_1=3*((d**2)/4)-(d/2)+1
      b=2*a_1-1
      h(a,b)=(lambda)*complex(0,1)
      h(b,a)=conjg(h(a,b))
      h(a+1,b+1)=conjg(h(a,b))
      h(b+1,a+1)=conjg(h(a+1,b+1))
      !print*,a,b,'1'
      ! endif
     elseif(a_1.gt.0)then
     h(a,a_1)=(lambda)*complex(0,1)
     h(a_1,a)=conjg(h(a,a_1))
     h(a+1,a_1+1)=conjg(h(a,a_1))
     h(a_1+1,a+1)=conjg(h(a+1,a_1+1))
     ! print*,a,a_1,'2'
     endif
     
     !print*,a,b,'3'
    endif
 
    if((y(l)==d))then
      a=2*l-1
      b=2*x(l)+1
      a_1=b-4
      h(a,b)=(lambda)*complex(0,1)
      h(b,a)=conjg(h(a,b))
      h(a+1,b+1)=conjg(h(a,b))
      h(b+1,a+1)=conjg(h(a+1,b+1))
       ! print*,a,b,'4'

        if(a_1.gt.0)then
        h(a,a_1)=(lambda)*complex(0,-1)
        h(a_1,a)=conjg(h(a,a_1))
        h(a+1,a_1+1)=conjg(h(a,a_1))
        h(a_1+1,a+1)=conjg(h(a+1,a_1+1))
        !print*,a,a_1,'5'
        endif
        

    endif
      
  enddo
      
  

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    



    
!***********************************************************************!	
    do i=1,d1

      if(mod(i,2).ne.0)	h(i,i)=-0.5*mag_field
      if(mod(i,2).eq.0)	h(i,i)=0.5*mag_field

    enddo
    
!*************************************************************************!    
    h_=h
    
    !print*,h(3,5),h(4,6)
    !do i=1,d1 
     ! do j=1,d1

      !if((mod(i,2).eq.0).and.(mod(j,2).eq.0))then
      ! h_up(i/2,j/2)=h(i,j)
         
      ! endif
      !if((mod(i,2).ne.0).and.(mod(j,2).ne.0))then        
      ! h_down((i+1)/2,(j+1)/2)=h(i,j)
          
      ! endif  
      ! write(34,*) i,j,h(i,j)
      ! if(h(i,j).ne.conjg(h(j,i)))  print*,i,j,'not hermitian'
     ! enddo
    
    !enddo
     !print*,h_up(2,5),h_down(2,5)
    !print*,  shape(evl_s)
      call zheev('v','u',d1,h,lda,evl_s,work,lwork,rwork,info)
    
      if(info.ne.0)then
         print*,'algorithm failed'  
      endif 
    
    do i=1,d1
      write(726,*) i, evl_s(i)
    enddo



    
    !  print*,evl_s(1),evl_s(d1)

endsubroutine matgen


!******************************************************************!
!************Function for chemical potential***********************!
!******************************************************************!	

double precision function get_mu_s(fill,temp2)
          use array
          !    use input
          use global
          implicit none

          double precision::f0,f,fL2,fR,mR,mL,m_d,temp2
          integer::fill
    !      double precision fill

          mR = maxval(evl_s)       !right-side chemical potential
          fr=0.0d0
          do i=1,d1
    fr=fr+(1.0d0/(exp((evl_s(i)-mR)/Temp2)+1.0d0))
          end do
          mL = minval(evl_s)       !left-side chemical potential
          fL2=0.0d0
          do i=1,d1
    fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/Temp2)+1.0d0))
          end do
          m_d = 0.5d0*(mL+mR)    !middle chemical potential
          f=0.0d0
          do i=1,d1
    f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
          end do
          !print*,f,fill
          do while(abs(f-fill).ge.1e-8)
    m_d = 0.5d0*(mL+mR)
          f=0.0d0
          do i=1,d1
    f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
          end do
          if(f.gt.fill)then
          !if middle filling is above target, make it the new right bound.
          mR = m_d
          fR = f
          elseif(f.lt.fill)then
          !if middle filling is below target, make it the new left bound.
          mL = m_d
          fR = f
          endif
          enddo        
          !Return the middle value
          get_mu_s = m_d
          return
end function get_mu_s	
          
!******************************************************************!
!************subroutine to calculate Spectral function*************!
!******************************************************************! 

subroutine spec_func
  use array
  use global
  implicit none
  double precision::mat_e,c1 !c2
  complex*16::c2
  complex*16::mat_sum(d1,d1)
  pi=4.0*atan(1.0)
  ! nk=0     
  !  ky=pi/2
  !  nk=nk+1
  !print*,kx,ky
  eta=0.025000d0
  
  !   print*,h(16,16)
  e=2.0d0  
  mat_sum=0.0d0
  do l=1,d2     
    do j=1,d2
     do b1=1,d1
       mat_e=((h(2*j,b1)*conjg(h(2*l,b1)))+(h(2*j-1,b1)*conjg(h(2*l-1,b1))))
       mat_sum(l,j)=mat_sum(l,j)+mat_e*((eta/(2*pi))/((e-evl_s(b1))**2+((eta/2)**2)))                     
      enddo
    enddo
  enddo
	    
  do k_x=0,d/2    !d2
    kx =-pi/2 +(pi/(d/2))*k_x
    ! spec=0.0d0 
    do k_y=0,d/2
      ky= -pi/2+(pi/(d/2))*k_y 
      spec_f=0.0d0 
      do l=1,d2     
        do j=1,d2
					            
          c1=(kx*(x(j)-x(l))+ky*(y(j)-y(l)))
          c2=complex(cos(c1),sin(c1))
          spec_f=spec_f+c2*mat_sum(l,j)     
					    
        enddo   
      enddo    
      write(120,*) kx,ky,spec_f/(d2)
      flush(120)
    enddo  ! ie loop       
  enddo
  !  enddo 
end subroutine spec_func


subroutine dos
  use array
  use global
  implicit none
  double precision::de,eta_
  integer::p
  pi=4.0*atan(1.0)
  e=-14.0d0
  de=0.01750d0
	eta_=0.050d0
  do ie=1,1600!no of interval bw bandwitdh of e
    e=e+de
    nsum=0.0d0
		
    do j=1,d1
      nsum=nsum+((eta_/pi)/((e-evl_s(j))**2+((eta_)**2)))
    enddo
		
    write(201,*) e,nsum/d1  
    flush(201)              				  
  enddo  
     
endsubroutine dos
    

subroutine bloch_states
  use global
  use array
  implicit none 
  integer::k,level,comp,p
  double precision::mat_e1,k_r,exp_k_r,norm,kx1(d1),ky1(d1),e_1(d1)
  double precision::bottom_band1((d/2)**2),flat_band1((d/2)**2),top_band1((d/2)**2),dummy_value,eps_
  double precision::bottom_band2((d/2)**2),flat_band2((d/2)**2),top_band2((d/2)**2)

  allocate(mat_1(((d/2)**2),d1),mat_2(((d/2)**2),d1),mat_3(((d/2)**2),d1),vec(d1))
  allocate(mat_4(((d/2)**2),d1),mat_5(((d/2)**2),d1),mat_6(((d/2)**2),d1))

  allocate(mat_1_up(((d/2)**2),d1),mat_2_up(((d/2)**2),d1),mat_3_up(((d/2)**2),d1))
  allocate(mat_4_up(((d/2)**2),d1),mat_5_up(((d/2)**2),d1),mat_6_up(((d/2)**2),d1))

  allocate(mat_1_down(((d/2)**2),d1),mat_2_down(((d/2)**2),d1),mat_3_down(((d/2)**2),d1))
  allocate(mat_4_down(((d/2)**2),d1),mat_5_down(((d/2)**2),d1),mat_6_down(((d/2)**2),d1))

  !vec=complex(0.0d0,0.0d0)
  eps_=1e-6
  dummy_value=-123456789.0
  bottom_band1=dummy_value
  bottom_band2=dummy_value
  flat_band1=dummy_value
  flat_band2=dummy_value
  top_band1=dummy_value  
  top_band2=dummy_value  

  ! print*,'here 1' 
  mat_1=complex(0.0d0,0.0d0)
  mat_2=complex(0.0d0,0.0d0)
  mat_3=complex(0.0d0,0.0d0) 
  mat_4=complex(0.0d0,0.0d0)
  mat_5=complex(0.0d0,0.0d0)
  mat_6=complex(0.0d0,0.0d0) 
  mat_1_up=complex(0.0d0,0.0d0) 
  mat_2_up=complex(0.0d0,0.0d0) 
  mat_3_up=complex(0.0d0,0.0d0) 
  mat_4_up=complex(0.0d0,0.0d0) 
  mat_5_up=complex(0.0d0,0.0d0) 
  mat_6_up=complex(0.0d0,0.0d0) 
  mat_1_down=complex(0.0d0,0.0d0) 
  mat_2_down=complex(0.0d0,0.0d0) 
  mat_3_down=complex(0.0d0,0.0d0) 
  mat_4_down=complex(0.0d0,0.0d0) 
  mat_5_down=complex(0.0d0,0.0d0) 
  mat_6_down=complex(0.0d0,0.0d0) 

  do m=1,d1
    e=evl_s(m)
    level=m
    i=0
    do k_y=1,(d/2)  !d2
     !  print*,'here 2'  
     ky =(2*pi/(d/2))*k_y
     do k_x=1,(d/2)
       kx=(2*pi/(d/2))*k_x
       ! print*,'for',p1,k1 
       k=k_x+(k_y-1)*(d/2)	 
       i=i+1   


       do lx=0,(d/2)-1  
         do ly=0,(d/2)-1 
           !print*,'lx,ly',lx,ly    
           do a1=0,2
             do sigma=0,1
              
               if(a1.le.1)  comp=(lx)*4+(a1)*2+sigma+(ly)*6*(d/2) 
               if(a1.eq.2)  comp=4*(d/2)+sigma+2*(lx)+(ly)*6*(d/2) 
                  write(121,*)lx,ly,a1+1,sigma+1,comp+1
                vec(comp+1)=complex(0.0d0,0.0d0)
                
               do jx=0,(d/2)-1  
                 do jy=0,(d/2)-1 
                    !    if(a1.le.2)  i_2=(jx-1+1)*4+(a1-1)*2+sigma+(jy-1+1)*6*(d/2)
                    !    if(a1.eq.3)  i_2=4*(d/2)+sigma+2*(jx-1+1)+(jy-1+1)*6*(d/2) 
                   nx=(mod((lx+jx),(d/2)))
                   ny=(mod((ly+jy),(d/2)))                
                   if(a1.le.1)  i_3=(nx)*4+(a1)*2+sigma+(ny)*6*(d/2)
                   if(a1.eq.2)  i_3=4*(d/2)+sigma+2*(nx)+(ny)*6*(d/2) 
                   !print*,'i_3',nx,ny,i_3              
                   k_r=(kx*(jx)+ky*(jy))
                   exp_k_r=complex(cos(k_r),-sin(k_r)) 
                   vec(comp+1)=vec(comp+1)+exp_k_r*(h( (i_3+1), m)) 
                 if(m==1)  write(122,*) comp+1,lx,ly,jx,jy,nx,ny,i_3,m
                   ! print*,vec(comp)
                  enddo  !jy loop         
                enddo    !jx loop          	      	        
              enddo      !sigma loop
            enddo        ! a1 loop
          enddo          ! ly loop
        enddo            ! lx loop     

        

  
       norm=0.0d0 
       do p=1,d1
          norm=norm+vec(p)*conjg(vec(p)) 
       enddo  ! p loop 
			 
        !print*,m,e,k_x,k_y,norm             
        if(norm.gt.eps_)then
          
          if((bottom_band1(k) .eq.dummy_value))then 
            bottom_band1(k)=e
            
            do p=1,d1
              mat_1(k,p)=vec(p)
             
              
              if(mod(p,2).ne.0) mat_1_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_1_down(k,p)=vec(p)

              
            enddo
            

          elseif((bottom_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_))then
            bottom_band2(k)=e
            !print*,k,e,m,'flat_band'
            do p=1,d1
              
              mat_2(k,p)=vec(p)
          
              if(mod(p,2).ne.0) mat_2_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_2_down(k,p)=vec(p)

            enddo

           ! print*,k,e,kx,ky
          elseif((flat_band1(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
            .and. (abs(e-bottom_band2(k)).gt.eps_))then
            flat_band1(k)=e
            !print*,k,e,m,'flat_band'
            do p=1,d1
              
              mat_3(k,p)=vec(p)
              
              if(mod(p,2).ne.0) mat_3_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_3_down(k,p)=vec(p)

            enddo
            
          elseif((flat_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
            .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_))then
            flat_band2(k)=e
            !print*,k,e,m,'flat_band'
            do p=1,d1
              
              mat_4(k,p)=vec(p)
              
              if(mod(p,2).ne.0) mat_4_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_4_down(k,p)=vec(p)

            enddo
            
          elseif((top_band1(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
            .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_)&
            .and. (abs(e-flat_band2(k)).gt.eps_))then
            top_band1(k)=e
            !print*,k,e,m,'flat_band'
            do p=1,d1
              
              mat_5(k,p)=vec(p)
              write(80,*)k,p,mat_5(k,p)
              if(mod(p,2).ne.0) mat_5_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_5_down(k,p)=vec(p)

            enddo

          elseif((top_band2(k) .eq.dummy_value).and. (abs(e-bottom_band1(k)).gt.eps_)&
            .and. (abs(e-bottom_band2(k)).gt.eps_).and. (abs(e-flat_band1(k)).gt.eps_)&
            .and. (abs(e-flat_band2(k)).gt.eps_).and. (abs(e-top_band1(k)).gt.eps_))then
            top_band2(k)=e
            !print*,k,e,m,'flat_band'
            do p=1,d1
              
              mat_6(k,p)=vec(p)
              write(81,*)k,p,mat_6(k,p)
              if(mod(p,2).ne.0) mat_6_up(k,p)=vec(p)
              if(mod(p,2).eq.0) mat_6_down(k,p)=vec(p)

            enddo
    
          endif        
       endif              
      enddo           
    enddo         
  enddo  

  
  do k_x=1,d/2
    kx =(2*pi/(d/2))*k_x
    do k_y=1,d/2

      ky=(2*pi/(d/2))*k_y
      k=k_x+(k_y-1)*(d/2)
      write(101,*)kx,ky,k,bottom_band1(k),bottom_band2(k),flat_band1(k),flat_band2(k),top_band1(k),top_band2(k)
    
    enddo
  enddo
  
 
end subroutine bloch_states

subroutine chern_number
  use global 
  use array
  implicit none
  integer::p,k,dkx_ky,count_k,kp1,p_,kyp1,kxp1,kp2,kp1p2
  complex*16::f12,u1,u2,u3,u4,u1k,u2k,u1k2,u2k1,u1k_,u2k_,u1k2_,u2k1_,f12_down
  complex*16::u_1,u_2,u_3,u_4,nk_up(d1),f1_2,nkp1_up(d1),nkp2_up(d1),nkp1p2_up(d1)
  double precision::k_r,kp1_r,kp2_r,kp1p2_r,norm_k,norm_kp1,norm_kp2,norm_kp1p2
  complex*16::nkp1_down(d1),nkp2_down(d1),nkp1p2_down(d1),nk_down(d1)
  double precision::norm_k_down,norm_kp1_down,norm_kp2_down,norm_kp1p2_down
  complex*16::u1_down,u2_down,u3_down,u4_down,c



  f12=complex(0.0d0,0.0d0)
  

  count_k=1

  do k_y=1,d/2
    ky =(2*pi/(d/2))*k_y
    do k_x=1,d/2
      kx=(2*pi/(d/2))*k_x
     
      kxp1=(mod(k_x,d/2))+1
      kyp1=(mod(k_y,d/2))+1

      kp1=kxp1+(k_y-1)*(d/2)    ! along kx direction
      kp2=k_x+(kyp1-1)*(d/2)    ! along ky direction
      kp1p2=kxp1+(kyp1-1)*(d/2) 
      
      
      
      norm_k=0.0d0
      norm_kp1=0.0d0
      norm_kp2=0.0d0
      norm_kp1p2=0.0d0
      
     

      do p=1,d1

        norm_k=norm_k+(mat_temp(count_k,p)*conjg(mat_temp(count_k,p)))
        norm_kp1=norm_kp1+(mat_temp(kp1,p)*conjg(mat_temp(kp1,p)))
        norm_kp2=norm_kp2+(mat_temp(kp2,p)*conjg(mat_temp(kp2,p)))
        norm_kp1p2=norm_kp1p2+(mat_temp(kp1p2,p)*conjg(mat_temp(kp1p2,p)))
        
        
      enddo
      
      
      do lx=1,d/2  
        do ly=1,d/2 
          !print*,'lx,ly',lx,ly
          k_r=kx*lx+ky*ly 
          kp1_r=  ((2*pi/(d/2))*kxp1)*lx + ((2*pi/(d/2))*k_y)*ly
          kp2_r=  ((2*pi/(d/2))*k_x)*lx + ((2*pi/(d/2))*kyp1)*ly
          kp1p2_r=  ((2*pi/(d/2))*kxp1)*lx + ((2*pi/(d/2))*kyp1)*ly


          do a1=1,3
            do sigma=1,2
             
              if(a1.le.2)  p_=(lx-1)*4+(a1-1)*2+sigma+(ly-1)*6*(d/2) 
              if(a1.eq.3)  p_=4*(d/2)+sigma+2*(lx-1)+(ly-1)*6*(d/2) 
              
              nk_up(p_)=complex(cos(k_r),sin(k_r))*(mat_temp(count_k,p_)/sqrt(norm_k))
              nkp1_up(p_)=complex(cos(kp1_r),sin(kp1_r))*(mat_temp(kp1,p_)/sqrt(norm_kp1))
              nkp2_up(p_)=complex(cos(kp2_r),sin(kp2_r))*(mat_temp(kp2,p_)/sqrt(norm_kp2))
              nkp1p2_up(p_)=complex(cos(kp1p2_r),sin(kp1p2_r))*(mat_temp(kp1p2,p_)/sqrt(norm_kp1p2))
              
             ! print*,p_,nkp1_up(p_)
            enddo
          enddo
        enddo
      enddo 
          
           
         
      u1=complex(0.0d0,0.0d0) 
      u2=complex(0.0d0,0.0d0) 
      u3=complex(0.0d0,0.0d0) 
      u4=complex(0.0d0,0.0d0) 


      do p=1,d1

        u1=u1+((conjg(nk_up(p))*nkp1_up(p)))     
        u2=u2+((conjg(nkp1_up(p))*nkp1p2_up(p)))       
        u3=u3+((conjg(nkp2_up(p))*nkp1p2_up(p)))      
        u4=u4+((conjg(nk_up(p))*nkp2_up(p)))

        
      enddo

      u1k=u1/sqrt(u1*conjg(u1))
      u2k1=u2/sqrt(u2*conjg(u2))
      u1k2=u3/sqrt(u3*conjg(u3))
      u2k=u4/sqrt(u4*conjg(u4))


      f12=f12+log(u1k*u2k1*(1/u1k2)*(1/u2k))
     
      count_k=count_k+1
      !write(110,*) kx,ky,real(f12/(2*pi*complex(0.0d0,1.0d0)))

    enddo
  enddo
  
  ch_no=f12/(2*pi*complex(0.0d0,1.0d0))
  
  
  
endsubroutine chern_number

subroutine check_states
  use global 
  use array
  implicit none
  integer::k
  complex*16::h3(d1,1),h1(d1,1),h2(1,d1),h4(1,1),h5(1,1)
  
  do k=1,(d**2)/4
    h1=complex(0.0d0,0.0d0)
    h2=complex(0.0d0,0.0d0)
    h3=complex(0.0d0,0.0d0)
    h4=complex(0.0d0,0.0d0)
    h5=complex(0.0d0,0.0d0)
    do i=1,d1
     
      h1(i,1)=mat_3_up(k,i)
       !if(k==6) write(99,*) mat_1(6,i),mat_2(6,i),mat_3(6,i)
    enddo
    
    do i=1,d1
     h2(1,i)=conjg(h1(i,1)) 
    enddo
   
    !h2=transpose(h1)
    h3=matmul(h_,h1)
    h4=matmul(h2,h3)
    h5=matmul(h2,h1)
    write(100,*)k,h4(1,1)/h5(1,1)
  enddo
endsubroutine check_states


subroutine chern_input
  use global
  use array
  implicit none

    if(i.eq.1) mat_temp=mat_1
    if(i.eq.2) mat_temp=mat_1_up      
    if(i.eq.3) mat_temp=mat_1_down     
    if(i.eq.4) mat_temp=mat_2
    if(i.eq.5) mat_temp=mat_2_up
    if(i.eq.6) mat_temp=mat_2_down
    if(i.eq.7) mat_temp=mat_3    
    if(i.eq.8) mat_temp=mat_3_up    
    if(i.eq.9) mat_temp=mat_3_down 
    if(i.eq.10) mat_temp=mat_4
    if(i.eq.11) mat_temp=mat_4_up     
    if(i.eq.12) mat_temp=mat_4_down
    if(i.eq.13) mat_temp=mat_5
    if(i.eq.14) mat_temp=mat_5_up    
    if(i.eq.15) mat_temp=mat_5_down   
    if(i.eq.16) mat_temp=mat_6  
    if(i.eq.17) mat_temp=mat_6_up    
    if(i.eq.18) mat_temp=mat_6_down


endsubroutine chern_input









