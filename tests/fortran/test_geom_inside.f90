! Exercise the real geom_tools module, without DRUtES main or a model solve.
program test_geom_inside
  use typy, only: rkind,qprec,ikind
  use globals, only: drutes_config
  use geom_tools, only: inside,inside_shoot
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(kind=rkind) :: triangle(3,2),reversed(3,2),p(2),offset(2),randoms(8),w(3),length
  real(kind=rkind) :: interval(2,1),polygon(4,2),nan
  real(kind=rkind), allocatable :: triangles(:,:,:)
  integer, allocatable :: seed(:)
  integer :: i,j,k,n,checks,unit,ne,np,expected,found,hits,permutation(3,6)
  logical :: boundary,answer,reference
  character(len=1024) :: path

  drutes_config%dimen=2
  checks=0
  call random_seed(size=n); allocate(seed(n)); seed=1729
  call random_seed(put=seed)
  permutation=reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[3,6])
  do i=1,5000
    call random_number(randoms)
    length=10.0_rkind**(-3+6*randoms(1))
    offset=0
    if(mod(i,2)==0) offset=[405591.710404_rkind,5352097.992453_rkind]
    triangle(1,:)=offset
    triangle(2,:)=offset+length*[0.5_rkind+randoms(2),randoms(3)-0.5_rkind]
    triangle(3,:)=offset+length*[randoms(4)-0.5_rkind,0.5_rkind+randoms(5)]
    ! Interior and exterior samples stay away from the roundoff boundary band.
    w=[0.1_rkind+0.3_rkind*randoms(6),0.1_rkind+0.3_rkind*randoms(7),0.0_rkind]
    w(3)=1-w(1)-w(2)
    do j=1,2
      if(j==2) w=[-0.1_rkind,0.4_rkind,0.7_rkind]
      p=triangle(1,:)+w(2)*(triangle(2,:)-triangle(1,:))+w(3)*(triangle(3,:)-triangle(1,:))
      reference=halfplane_reference(triangle,p)
      call require(reference.eqv.(j==1),'independent random-point oracle')
      do k=1,6
        reversed=triangle(permutation(:,k),:)
        answer=inside(reversed,p,atboundary=boundary)
        call require(answer.eqv.reference,'random triangle / all vertex orders')
        call require(.not.boundary,'random point not on boundary')
      end do
    end do
    p=triangle(1,:)
    call require(inside(triangle,p,atboundary=boundary),'vertex included')
    call require(boundary,'vertex marked boundary')
    p=triangle(1,:)+0.5_rkind*(triangle(2,:)-triangle(1,:))
    call require(inside(triangle,p,atboundary=boundary),'edge included')
    call require(boundary,'edge marked boundary')
  end do
  triangle(1,:)=[0.0_rkind,0.0_rkind]
  triangle(2,:)=[1.0_rkind,0.0_rkind]
  triangle(3,:)=[0.0_rkind,1.0_rkind]
  p=[0.2_rkind,0.3_rkind]
  do i=1,100
    call require(inside(triangle,p),'repeat deterministic interior')
    call require(.not.inside(triangle,[0.7_rkind,0.7_rkind]),'repeat deterministic exterior')
  end do
  p=[0.5_rkind,0.5_rkind]
  call require(inside(triangle,p,atboundary=boundary),'shared edge first triangle')
  call require(boundary,'shared edge first boundary flag')
  reversed(1,:)=[1.0_rkind,1.0_rkind]
  reversed(2,:)=triangle(3,:); reversed(3,:)=triangle(2,:)
  call require(inside(reversed,p,atboundary=boundary),'shared edge second triangle')
  call require(boundary,'shared edge second boundary flag')
  offset=[405591.710404_rkind,5352097.992453_rkind]
  do i=1,3
    reversed(i,:)=triangle(i,:)+offset
  end do
  call require(inside(reversed,offset+[0.5_rkind,0.5_rkind],atboundary=boundary),'UTM shared edge')
  call require(boundary,'UTM shared edge boundary flag')
  call require(.not.inside(reversed,offset+[0.5_rkind,0.500001_rkind]),'UTM just outside sloping edge')
  call require(inside(reversed,offset+[0.5_rkind,0.499999_rkind],atboundary=boundary),'UTM just inside edge')
  call require(.not.boundary,'UTM interior not boundary')
  p=[0.2_rkind,0.3_rkind]
  call require(inside(triangle*1.e-100_rkind,p*1.e-100_rkind),'very small scale')
  call require(inside(triangle*1.e100_rkind,p*1.e100_rkind),'very large scale')
  triangle(3,:)=[2.0_rkind,0.0_rkind]
  call require(.not.inside(triangle,[0.5_rkind,0.0_rkind],atboundary=boundary),'collinear rejected')
  call require(.not.boundary,'degenerate not marked boundary')
  triangle=0
  call require(.not.inside(triangle,[0.0_rkind,0.0_rkind]),'collapsed rejected')
  nan=ieee_value(0.0_rkind,ieee_quiet_nan)
  triangle(3,2)=nan
  call require(.not.inside(triangle,p),'nonfinite triangle rejected')
  triangle(1,:)=[0.0_rkind,0.0_rkind]
  triangle(2,:)=[1.0_rkind,0.0_rkind]
  triangle(3,:)=[0.0_rkind,1.0_rkind]
  call require(.not.inside(triangle,[nan,0.0_rkind]),'nonfinite point rejected')
  ! Explicit dimension overrides the ambient one, and optional arguments forward.
  drutes_config%dimen=1
  call require(inside(triangle,p,dimen_input=2_ikind),'explicit 2D dimension')
  interval(:,1)=[0.0_rkind,1.0_rkind]
  call require(inside(interval,[0.5_rkind]),'legacy interval interior')
  call require(.not.inside(interval,[1.5_rkind]),'legacy interval exterior')
  call require(inside_shoot(interval,[0.0_rkind],boundary,1_ikind),'legacy name accessible')
  call require(boundary,'legacy interval endpoint')
  drutes_config%dimen=2
  polygon(1,:)=[0.0_rkind,0.0_rkind]; polygon(2,:)=[1.0_rkind,0.0_rkind]
  polygon(3,:)=[1.0_rkind,1.0_rkind]; polygon(4,:)=[0.0_rkind,1.0_rkind]
  call require(inside(polygon,polygon(1,:),atboundary=boundary),'legacy polygon fallback')
  call require(boundary,'legacy polygon vertex')

  call get_command_argument(1,path)
  if(len_trim(path)>0) then
    open(newunit=unit,file=trim(path),status='old',action='read')
    read(unit,*) ne,np
    allocate(triangles(3,2,ne))
    do i=1,ne
      read(unit,*) ((triangles(j,k,i),k=1,2),j=1,3)
    end do
    do i=1,np
      read(unit,*) expected,p
      found=0; hits=0
      do j=1,ne
        if(inside(triangles(:,:,j),p)) then
          if(found==0) found=j
          hits=hits+1
        end if
      end do
      if(found/=expected .or. hits/=1) then
        print *, 'FAIL mesh query/expected/first/hits/point:',i,expected,found,hits,p
        error stop 'wrong containing triangle'
      end if
      call require(halfplane_reference(triangles(:,:,found),p),'mesh high-precision oracle')
    end do
    close(unit)
    print *, 'Rhine unique triangle queries passed:',np,' triangles:',ne
  end if
  print *, 'Deterministic inside checks passed:',checks
contains
  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if(.not.ok) then
      print *, 'FAIL ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine

  ! Independent higher-precision oracle: oriented edge half-planes, no division
  ! and no production barycentric helper or tolerance is reused.
  logical function halfplane_reference(t,point) result(ok)
    real(kind=rkind), intent(in) :: t(3,2),point(2)
    real(kind=qprec) :: a(2),b(2),r(2),cross(3)
    integer :: e,next
    do e=1,3
      next=mod(e,3)+1
      a=real(t(e,:),qprec); b=real(t(next,:),qprec); r=real(point,qprec)-a
      cross(e)=(b(1)-a(1))*r(2)-(b(2)-a(2))*r(1)
    end do
    ok=all(cross>0) .or. all(cross<0)
  end function
end program
