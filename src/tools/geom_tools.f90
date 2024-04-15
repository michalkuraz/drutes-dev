! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file geom_tools.f90
!! \brief Geometrical tools mainly for 2D.
!<




module geom_tools

  public :: init_transform
  public :: inside
  public :: inside3D
  public :: triarea
  public :: tetravol
  public :: dist
  public :: angle
  public :: inline
  public :: solve_bisect
  public :: aboveline
  public :: najdi_bod
  public :: max_height
  public :: interpolate
  public :: get_nz
  public :: find_neighbours
  public :: getcoor
  private :: shoot
  public :: map1d2d,map1d2dJ !later modified by J due to laziness
  public :: getnormal
  public :: get_layer
  public :: isboundary
  public :: get_normals3D

  
  contains
  
  function isboundary(quadpnt) result(answer)
    use typy
    use global_objs
    use globals
    
    type(integpnt_str), intent(in) :: quadpnt
    logical :: answer
    
    select case(quadpnt%type_pnt)
      case("obpt")
        if (observation_array(quadpnt%order)%boundary) then
          answer = .true.
        else
          answer = .false.
        end if
      case("ndpt")
        if (nodes%boundary(quadpnt%order)) then
          answer = .true.
        else
          answer = .false.
        end if
      case default
        answer = .false.
    end select
    
    answer = .false.
    
  end function isboundary 
  
  function get_layer(quadpnt) result(layer)
    use typy
    use globals
    use global_objs
    
    type(integpnt_str), intent(in) :: quadpnt
    integer(kind=ikind) :: layer
    integer(kind=ikind) :: el
    
    select case(quadpnt%type_pnt)
      case("ndpt")
        el = nodes%element(quadpnt%order)%data(1)
        layer=elements%material(el)
      case("obpt", "gqnd","xypt", "numb")
        el = quadpnt%element
        if (.not. (el >= 1 .and. el<=elements%kolik)) then
          print *, "error in quadpnt data", el, quadpnt%type_pnt
          print *, "called from geom_tools::get_layer"
          print *, "contact Michal -> michalkuraz@gmail.com"
          ERROR STOP
        end if
        layer = elements%material(el)
      case default
          print *, "err or in quadpnt data"
          print *, "called from geom_tools::get_layer"
          print *, "contact Michal -> michalkuraz@gmail.com"
          ERROR STOP
    end select
  
  end function get_layer
  
  
  subroutine getcoor(quadpnt, array)
    use typy
    use pde_objs
    use globals
    use debug_tools
    
    type(integpnt_str), intent(in) :: quadpnt
    real(kind=rkind), dimension(:), intent(out) :: array
    integer(kind=ikind), dimension(3) :: nd
    integer(kind=ikind) :: el, i
    real(kind=rkind) :: dist1, dist2, r1, r2, k1, k2, q1, q2
    real(kind=rkind), dimension(2) :: r,k,q,d
    real(kind=rkind), dimension(3,2) :: triagpnt
    real(kind=rkind), dimension(2,2) :: midpnt
    logical, dimension(2) :: vert 
     
    
    select case(quadpnt%type_pnt)
      case("gqnd")
      select case (drutes_config%dimen)
        case(1)
          el = quadpnt%element
          nd(1:2) = elements%data(el,1:2)
          dist1 = nodes%data(nd(2),1) - nodes%data(nd(1),1)
          array(1) = nodes%data(nd(1),1) + dist1*gauss_points%point(quadpnt%order,1)
        case(2)
         
          ! triag scheme	  
          ! 	   triag(3)
          !          |\ 
          ! 	       | \ 
          ! 	       |  \
          ! 	       |   \
          ! 	       |    \
          ! 	       |     \
          ! 	       |      \
          ! 	       |       \
          ! midpnt(1)|...x    \
          ! 	       |   .     \
          !      d(1)|   .      \
          ! 	       |   . d(2)  \
          ! triag(1) -------------\ triag(2)
          ! 	       midpnt(2)
          
        
          vert = .false.
          el = quadpnt%element
          do i=1,3
            triagpnt(i,:) = nodes%data(elements%data(el,i),:)
          end do
          d(1) = dist(triagpnt(1,:), triagpnt(3,:))
          d(2) = dist(triagpnt(1,:), triagpnt(2,:))
          !z distance on unit triangle
          r(1) = gauss_points%point(quadpnt%order,2)
          !x distance on unit triangle
          r(2) = gauss_points%point(quadpnt%order,1)
          
          midpnt(1,:) = triagpnt(1,:) + (triagpnt(3,:) - triagpnt(1,:)) * r(1)
          
          midpnt(2,:) = triagpnt(1,:) + (triagpnt(2,:) - triagpnt(1,:)) * r(2)
          
          if (abs(midpnt(2,1)-triagpnt(1,1)) > 100*epsilon(k1) ) then
            k(1) = (midpnt(2,2) - triagpnt(1,2))/(midpnt(2,1)-triagpnt(1,1))
          else
            vert(1) = .true.
          end if
          
          if (abs(midpnt(1,1)-triagpnt(1,1)) > 100*epsilon(k1) ) then
            k(2) = (midpnt(1,2) - triagpnt(1,2))/(midpnt(1,1)-triagpnt(1,1))
          else
            vert(2) = .true.
          end if
          
          where (.not.(vert))
            q(:) = midpnt(:,2) - k(:)*midpnt(:,1)
          end where
          
    ! 	    print *, "----beg-------"
    ! 	    print *, midpnt(1,:)
    ! 	    print *, midpnt(2,:)
    ! 	    print *, k
    ! 	    print *, vert
    ! 	    print *, "-----end-----"
          
          
          if (.not. vert(1) .and. .not. vert(2)) then
            array(1) = (q(1)- q(2))/(k(2)-k(1))
            if (abs(k(1) - k(2)) <= epsilon(k1)) then
              print *, "there must be a bug, both lines have equal tangent, or your mesh is really bad"
              print *, "interrupted from geom_tools::getcoor"
              ERROR STOP
            end if
            array(2) = k(1)*array(1) + q(1)
          else
            if (vert(1) .and. .not. (vert(2))) then
              array(1) = midpnt(1,1)
              array(2) = midpnt(2,2) + k(2)*(midpnt(1,1)-midpnt(2,1))
            else if (.not.vert(1) .and. (vert(2))) then
              array(1) = midpnt(2,1)
              array(2) = midpnt(1,2) + k(1)*(midpnt(2,1) - midpnt(1,1))
            else if ( vert(1) .and. vert(2) ) then
              print *, "there must be a bug, both lines are vertical"
              print *, "interrupted from geom_tools::getcoor"
              ERROR STOP
            end if
          end if   
        end select
      case("obpt")
        array = observation_array(quadpnt%order)%xyz
      case("ndpt")
        array = nodes%data(quadpnt%order,:)
    end select
   end subroutine getcoor

  !> returns the solution values on gaussian nodes positions
  subroutine interpolate(vals, points)
    use typy
    use globals
    use global_objs
    use globals2D
    use core_tools

    real(kind=rkind), dimension(:), intent(in) :: vals
    real(kind=rkind), dimension(:), intent(out) :: points
    integer(kind=ikind) :: i 
    real(kind=rkind) :: xder, yder

    select case(drutes_config%dimen)
      case(1)
        do i=1, ubound(gauss_points%weight,1)
          points(i) = (vals(2) - vals(1))*gauss_points%point(i,1) + vals(1) 
        end do
      case(2)
        call plane_derivative((/0.0_rkind, 0.0_rkind, vals(1)/), (/1.0_rkind, 0.0_rkind, vals(2)/), (/0.0_rkind, &
                               1.0_rkind, vals(3)/), xder, yder)

        do i=1, ubound(gauss_points%weight,1)
          points(i) = vals(1) + xder*gauss_points%point(i,1) + yder*gauss_points%point(i,2)
        end do

    end select

  end subroutine interpolate



  !> procedure to create transformation matrix 
  !! \f[  \left( \begin{array}{c}  e \\ f \end{array} \right) + \left( \begin{array}{cc}   a & b \\ c & d  \end{array}\right) \left( \begin{array}{c}  x \\ y  \end{array} \right)  = \left( \begin{array}{c}  x^T \\ y^T \end{array} \right) \f]
  !! where (e,f) is the vector of move and the matrix \f[ \left( \begin{array}{cc}   a & b \\ c & d  \end{array}\right) \f] is the matrix of rotation
  !! each line represents the element id, the first column in row represents numbers e,f,a,b the second column reprerents numbers f,c,d
  !<
   subroutine init_transform(aT,bT,cT,matrix)
    use typy
    !> nodes of triangle
    real(kind=rkind), dimension(:), intent(in) :: aT
    real(kind=rkind), dimension(:), intent(in) :: bT
    real(kind=rkind), dimension(:), intent(in) :: cT
    !> to store the transformation matrix
    real(kind=rkind), dimension(:,:), intent(out) :: matrix
    real(kind=rkind) :: a,b,c,d,e,f
    integer(kind=ikind) :: i
    
    

    !shape of the unit triangle

! [0,1]
! |\
! |  \
! |     \
! |       \
! |         \
! |           \
! |               \ 
! |----7|---8|---9|--\ 
! |     |    |    |   \
! |----4|---5|---6|---  \
! |     |    |    |       \
! |----1|---2|---3|-----   -\       
! |     |    |    |           \ 
! |     |    |    |             \
! -----------------------------------------
! [0,0]                          [1,0]

    e = aT(1)
    f = aT(2)

    a = bT(1) - e
    c = bT(2) - f

    b = cT(1) - e 
    d = cT(2) - f 
      
    matrix(1,1) = e
    matrix(1,2) = f
    matrix(2,1) = a
    matrix(2,2) = c
    matrix(3,1) = b
    matrix(3,2) = d

  end subroutine  init_transform


  !> function that provides z coordinate of inner boundary normal vector, the boundary is defined by two points: bod1, bod2
  function get_nz(el_id, order1, order2) result(zcoord)
    use typy
    use globals
    use global_objs

    integer(kind=ikind), intent(in) :: el_id, order1, order2
    real(kind=rkind), dimension(2) :: bod1, bod2
    real(kind=rkind) :: zcoord

    real(kind=rkind) :: tg
    integer(kind=ikind) :: order3


    bod1 = nodes%data(elements%data(el_id,order1),:)
    bod2 = nodes%data(elements%data(el_id,order2),:)

    if (abs(bod1(1)-bod2(1)) > 100*epsilon(tg)) then
      tg = (abs(bod2(2)-bod1(2))/abs(bod1(1)-bod2(1)))
      zcoord = 1.0_rkind/sqrt(1+tg*tg)
    else
      zcoord = 0
    end if

    select case(order1+order2)
      case(3)
        order3=3
      case(4)
        order3=2
      case(5)
        order3=1
    end select

    if (.not. aboveline(bod1, bod2, nodes%data(elements%data(el_id, order3),:))) then
      zcoord = -zcoord
    end if

  end function get_nz
  
  !> returns normal vector coordinates of an inner normal boundary vector to tetrahedron surface
  !! plane defined by given tetrahedron element wall is defined by
  !! \f[ ax +by +cz + d = 0 \f]
  !! function plane_coeff returns vector 
  !! \f[ coeffs = (a,b,c,d)^T \f]
  !! line perpendicular to the surface passing through the 4th tetrahedron point with coordinates \f[ (x_4, y_4, z_4)^T \f]
  !! is given by
  !! \f[ \begin{split} x &= x_4 + at \\ y &= y_4 + bt \\ z &= z_4 + ct \end{split} \f]
  !! \f[ t \f] is obtained from
  !! \f[ a(x_4 + at) + b(y_4 + bt) + c(z_4 + ct) + d = 0 \f]
  !! perpendicular surface projection of point 4 to surface has z coordinate
  !! \f[ z_i = z_4 + ct \f]
  !! resulting sin of z-coordinate is obtained from
  !! \f sinz =  \frac{||c||}{\sqrt{a^2 + b^2 + c^2}}
  !! if \f[ z_4 > z_i \f] inner vector points upwards else \f[ sinz = -sinz \f]
  !>
  function get_normals3D(a,b,c, exter) result(sins)
    use typy
    use debug_tools
    use core_tools
    
    
    real(kind=rkind), dimension(:), intent(in) :: a,b,c, exter
    real(kind=rkind), dimension(3) :: sins
    
    real(kind=rkind), dimension(4) :: coeffs
    real(kind=rkind), dimension(3,3) :: pts
    real(kind=rkind) :: length, t
    real(kind=rkind), dimension(3) :: inters
    integer(kind=ikind) :: i
    
    pts(1,:) = a
    pts(2,:) = b
    pts(3,:) = c
     
    
    call plane_coeff(pts, coeffs)
    
    length = 0
    
    do i=1,3
      length = length + coeffs(i)*coeffs(i)
    end do
    length = sqrt(length)
    
    if ( length > 100*epsilon(length)) then
      t = (-coeffs(4) - coeffs(1)*exter(1) - coeffs(2)*exter(2) - coeffs(3)*exter(3))/(length*length)
      inters = exter + coeffs(1:3)*t
      sins = abs(coeffs(1:3))/length
      if (exter(3) < inters(3)) then
        sins(3) = -sins(3)
      end if
      if (exter(2) < inters(2) ) then
        sins(2) = -sins(2)
      end if
      if (exter(1) < inters(1)) then
        sins(1) = -sins(1)
      end if
    else
      sins = 0
    end if
      
    
    

  
  end function get_normals3D
  
  
  function get_normals2D(a,b,exter) result(nvect)
  	use typy
  	real(kind=rkind), dimension(:), intent(in):: a,b,exter
  	real(kind=rkind), dimension(2) :: nvect
  	
  	real(kind=rkind) :: k
  	
  	if (abs(a(1)-b(1)) > 100*epsilon(a(1))) then
  	  k = (a(2)-b(2))/(a(1) - b(1))
  	  nvect = [k,1.0_rkind]/sqrt(k*k+1)
  	  if (.not. aboveline(a,b,exter)) then
  	  	nvect = -nvect
  	  end if
    else
      nvect = [0.0_rkind,1.0_rkind]
      if (aboveline(a,b,exter)) then
      	nvect = -nvect
      end if
  	end if
  	
  
  end function get_normals2D
  
  subroutine get_normals(el_id, bcpts, nvect) 
    use typy
    use pde_objs
    integer(kind=ikind), intent(in) :: el_id
    type(bcpts_str), intent(in) :: bcpts
    real(kind=rkind), dimension(:), allocatable, intent(out) :: nvect
    
    integer(kind=ikind) :: D
    
    D = drutes_config%dimen
    
    if (.not. allocated(nvect)) allocate(nvect(D))
    
    if (ubound(nvect,1) /= D) then
      print *, "runtime error, contact developer michalkuraz@gmail.com"
      print *, "exited from geom_tools::getnormals"
      ERROR STOP
    end if
    
    select case(D)
      case(1)
        if (nodes%data(elements%data(el_id, bcpts%border(1)),1) > nodes%data(elements%data(el_id, bcpts%extpt),1 ) ) then
          nvect(1) = -1
        else
          nvect(1) = 1
        end if
      case(2)
        nvect = get_normals2D(nodes%data(elements%data(el_id,bcpts%border(1)),:), &
                										  nodes%data(elements%data(el_id,bcpts%border(2)),:), &
                										  nodes%data(elements%data(el_id,bcpts%extpt),:) )
                                      
      case(3)
        nvect = get_normals3D(nodes%data(elements%data(el_id,bcpts%border(1)),:), &
                      nodes%data(elements%data(el_id,bcpts%border(2)),:), &
                      nodes%data(elements%data(el_id,bcpts%border(3)),:), &
                      nodes%data(elements%data(el_id,bcpts%extpt),:) )
    end select
    
    
    
  end subroutine get_normals
  
  function get_fluxsgn(el_id, bcpts, flux) result(sgn)
    use typy
    use pde_objs
    use globals
    use core_tools
    
    integer(kind=ikind), intent(in) :: el_id
    type(bcpts_str), intent(in) :: bcpts
    real(kind=rkind), dimension(:), intent(in) :: flux
    
    integer(kind=ikind) :: sgn
    
    integer(kind=ikind) :: D
    real(kind=rkind), dimension(3) :: line
    real(kind=rkind), dimension(4) :: plane
    real(kind=rkind), dimension(:,:), allocatable, save :: pts
    real(kind=rkind), dimension(:), allocatable, save :: ptsext
    real(kind=rkind) :: tinter
    
    
    D = drutes_config%dimen
    
    if (D > 1 .and. .not.(allocated(pts))) then
      allocate(pts(D,D))
      allocate(ptsext(D))
    end if
    
    select case(D)
      case(1)
        if (nodes%data(elements%data(el_id, bcpts%border(1)),1) > nodes%data(elements%data(el_id, bcpts%extpt),1 ) ) then
          if (flux(1) >= 0) then
            sgn = -1
          else
            sgn = 1
          end if
        else
         if (flux(1) < 0) then
            sgn = -1
          else
            sgn = 1
          end if
        end if
          
        
      case(2)
        if (norm2(flux) > 1e3*epsilon(1.0_rkind)) then
          pts(1,:) = nodes%data(elements%data(el_id, bcpts%border(1)),:)
          pts(2,:) = nodes%data(elements%data(el_id, bcpts%border(2)),:)
          
          ptsext = nodes%data(elements%data(el_id, bcpts%extpt),:) 
          
          call line_coeff(pts(1,:), pts(2,:), line)
          
          if ( abs(flux(1)*line(1) + flux(2)*line(2)) > 1e3*epsilon(1.0_rkind)) then
            tinter = -(line(3) + line(1)*ptsext(1) + line(2)*ptsext(2))/(flux(1)*line(1) + flux(2)*line(2))
          else
            sgn = 0
            return
          end if
          if (tinter > 0 ) then
            sgn = -1
          else
            sgn = 1
          end if
        else
          sgn = 0
        end if
          
        
      case(3)
      	
        if (norm2(flux) > 1e3*epsilon(1.0_rkind)) then
                  
          pts(1,:) = nodes%data(elements%data(el_id, bcpts%border(1)),:)
          pts(2,:) = nodes%data(elements%data(el_id, bcpts%border(2)),:)
          pts(3,:) = nodes%data(elements%data(el_id, bcpts%border(3)),:)
          
          ptsext = nodes%data(elements%data(el_id, bcpts%extpt),:) 
                    
          call plane_coeff(pts, plane) 
          
          if (abs(plane(1)*flux(1) + plane(2)*flux(2) + plane(3)*flux(3)) > 1e3*epsilon(1.0_rkind)) then
            tinter = -(plane(1)*ptsext(1) + plane(2)*ptsext(2) + plane(3)*ptsext(3) + plane(4)) / &
                        (plane(1)*flux(1) + plane(2)*flux(2) + plane(3)*flux(3))
          else
            sgn = 0
            return
          end if
    
          if (tinter > 0 ) then
            sgn = -1
          else
            sgn = 1
          end if
        else
          sgn = 0
        end if
    end select
  
  end function get_fluxsgn  
  
  subroutine line_coeff(A,B, coeffs)
    use typy
    
    real(kind=rkind), dimension(:), intent(in) :: A,B
    real(kind=rkind), dimension(:), intent(out) :: coeffs
    
    real(kind=rkind) :: k
    
    if (abs(A(1) - B(1)) > 1e3*epsilon(k)) then
      k = (B(2) - A(2))/(B(1) - A(1))
    else
      coeffs = [1.0_rkind, 0.0_rkind, A(1)]
      return
    end if
    
    coeffs(3) = A(2) - k*A(1)
    coeffs(1) = k
    coeffs(2) = -1.0_rkind
    
  end subroutine line_coeff

  !> function that provides x coordinate of inner boundary normal vector, the boundary is defined by two points: bod1, bod2
  function get_nx(el_id, order1, order2) result(xcoord)
    use typy
    use globals
    use global_objs

    integer(kind=ikind), intent(in) :: el_id, order1, order2
    real(kind=rkind), dimension(2) :: bod1, bod2
    real(kind=rkind) :: xcoord

    real(kind=rkind) :: tg, zcoord, slope
    integer(kind=ikind) :: order3


    bod1 = nodes%data(elements%data(el_id,order1),:)
    bod2 = nodes%data(elements%data(el_id,order2),:)

    if (abs(bod1(1)-bod2(1)) > epsilon(tg)) then
      tg = (abs(bod2(2)-bod1(2))/abs(bod1(1)-bod2(1)))
      zcoord = 1/sqrt(1+tg*tg)
    else
      zcoord = 0
    end if

    select case(order1+order2)
      case(3)
        order3=3
      case(4)
        order3=2
      case(5)
        order3=1
    end select
    
    xcoord = sqrt(1-zcoord*zcoord)
    
    if (abs(bod1(1)-bod2(1)) < 10*epsilon(bod1(1)) ) then
      if (nodes%data(elements%data(el_id, order3),1) < bod1(1) ) then
        xcoord = -xcoord
      end if
    else
      slope = (bod2(2) - bod1(2))/(bod2(1) - bod1(1))
      if (slope < 0) then
        if (.not. aboveline(bod1, bod2, nodes%data(elements%data(el_id, order3),:))) then
          xcoord = -xcoord
        end if
      else
        if ( aboveline(bod1, bod2, nodes%data(elements%data(el_id, order3),:))) then
          xcoord = -xcoord
        end if
      end if
    end if

  end function get_nx

 function inside(domain,bod, atboundary, dimen_input) result(true)
    use typy
    use globals
    use global_objs
    use debug_tools

    
    real(kind=rkind), dimension(:), intent(in) :: bod
    real(kind=rkind), dimension(2) :: a,b
    real(kind=rkind), dimension(:,:,:), allocatable :: inter
    real(kind=rkind), dimension(:,:), intent(in)  :: domain
    integer(kind=ikind) :: zasah, i, l, p, j, k
    integer(kind=ikind), dimension(:), allocatable :: checked
    logical, dimension(:), allocatable :: valid
    logical :: true
    logical, intent(out), optional :: atboundary
    integer(kind=ikind), intent(in), optional :: dimen_input
    real(kind=rkind), dimension(:), allocatable :: smer
    integer(kind=ikind) :: dimen_loc
    
    
    
    if (present(dimen_input)) then
      dimen_loc = dimen_input
    else
      dimen_loc = drutes_config%dimen
    end if
    
    select case(dimen_loc)
      case(1)
        if ((domain(1,1) <= bod(1) .and. domain(2,1) >= bod(1)) .or. &
            (domain(1,1) >= bod(1) .and. domain(2,1) <= bod(1)) ) then
          true = .true.
          
          if (present(atboundary)) then
            if (abs( bod(1) - domain(1,1)) < epsilon(bod(1)) .or. abs( bod(1) - domain(2,1)) < epsilon(bod(1))) then
               atboundary = .true.
            else
               atboundary = .false.
            end if
          end if
        else
          true = .false.
        end if

      case(2)

        allocate(smer(ubound(domain,1)+1))
        allocate(inter(ubound(smer,1),ubound(domain,1),2))
        allocate(valid(ubound(domain,1)))
        allocate(checked(ubound(smer,1)))
        checked = 0


        do i=1, ubound(domain,1)
          if (i < ubound(domain,1)) then
            l = i
            p = i+1
          else
            l = i
            p = 1
          end if
	  
          a = domain(l,:)
          
          b = domain(p,:)
          

          if (inline(domain(l,:), domain(p,:), bod)) then
            true = .TRUE.
            if (present(atboundary)) atboundary = .TRUE.
            RETURN
          end if

          do j=1, ubound(smer,1)
            smer(j) = rand() + 0.25
            call shoot(bod, smer(j), domain,  inter(j,:,:), valid)
            zasah = 0
            do k=1, ubound(domain,1)
              if (valid(k)) then
                if (dist(a,inter(j,k,:), dimen_loc) > 10*epsilon(smer(1)) .and. dist(b,inter(j,k,:), dimen_loc) & 
                                                                                                      > 10*epsilon(smer(1))) then
                  zasah = zasah + 1
                else
                  checked(j) = -1
                end if
              end if
            end do
            
            if (checked(j) /= -1) then
              if (modulo(zasah,2) /= 0 .and. zasah /= ubound(domain,1)) then
                checked(j) = 1
              end if
            end if
  
	      
          end do
        end do
        j=0
        do i=1, ubound(checked,1) 
          if (checked(i) == 0) then
            true = .false.
            if (present(atboundary)) atboundary = .false.
            RETURN
          end if
          j = j + checked(i)
        end do
        
        if (j == -ubound(checked,1)) then
          print *, "bug in geom_tools::inside, contact developer with a full bug report"
          ERROR stop
        end if
        
        true = .true.
        if (present(atboundary)) atboundary = .false.
      
      case default
        ERROR stop "generated from geom_tools::inside"
    end select

  end function inside
  
  
  function inside3D(volume, bod) result(true)
    use typy
    use linalg
    use debug_tools
    use core_tools
    
    real(kind=rkind), dimension(:,:,:), intent(in) :: volume
    real(kind=rkind), dimension(:), intent(in) :: bod
    logical :: true

    real(kind=rkind), dimension(3) :: inter
    real(kind=rkind), dimension(4) :: coeff
    real(kind=rkind), dimension(3,2) :: domain2D
    logical :: fail, hit
    integer(kind=ikind) :: pt, surf, bang
    real(kind=rkind) :: t, xder, yder
    real(kind=rkind), dimension(2) :: A, B, C
    
    integer(kind=ikind) :: i
    
    if (ubound(volume,2) /= 3 ) then
      print *, "the volume must be surrounded by triangles, typicaly Delauney triangulation"
      print *, "exited from geom_tools::inside3D"
      ERROR STOP
    end if
    
    if (ubound(volume,3) /= 3 ) then
      print *, "points in 3D must have three coordinates"
      print *, "exited from geom_tools::inside3D"
      ERROR STOP
    end if
    

    bang = 0
    do surf=1, ubound(volume,1)

      call plane_coeff(volume(surf, :, :), coeff)
      
      
  
      inter(1) = bod(1)
      inter(3) = bod(3)
      if (abs(coeff(2)) > 1e3*epsilon(coeff(2))) then
        t = (coeff(1)*bod(1) + coeff(2)*bod(2) + coeff(3)*bod(3)  + coeff(4))/coeff(2)
        inter(2) = bod(2) + t
          
          if (inter(2) >= bod(2)) then
         
            do pt=1,3
              domain2D(pt,:) = [volume(surf, pt, 1), volume(surf, pt, 3)]
            end do
            if (inside(domain2D, [inter(1), inter(3)], dimen_input = 2_ikind)) bang = bang + 1
          end if
        else
          if (abs(inter(1) - volume(surf, 1, 1)) < epsilon(t) .and. abs(inter(1) - volume(surf, 2, 1)) < epsilon(t) .and. &
              abs(inter(1) - volume(surf, 3, 1)) < epsilon(t)) then
           
            do pt=1,3
              domain2D(pt,:) = [volume(surf, pt, 2), volume(surf, pt, 3)]
            end do 
         
            if (inside(domain2D, [inter(2), inter(3)], dimen_input = 2_ikind)) then
              true = .true.
              RETURN
            end if
          end if
        end if
      end do
            
        
        
        
        
        
        
          
    
    if (modulo(bang, 2_ikind) /= 0) then
      true = .true.
    else
      true = .false.
    end if
    
  
  end function inside3D
  
  subroutine shoot(bod, k, domain, inter, valid)
    use typy
    real(kind=rkind), dimension(:), intent(in) :: bod
    real(kind=rkind), intent(in) :: k
    real(kind=rkind), dimension(:,:), intent(in) :: domain
    real(kind=rkind), dimension(:,:), intent(out) :: inter
    logical, dimension(:), intent(out) :: valid
    
    real(kind=rkind) :: q, qq, kk
    real(kind=rkind), dimension(2) :: a,b
    integer(kind=ikind) :: l,p, i
    
    q = bod(2) - k*bod(1)
    
    do i=1, ubound(domain,1)
      if (i < ubound(domain,1)) then
        l = i
        p = i+1
      else
        l = i
        p = 1
      end if 
      
      a = domain(l,:)
    
      b = domain(p,:)
      
      if (abs(a(1)-b(1)) > 10*epsilon(a(1))) then
        kk = (b(2) - a(2))/(b(1)-a(1))
        qq = b(2) - b(1)*kk
        if (abs(k-kk) > 10*epsilon(k)) then
          inter(i,1) = (qq-q)/(k-kk)
          inter(i,2) = k*inter(i,1) + q
        else
          if (abs(q-qq) > 10*epsilon(q)) then
            inter(i,:) = (/huge(k), huge(k)/)
          else
            if (a(1) > b(1)) then
              inter(i,:) = a
            else
              inter(i,:) = b
            end if
          end if
        end if
      else
        inter(i,1) = a(1)
        inter(i,2) = k*inter(i,1) + q
      end if
      if (inline(a,b, inter(i,:)) .and. inter(i,1) > bod(1)) then
        valid(i) = .true.
      else
        valid(i) = .false.
      end if
    end do
    
 
	
  end subroutine shoot


  function triarea(a,b,c) result(area)
    use typy
    !> triangle points
    real(kind=rkind), dimension(:), intent(in) :: a,b,c
    !> result
    real(kind=rkind) :: area
    real(kind=rkind) :: la,lb,lc
    
    lc = dist(a,b)

    lb = dist(a,c)
    
    la = dist(b,c)
    
    area = 0.25*sqrt((la+lb+lc)*(lb + lc - la)*(lc + la - lb)*(la + lb - lc))
  
  end function triarea
  
  function tetravol(a,b,c,d) result(vol)
    use typy
    use simplelinalg
    use core_tools
    
    real(kind=rkind), dimension(:), intent(in) :: a,b,c,d
    real(kind=rkind) :: vol
    
    real(kind=rkind), dimension(4,4) :: mtx
    integer(kind=ikind) :: i
    
  
    mtx(1,:) = [a, 1.0_rkind]
    mtx(2,:) = [b, 1.0_rkind]
    mtx(3,:) = [c, 1.0_rkind]
    mtx(4,:) = [d, 1.0_rkind]
    
    vol = 1.0_rkind/6.0_rkind*abs(determinant(mtx))
    
  end function tetravol

  
  function dist(A,B, dimen_input) result(l)
    use typy
    use globals
    
    real(kind=rkind), dimension(:), intent(in) :: A,B
    integer(kind=ikind), intent(in), optional :: dimen_input
    real(kind=rkind) :: l
    integer(kind=ikind) :: dimen_loc
    
    if (present(dimen_input)) then
      dimen_loc = dimen_input
    else
      dimen_loc = drutes_config%dimen
    end if
    
    select case(dimen_loc)
      case(1)
        l = abs(A(1) - B(1))
      case(2)
        l = sqrt((A(1)-B(1))*(A(1)-B(1)) + (A(2)-B(2))*(A(2)-B(2)))
      case(3)
        l = sqrt( (A(1)-B(1))*(A(1)-B(1)) + (A(2)-B(2))*(A(2)-B(2)) + (A(3)-B(3))*(A(3)-B(3)) )
    end select
    
  
  
  end function dist


  function angle(A, B, C) result(gamma)
    use typy
    
    real(kind=rkind), dimension(:), intent(in) :: A, B, C
    real(kind=rkind) :: gamma
    real(kind=rkind) :: la, lb, lc
    
    lc = dist(A,B)
    lb = dist(A,C)
    la = dist(B,C)
    
    if (la > 10*epsilon(la) .and. lb > 10*epsilon(la)) then
      gamma = acos(min((lc**2-la**2-lb**2)/(-2.0_rkind*la*lb),1.0_rkind))
    else
      gamma = acos(1.0_dprec)
    end if

  end function angle


    !> function to determine whether the third point is on line formed by the first two points
  function inline(A,B,bod) result(true)
    use typy
    real(kind=rkind), dimension(:), intent(in) :: bod, a,b
    real(kind=rkind) :: precision
    logical :: true
    

    precision = 1000*epsilon(precision)


    if (abs(dist(a,b, 2_ikind) - dist(a,bod, 2_ikind) - dist(b,bod, 2_ikind)) < precision) then
      true = .true.
    else
      true = .false.
    end if


    
  end function inline


  !> bisecton method searching for a root of function f on defined intervals xmin and xmax
  subroutine solve_bisect(xmin, xmax, f, reps, ierr, solution)
    use typy 
    use globals
    use core_tools

    !> a interval where the root is searched
    real(kind=rkind), intent(in) :: xmin, xmax
    !> function which is processed
    interface
      pure function f(x) result(y)
        use typy
        real(kind=rkind), intent(in) :: x
        real(kind=rkind) :: y
      end function
    end interface
    !> accuracy
    real(kind=rkind), intent(in) :: reps
    !> error code output
    integer, intent(out) :: ierr
    !> the result
    real(kind=rkind), intent(out) :: solution

    real(kind=rkind) :: xmid, x1, x2

    if ((f(xmin) > 0 .and. f(xmax) > 0) .or. (f(xmin) < 0 .and. f(xmax) < 0)) then
      ierr = -1
      call  write_log("WARNING! geom_tools::solve_bisect returned error code -1")
      RETURN
    end if

    x1=xmin
    x2=xmax
    xmid = (xmin+xmax)/2.0

    do 
      if ((f(x1) > 0 .and. f(x2) > 0) .or. (f(x1) < 0 .and. f(x2) < 0)) then
        ierr = -2
        call  write_log("WARNING! geom_tools::solve_bisect returned error code -2")
        RETURN
            end if
            if (f(x1) > 0 .and. f(x2) < 0) then
        if (f(xmid) > 0) then
          x1 = xmid
        else
          x2 = xmid
        end if
            else
        if (f(xmid) > 0) then
          x2 = xmid
        else
          x1 = xmid
        end if
            end if
            xmid = (x1+x2)/2.0
            if (abs(x1-x2) < reps) then
        solution = (x1+x2)/2.0
        ierr = 1
        EXIT
      end if
    end do

  end subroutine solve_bisect


  !> function to determine if the point C lies above or below the line defined by points A and B
  !! result is logical value nad, if true point lies above if false it lies below
  !! if the line is vertical, then if true the point lies on the left hand side from the line, if false the point lies on the right hand side from the line.
  !<
  function aboveline (a,b,c) result(nad)
    use typy
    use globals
    use globals2d
    real(kind=rkind), dimension(:), intent(in) :: a,b,c
    logical :: nad
    real(kind=rkind) :: xa1, yb1, q1, xa2, yb2, q2, x_int, y_int
    !first define the coefficient of the line equation defined by A and B
    
    if (abs(b(1)-a(1)) > 100*epsilon(xa1) .and. abs(b(2)-a(2)) > 100*epsilon(xa1)) then
      xa1 = (b(2) - a(2))/(b(1) - a(1))
      
      yb1 = -1.0_rkind
      
      q1 = -xa1*a(1) - yb1*a(2)

      xa2 = yb1
      
      yb2 = -xa1
      
      q2 = -xa2*c(1) - yb2*c(2) 
      
      x_int = (-q1 + yb1*q2/yb2)/(xa1 - yb1*xa2/yb2)
      
      y_int = (-q2 - xa2*x_int)/yb2
      
      
      if (y_int >= c(2)) then
        nad = .false.
      else
        nad = .true.
      end if
    else
      if (abs(b(1)-a(1)) < 100*epsilon(xa1)) then ! vertical line
 
        if (c(1) >= b(1)) then 
          nad = .false.
        else 
          nad = .true.
        end if
        
      end if
      
      if (abs(b(2)-a(2)) < 100*epsilon(xa1)) then ! horizontal line
        
        if (c(2) >= b(2)) then
          nad = .true.
        else
          nad = .false.
        end if
      end if
    end if
    
  end function aboveline
  
  
  !> this function returns outer normal boundary vector
  !! bcpoints - nodes ids at the element boundary, integer, dimension(2)
  !! thirdpt - coordinates of a point at at element, which lies out of the boundary. It is important in order to obtain proper positive direction of the normal vector
  !<    
  subroutine getnormal(bcpoints, thirdpt, nvect)
    use typy
    use globals
    
    real(kind=rkind), dimension(:,:), intent(in) :: bcpoints
    real(kind=rkind), dimension(:), intent(in) :: thirdpt
    real(kind=rkind), dimension(:), allocatable, intent(out) :: nvect
    
    real(kind=rkind) :: k
    
    if (drutes_config%dimen == 1) then
      print *, "RUNTIME ERROR, function geomtools::getnormal is called for 1D problem"
      print *, "if (desporate) contact Michal -> michalkuraz@gmail.com"
      ERROR STOP
    end if    
    
    if (.not. allocated(nvect)) allocate(nvect(drutes_config%dimen))
  
    if (abs(bcpoints(1,1)-bcpoints(2,1)) > 100*epsilon(bcpoints(1,1)) .and. &
      abs(bcpoints(1,2)-bcpoints(2,2)) > 100*epsilon(bcpoints(1,1))) then
     
      k = (bcpoints(1,2)-bcpoints(2,2))/(bcpoints(1,1)-bcpoints(2,1))
     
      if (aboveline(bcpoints(1,:), bcpoints(2,:), thirdpt)) then
        nvect(2) = -1
      else
        nvect(2) = 1
      end if
      
      ! change x for y
      if (aboveline((/bcpoints(1,2), bcpoints(1,1)/), (/bcpoints(2,2), bcpoints(2,1)/), (/thirdpt(2), thirdpt(1)/))) then
        nvect(1) =  -abs(1.0_rkind/k)
      else
        nvect(1) = abs(1.0_rkind/k)
      end if
      
      nvect = nvect/(sqrt(nvect(1)*nvect(1)+nvect(2)*nvect(2)))
     
    !vertical line
    else if (abs(bcpoints(1,1)-bcpoints(2,1)) < 100*epsilon(bcpoints(1,1)) .and. &
       abs(bcpoints(1,2)-bcpoints(2,2)) > 100*epsilon(bcpoints(1,1))) then
       
       if (aboveline(bcpoints(1,:), bcpoints(2,:), thirdpt)) then
         nvect(1) = 1
       else
         nvect(1) = -1
       end if
       
       nvect(2) = 0
       
     !horizontal line
    else if (abs(bcpoints(1,1)-bcpoints(2,1)) > 100*epsilon(bcpoints(1,1)) .and. &
       abs(bcpoints(1,2)-bcpoints(2,2)) < 100*epsilon(bcpoints(1,1))) then
       
       if (aboveline(bcpoints(1,:), bcpoints(2,:), thirdpt)) then
         nvect(2) = -1
       else
         nvect(2) = 1
       end if
       
       nvect(1) = 0
     else
       print *, "runtime error, exited from geomtools::getnormal seems like your boundary is defined by a single point"
       print *, "check your mesh, if your mesh is ok, then contact Michal -> michalkuraz@gmail.com"
       ERROR STOP
     end if
  
  end subroutine getnormal


 !> subroutine seeking for the point which is on line identified by points a and b, and line between this point and c is rectangular to the line <a,b>
  subroutine najdi_bod(a,b,c,bod)
    use typy
    !> element points
    real(kind=rkind), dimension(:), intent(in) :: a,b,c
    !> the rectangular pointer
    real(kind=rkind), dimension(:), intent(out) :: bod

    real(kind=rkind) :: k, q, q2


    if (abs(b(1) - a(1)) > 100*epsilon(a(1))) then
      k = (b(2)-a(2))/(b(1)-a(1))

      q = a(2) - k*a(1)

      q2 = -c(1) - k*c(2)

      bod(2) = (-k*q2+q)/(k*k+1)

      bod(1) = -k*bod(2)-q2

    else
      bod(1) = a(1)
      bod(2) = c(2)
    end if
      
  end subroutine najdi_bod

  !> subroutine searching for the maximum height on triangle element
  subroutine max_height(a,b,c, height)
    use typy
    !>input element edges
    real(kind=rkind), dimension(:), intent(in) :: a,b,c
    !> the maximal element height
    real(kind=rkind), intent(out) :: height
    
    real(kind=rkind), dimension(2) :: bod
    real(kind=rkind) :: ha,hb,hc
    real(kind=rkind) :: la, lb, lc
    
    
    call najdi_bod(a,b,c, bod)
    
    hc = dist(c,bod)
    
    lc = dist(a,b)
    la = dist(c,b)
    lb = dist(a,c)
    
    ha = la/lc*hc
    hb = lb/lc*hc
    
    height = max(ha,hb,hc)
!     height = (ha + hb + hc)/3
    
  
  end subroutine max_height


  function gravity_center(a,b,c) result(gc)
    use typy
    use linalg
    real(kind=rkind), dimension(:), intent(in) :: a,b,c
    real(kind=rkind), dimension(2) :: gc

    real(kind=rkind), dimension(2,2) :: line_pars, matrix
    real(kind=rkind), dimension(2) :: point, bside
    integer(kind=ikind) :: i
    logical :: vertical

    
    point = (b+c)/2.0_rkind
    


    if (abs(a(1)-point(1)) > 1000*epsilon(a(1))) then
      bside(1) = a(2)
      bside(2) = point(2)

      matrix(:,2) = 1.0_rkind
      matrix(1,1) = a(1)
      matrix(2,1) = point(1)

      call gem(matrix, bside, line_pars(1,:))

      vertical = .false.

    else
      gc(1) = (a(1)+point(1))/2.0_rkind
      vertical = .true.
      i = 2
    end if
      
    
    point = (a + c)/2.0_rkind
   
    if (abs(b(1)-point(1)) > 1000*epsilon(a(1))) then
      bside(1) = b(2)
      bside(2) = point(2)
      matrix(:,2) = 1
      matrix(1,1) = b(1)
      matrix(2,1) = point(1)
    
      call gem(matrix, bside, line_pars(2,:))

    else
      gc(1) = (b(1)+point(1))/2.0_rkind
      vertical = .true.
      i = 1
    end if
      
    if (.not.(vertical)) then
      bside = -line_pars(:,2)

      matrix(:,1) = line_pars(:,1)

      matrix(:,2) = -1.0_rkind

      call gem(matrix, bside, gc)
    else
      gc(2) = line_pars(i,1)*gc(1) + line_pars(i,2)
    end if
    

  end function gravity_center
  
  
    subroutine find_neighbours(el, nd)
      use typy
      use globals
      use global_objs
      use printtools
      use core_tools
      use debug_tools
      use debug_tools

      type(element), intent(in out) :: el
      type(node), intent(in) :: nd
      
      integer(kind=ikind) :: i,j,k,l,m,upward,downward,pos, n, o
      
      el%neighbours = 0_ikind
      
      select case(drutes_config%dimen)
        case(1)
          do i=1, el%kolik - 1
            el%neighbours(i,1) = i-1
            el%neighbours(i,2) = i+1
           end do
           el%neighbours(elements%kolik,1) = elements%kolik - 1
           el%neighbours(elements%kolik,2) = 0
        case(2:3)
            !set neighbours
            do i=1, el%kolik
              pos = 0
              j = i
              k = i

              call progressbar(int(100*i/el%kolik))

              okoli: do 
                j = min(j+1, el%kolik+1)
                k = max(k-1, 0_ikind)

                if (j <= el%kolik) then
                  upward = 0
                  moje1: do l=1,ubound(el%data,2)
                      nasel1: do m=1,ubound(el%data,2)
                        if (el%data(i,l) == el%data(j,m) .and. i/=j) then
                          upward = upward + 1
                          if (upward == drutes_config%dimen) then 
                      pos = pos + 1
                      el%neighbours(i,pos) = j
                      EXIT nasel1
                          end if
                        end if
                    end do nasel1
                  end do moje1
                end if

                if (k > 0) then
                  downward = 0
                  moje2: do l=1,ubound(el%data,2)
                    nasel2: do m=1,ubound(el%data,2)
                        if (el%data(i,l) == el%data(k,m) .and. i /= k) then
                          downward = downward + 1
                          if (downward == drutes_config%dimen) then 
                      pos = pos + 1
                      el%neighbours(i,pos) = k
                      EXIT nasel2
                          end if
                        end if
                    end do nasel2 
                  end do moje2
                end if

                if (pos == ubound(el%data,2) .or. (j == el%kolik+1 .and. k == 0_ikind)) then
                  EXIT okoli
                end if
              end do okoli
            end do
            
      end select
      

    end subroutine find_neighbours

  
    subroutine map1d2d(filename)
      use typy
      use global_objs
      use globals
      use pde_objs
      use core_tools
      use readtools
      use debug_tools
      
      character(len=*), intent(in) :: filename
      
      integer :: fileid, ierr
      integer(kind=ikind) :: i, counter, k, l, m, proc, j, ii
      real(kind=rkind), dimension(:,:), allocatable :: input
      real(kind=rkind) :: tmp, value, grad
      logical :: val_defined
      
      call find_unit(fileid)
      
      
      open(unit=fileid, file=cut(filename), action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open file with the vertical distribution of the initial condition, exiting...."
        print *, "called from map1d2d::drutes_init"
        ERROR STOP
      end if
      
      counter = 0
      
      do
        call comment(fileid)
        read(fileid, fmt=*, iostat=ierr) tmp
        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
      end do
      
      allocate(input(counter, ubound(pde,1)+1))
      
      close(fileid)
      
      open(unit=fileid, file=cut(filename), action="read")
      
      do i=1, counter
        call fileread(input(i,:), fileid)
      end do
            
      do i=1, elements%kolik
        do j=1, ubound(elements%data,2)
          k = elements%data(i,j)
          l = nodes%edge(k)
          do proc=1, ubound(pde,1)
            m = pde(proc)%permut(k)
            if (m == 0) then
              if (l>100) then
                call pde(proc)%bc(l)%value_fnc(pde(proc), i, j, value)
                val_defined=.true.
              else
                print *, "runtime error, exited from map1d2d::drutes_init"
                ERROR STOP
              end if             
            else
              val_defined=.false.
              do ii=1, ubound(input,1)-1
                if (nodes%data(k,2) >= input(ii,1) .and. nodes%data(k,2) <= input(ii+1,1)) then
                  grad = (input(ii+1,proc+1) - input(ii,proc+1))/(input(ii+1,1) - input(ii,1))
                  value =  input(ii,proc+1) + grad*(nodes%data(k,2) - input(ii,1))
                  val_defined=.true.
                  if (isnan(value)) then
                    call file_error(fileid, "bad input file with initial conditions")
                  end if
                  EXIT
                end if
              end do
              if (.not. val_defined) then
                if (nodes%data(k,2) < minval(input(:,1))) then
                  value = input(minloc(input(:,1),1), proc)
                  val_defined=.true.
                end if
                if (nodes%data(k,2) < maxval(input(:,1))) then
                  value = input(maxloc(input(:,1),1), proc)
                  val_defined=.true.
                end if
                if (.not. val_defined) then
                  call file_error(fileid, "Have you covered the entire vertical profile with data for the initial condition?")
                end if
              end if
            end if
            pde(proc)%solution(k) =  value 
          end do
        end do
      end do
            
      close(fileid)
 
      
    end subroutine map1d2d
    
    
    subroutine map1d2dJ(pde_loc,filename, correct_h)
      use typy
      use global_objs
      use globals
      use pde_objs
      use core_tools
      use readtools
      use debug_tools
      
      character(len=*), intent(in) :: filename
      class(pde_str), intent(in out) :: pde_loc
      logical, intent (in) :: correct_h
      integer :: fileid, ierr
      integer(kind=ikind) :: i, counter, k, l, m, proc, j, ii, D
      real(kind=rkind), dimension(:,:), allocatable :: input
      real(kind=rkind) :: tmp, value, grad
      logical :: val_defined
      
      
      D = drutes_config%dimen
      call find_unit(fileid)
      
      
      open(unit=fileid, file=cut(filename), action="read", iostat=ierr)
      if (ierr /= 0) then
        print *, "unable to open file with the vertical distribution of the initial condition, exiting...."
        print *, "called from map1d2dJ::drutes_init"
        ERROR STOP
      end if
      
      counter = 0
      
      do
        call comment(fileid)
        read(fileid, fmt=*, iostat=ierr) tmp
        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
      end do
      
      allocate(input(counter, 3))
      
      close(fileid)
      
      open(unit=fileid, file=cut(filename), action="read")
      
      do i=1, counter
        call fileread(input(i,:), fileid)
      end do

      do i=1, elements%kolik
        do j=1, ubound(elements%data,2)
          k = elements%data(i,j)
          l = nodes%edge(k)
            m = pde_loc%permut(k)
            if (m == 0) then
              if (l>100) then
                do ii=1, ubound(input,1)-1
    
                  if (nodes%data(k,D) >= input(ii,2) .and. nodes%data(k,D) <= input(ii+1,2)) then
                    grad = (input(ii+1,3) - input(ii,3))/(input(ii+1,2) - input(ii,2))
                    value =  input(ii,3) + grad*(nodes%data(k,D) - input(ii,2))
                    val_defined=.true.
                    if (isnan(value)) then
                      call file_error(fileid, "bad input file with initial conditions")
                    end if
                    EXIT
                   end if
                 end do
                 if(.not. val_defined) then
                   call pde(proc)%bc(l)%value_fnc(pde(proc), i, j, value)
                   val_defined=.true.
                 end if 
              else
                print *, "runtime error, exited from map1d2dJ::drutes_init"
                ERROR STOP
              end if             
            else
              val_defined=.false.
              do ii=1, ubound(input,1)-1
                if (nodes%data(k,D) >= input(ii,2) .and. nodes%data(k,D) <= input(ii+1,2)) then
                  
                  grad = (input(ii+1,3) - input(ii,3))/(input(ii+1,2) - input(ii,2))
                  value =  input(ii,3) + grad*(nodes%data(k,D) - input(ii,2))

                  val_defined=.true.
                  if (isnan(value)) then
                    call file_error(fileid, "bad input file with initial conditions")
                  end if
                  EXIT
                end if
              end do
              if (.not. val_defined) then
                if (nodes%data(k,D) < minval(input(:,2))) then
                  value = input(minloc(input(:,2),1),3)
                  val_defined=.true.
                end if
                if (nodes%data(k,D) > maxval(input(:,2))) then
                  value = input(maxloc(input(:,2),1), 3)
                  val_defined=.true.
                end if
                if (.not. val_defined) then
                  print *, nodes%data(k,D) ; stop
                  call file_error(fileid, "Have you covered the entire vertical profile with data for the initial condition?")
                end if
              end if
            end if
            pde_loc%solution(k) =  value 
            if(correct_h) then
              pde_loc%solution(k) =  pde_loc%solution(k)+nodes%data(k,D) 
            end if 
        end do
      end do
        
      close(fileid)
      
    end subroutine map1d2dJ
    

    

end module geom_tools
