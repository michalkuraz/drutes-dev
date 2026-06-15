module ncmap
  use typy
  use global_objs
  use globals
  use ncglobvars
  use ncmesh

  implicit none

  private

  public :: mapel
  public :: find_nc_element_for_point
  public :: fem_element_center

  contains


  subroutine mapel()



    integer(kind=ikind) :: el
    real(kind=rkind) :: xc, yc
    logical :: ok_center

    if (allocated(el2ncgrid)) then
      deallocate(el2ncgrid)
    end if

    allocate(el2ncgrid(elements%kolik))

    el2ncgrid(:) = -1_ikind

    do el = 1_ikind, elements%kolik

      call fem_element_center(nodes, elements, el, xc, yc, ok_center)

      if (.not. ok_center) then
        el2ncgrid(el) = -1_ikind
        cycle
      end if

      el2ncgrid(el) = find_nc_element_for_point(xc, yc, ncnodes, ncelements)

    end do

  end subroutine mapel


  function find_nc_element_for_point(xp, yp, ncnodes, ncelements) result(ncel_found)
    real(kind=rkind), intent(in) :: xp, yp
    type(node), intent(in) :: ncnodes
    type(element), intent(in) :: ncelements

    integer(kind=ikind) :: ncel_found

    integer(kind=ikind) :: ncel
    integer(kind=ikind) :: i
    integer(kind=ikind) :: inode

    real(kind=rkind), dimension(4) :: xq
    real(kind=rkind), dimension(4) :: yq

    real(kind=rkind) :: xmin, xmax
    real(kind=rkind) :: ymin, ymax

    logical :: inside
    logical :: valid_nc_element

    ncel_found = -1_ikind

    do ncel = 1_ikind, ncelements%kolik

      valid_nc_element = .true.

      do i = 1_ikind, 4_ikind

        inode = ncelements%data(ncel, i)

        if (inode <= 0_ikind .or. inode > ncnodes%kolik) then
          valid_nc_element = .false.
          exit
        end if

        xq(i) = ncnodes%data(inode, 1)
        yq(i) = ncnodes%data(inode, 2)

      end do

      if (.not. valid_nc_element) then
        cycle
      end if

      xmin = minval(xq)
      xmax = maxval(xq)
      ymin = minval(yq)
      ymax = maxval(yq)

      if (xp < xmin .or. xp > xmax) then
        cycle
      end if

      if (yp < ymin .or. yp > ymax) then
        cycle
      end if

      call point_in_quad(xp, yp, xq, yq, inside)

      if (inside) then
        ncel_found = ncel
        return
      end if

    end do

  end function find_nc_element_for_point


  subroutine fem_element_center(nodes, elements, el, xc, yc, ok)
    type(node), intent(in) :: nodes
    type(element), intent(in) :: elements
    integer(kind=ikind), intent(in) :: el

    real(kind=rkind), intent(out) :: xc
    real(kind=rkind), intent(out) :: yc
    logical, intent(out) :: ok

    integer(kind=ikind) :: i
    integer(kind=ikind) :: inode
    integer(kind=ikind) :: nnodes_el

    xc = 0.0_rkind
    yc = 0.0_rkind
    ok = .false.

    if (el <= 0_ikind .or. el > elements%kolik) then
      return
    end if

    nnodes_el = 0_ikind

    do i = 1_ikind, size(elements%data, 2, kind=ikind)

      inode = elements%data(el, i)

      if (inode <= 0_ikind) then
        cycle
      end if

      if (inode > nodes%kolik) then
        cycle
      end if

      nnodes_el = nnodes_el + 1_ikind

      xc = xc + nodes%data(inode, 1)
      yc = yc + nodes%data(inode, 2)

    end do

    if (nnodes_el <= 0_ikind) then
      return
    end if

    xc = xc / real(nnodes_el, kind=rkind)
    yc = yc / real(nnodes_el, kind=rkind)

    ok = .true.

  end subroutine fem_element_center


end module ncmap
