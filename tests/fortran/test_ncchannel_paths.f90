program test_ncchannel_paths
  use typy
  use ncglobvars
  use init_netcdf, only: readchannel
  implicit none

  integer(kind=ikind) :: i, node_start, edge_index
  integer :: status
  character(len=32) :: argument

  call get_command_argument(1, argument, status=status)
  if (status /= 0) error stop "missing channel-count argument"
  read(argument, *) channel_count

  call readchannel()

  if (channel_nd%kolik /= 3_ikind + 2_ikind*(channel_count-1_ikind)) &
    error stop "wrong combined channel node count"
  if (channel_el%kolik /= 2_ikind + channel_count-1_ikind) &
    error stop "wrong combined channel segment count"

  call check_edge(1_ikind, 1_ikind, 2_ikind)
  call check_edge(2_ikind, 2_ikind, 3_ikind)

  do i = 2_ikind, channel_count
    node_start = 4_ikind + 2_ikind*(i-2_ikind)
    edge_index = i + 1_ikind
    call check_edge(edge_index, node_start, node_start+1_ikind)
    if (abs(channel_nd%data(node_start,1) - 100.0_rkind*i) > 1.0e-12_rkind) &
      error stop "channel file order or coordinates are wrong"
  end do

  print *, "multiple channel path checks: OK"

contains

  subroutine check_edge(index, first_node, second_node)
    integer(kind=ikind), intent(in) :: index, first_node, second_node

    if (channel_el%data(index,1) /= first_node .or. &
        channel_el%data(index,2) /= second_node) then
      error stop "wrong channel connectivity or artificial inter-path edge"
    end if
  end subroutine check_edge

end program test_ncchannel_paths
