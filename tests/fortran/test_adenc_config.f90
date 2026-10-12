! Read the actual initialization settings without invoking DRUtES main.
program test_adenc_config
  use globals, only: drutes_config,end_time
  use ncglobvars, only: LSbank_noflow,LSconservative,LSbalance,LSsupg,LSsupg_factor,LSshock,LSshock_factor
  use ncboundary, only: read_adenc_banks
  use ncconservative, only: read_adenc_conservative
  use nchydroflow, only: read_adenc_hydro,LShydro
  use ncsupg, only: read_adenc_supg
  implicit none
  drutes_config%dimen=2
  drutes_config%it_method=0
  drutes_config%run_from_backup=.false.
  end_time=172800
  call read_adenc_supg()
  call read_adenc_banks()
  call read_adenc_conservative()
  call read_adenc_hydro()
  write(*,'(A,6(1X,I1),2(1X,ES24.16))') 'CONFIG', &
    merge(1,0,LSbank_noflow),merge(1,0,LSconservative),merge(1,0,LSbalance),merge(1,0,LShydro), &
    merge(1,0,LSsupg),merge(1,0,LSshock),LSsupg_factor,LSshock_factor
end program
