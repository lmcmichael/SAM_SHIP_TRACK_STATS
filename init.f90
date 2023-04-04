! Initialize some arrays from vars module:

subroutine init()

  use vars
  use tracers
  use rad, only: qrad

  implicit none

  u = 0.
  v = 0.
  w = 0.
  t = 0.
  qrad = 0.

  zb_old = 0.
  qt_old = 0.
  tl_old = 0.
  cl_old = 0.
  cb_old = 0.
  lwp_old = 0.
  tlcl_old = 0.
  tlsc_old = 0.
  qtcl_old = 0.
  qtsc_old = 0.

  fzero = 0.
 
  ttend = 0.
  qtend = 0.
  wsub = 0.
  unudge = 0.
  vnudge = 0.
  tnudge = 0.
  qnudge = 0.
  qlsvadv = 0.
  tlsvadv = 0.
  ulsvadv = 0.
  vlsvadv = 0.
  qstor = 0.
  tstor = 0.
  ustor = 0.
  vstor = 0.

  radlwup = 0.
  radlwdn = 0.
  radswup = 0.
  radswdn = 0.
  radqrlw = 0.
  radqrsw = 0.

  tlat = 0. 
  tlatqi = 0.
  tadv = 0.
  tdiff = 0.
  qifall = 0.
  qpfall = 0.

  trwle = 0.
  trwsb = 0.
  tradv = 0.
  trdiff = 0.
  trphys = 0.

  gamt0 = 0.
  gamq0 = 0.

end subroutine init



