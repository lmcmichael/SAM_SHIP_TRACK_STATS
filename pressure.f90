subroutine pressure

! call a pressure solver

use grid
implicit none

call t_startf ('pressure')

! for big runs only (nx >> nzm)
if(RUN3D) then
  call pressure_big
else
  call pressure_orig
end if

call t_stopf ('pressure')

end
