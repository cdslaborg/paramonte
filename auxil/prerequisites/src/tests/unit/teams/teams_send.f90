program test_teams1
  use iso_fortran_env, only : team_type
  implicit none

  type(team_type) :: team, initial
  integer :: x[*], me

  me = this_image()
  x = me
  
  form team (mod(me,2)+1,team)

  change team (team)

  if(me == 1) x[2, team = initial] = me

  end team

  if(me == 2 .and. x == 1) print *,"Test passed."
  
end program
