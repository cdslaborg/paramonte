program teams_subset
  !! author: Nathan Weeks
  !!
  !! Regression test for prior potential deadlock in change team/end team
  use iso_fortran_env, only : team_type
  use oc_assertions_interface, only : assert
  implicit none

  type(team_type) team

  associate( initial_team_image => this_image())
  associate( myteam => merge(1,2,any(initial_team_image==[1,num_images()])) )

    form team (myteam, team)

    if (myteam == 1) then
      block
        integer max_image(2), min_image(2)
        change team(team)
          max_image = [initial_team_image, this_image()]
          call co_max(max_image)
          min_image = [initial_team_image, this_image()]
          call co_min(min_image)
        end team
        call assert( all(min_image==1), "minimum image number in team 1 is 1" )
        call assert( all(max_image==[num_images(), 2]), "maximum image numbers on team 1 is num_images() or 2")
      end block
    end if

    sync all

    if (initial_team_image == 1) print *,"Test passed."

  end associate; end associate

end program
