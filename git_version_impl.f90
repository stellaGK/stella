! SPDX-License-Identifier: MIT


submodule (git_version) git_version_impl
  implicit none
contains
  module procedure get_git_version
    integer, parameter :: max_length = 40
    integer :: length

    length = min(max_length, len("v0.5.1-370-g195f3b679-dirty"))
    allocate(character(length)::get_git_version)
    get_git_version = "v0.5.1-370-g195f3b679-dirty"(1:length)
    get_git_version = trim(get_git_version)
  end procedure get_git_version

  module procedure get_git_hash
    integer :: length

    length = 7
    if (present(length_in)) then
      if (length_in <= 40) then
        length = length_in
      end if
    end if

    allocate(character(length)::get_git_hash)
    get_git_hash = "195f3b6799303aca690d9a6fc0fd4a51749d27d6"(1:length)
  end procedure get_git_hash

  module procedure get_git_state
    if ("dirty" == "clean") then
      get_git_state = ""
    else
      get_git_state = "-dirty"
    endif
  end procedure get_git_state

  module procedure get_git_date
    get_git_date = "2024-07-26"
  end procedure get_git_date
end submodule git_version_impl
