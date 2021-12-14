! SPDX-License-Identifier: MIT
submodule (git_version) git_version_impl
  implicit none
#ifndef GIT_SHA1
#define GIT_SHA1 "unknown"
#endif

#ifndef GIT_STATE
#define GIT_STATE "unknown"
#endif

#ifndef GIT_VERSION
#define GIT_VERSION "unknown"
#endif
contains
  module procedure get_git_version
    integer, parameter :: max_length = 40
    integer :: length

    length = min(max_length, len(GIT_VERSION))
    allocate(character(length)::get_git_version)
    get_git_version = GIT_VERSION(1:length)
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
    get_git_hash = GIT_SHA1(1:length)
  end procedure get_git_hash

  module procedure get_git_state
    if (GIT_STATE == "clean") then
      get_git_state = ""
    else
      get_git_state = "-dirty"
    endif
  end procedure get_git_state
end submodule git_version_impl
