::: This calls setup_win after aliasing devenv.exe to vcexpress.exe :::
doskey devenv=vcexpress $1 $2 $3
set CMAKE_GENERATOR=-G "Visual Studio 9 2008"
call setup_win_cl_only.bat