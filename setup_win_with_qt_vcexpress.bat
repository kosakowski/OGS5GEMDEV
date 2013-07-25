::: This calls setup_win after aliasing devenv.exe to vcexpress.exe :::
set USE_VC_EXPRESS=TRUE
set CMAKE_GENERATOR=-G "Visual Studio 9 2008"
call setup_win_with_qt.bat