# update GEM sources
#cp -r ../../sources/ThirdParty/GEM/*.h ./GEM/
#cp -r ../../sources/ThirdParty/GEM/*.cpp ./GEM/
#make GEM library
#cd GEM
#make
#cd ..
#compile 
g++ -g gems2ogs.cpp `wx-config --libs` `wx-config --cxxflags` -I./GEM -DIPMGEMPLUGIN  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -D__WXGTK__ -L/usr/lib64 -pthread  -L./GEM -lgem  -o gems2gsrf 


##g++ HelloWorldApp.cpp -I/usr/lib64/wx-2.8-stl/wx/include/gtk2-unicode-release-2.8 -I/usr/include/wx-2.8 -I./GEM -DIPMGEMPLUGIN -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -D__WXGTK__ -L/usr/lib64 -pthread   -L/usr/lib64/wx-2.8-stl   -lwx_gtk2u_richtext-2.8 -lwx_gtk2u_aui-2.8 -lwx_gtk2u_xrc-2.8 -lwx_gtk2u_qa-2.8 -lwx_gtk2u_html-2.8 -lwx_gtk2u_adv-2.8 -lwx_gtk2u_core-2.8 -lwx_baseu_xml-2.8 -lwx_baseu_net-2.8 -lwx_baseu-2.8 -L./GEM -lgemipm2k  -o gems2gsrf 
