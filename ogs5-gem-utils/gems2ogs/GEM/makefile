src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
CXX = g++
LDFLAGS =
CPPFLAGS = 
CXXFLAGS = -c -g -O3 -DIPMGEMPLUGIN -D__unix
LIB =libgem.a
LINK = ar

$(LIB): $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^
	$(LINK) -rvs $(LIB) $(obj)

clean:
	rm -f $(obj) $(LIB)


