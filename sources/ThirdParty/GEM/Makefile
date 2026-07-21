CXX = g++
CXXFLAGS =   -c -g -O3 -DIPMGEMPLUGIN -DSPDLOG_FMT_EXTERNAL=ON -DNO_THERMOFUN -D__unix 
AR = ar
ARFLAGS=rcs
TARGET = libgem.a

# Sucht automatisch alle .cpp Dateien im aktuellen Verzeichnis
SRCS = $(wildcard *.cpp)

# Erzeugt aus den .cpp Dateinamen die entsprechenden .o (Objekt) Dateinamen
OBJS = $(SRCS:.cpp=.o)

# Standard-Ziel (wird ausgeführt, wenn du nur 'make' tippst)
all: $(TARGET)

# Regel zum Erstellen der statischen Bibliothek aus den .o Dateien
$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)

# Regel zum Kompilieren der .cpp Dateien in .o Dateien
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Bereinigt das Verzeichnis von Kompilaten
clean:
	rm -f $(OBJS) $(TARGET)

# Verhindert Konflikte, falls Dateien existieren, die wie die Targets heissen
.PHONY: all clean
