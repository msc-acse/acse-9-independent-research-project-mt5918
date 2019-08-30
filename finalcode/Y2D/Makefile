all: $(SRCS)


# Source Files
SRCS = frame.c  Y.c  Ycd.c  Yfd.c  Yfrd.c\
       Yad.c Yid.c  Ymd.c  Yod.c  Yrd.c  Ysd.c  Ywd.c

# Object Files
OBJS=$(SRCS:.c=.o)

# Executable Program Name
PROG=Yf

# Compiler
CXX = gcc

# Include Directories
INCLUDES = -I. -I../include/ -I/usr/local/include

# Compiler Flags
CXXFLAGS = -O3 -ffast-math -funroll-loops -ftree-vectorize  


# Objects
%.o : %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $^

# Library Paths
LIBS = -L../ -L/usr/local/lib64 -lm\
       -L/usr/X11/lib64 -lX11 -lXext -lXt\
       -L/usr/X11R6/lib64 -lGL -lGLU


all: $(OBJS) $(PROG)

# Compilation command
$(PROG) : $(OBJS)
	$(CXX) $(INCLUDES) $(CXXFLAGS) $(OBJS) $(LIBS) -o $(PROG)

# Depend command
depend : $(SRCS)
	makedepend -- $(CXXFLAGS) -- $(INCLUDES) $(SRCS)

# Clean command
clean:
	rm -f $(PROG) $(OBJS) core
	
# DO NOT DELETE
