CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lGL -lGLU -lglut -lpng
endif
	
RM = /bin/rm -f 
all: main 
main: as1.o 
	$(CC) $(CFLAGS) -o as1 as1.o $(LDFLAGS) 
as1.o: as1.cpp
	$(CC) $(CFLAGS) -c as1.cpp -o as1.o
clean: 
	$(RM) *.o as1
 


