all: Alpha

CC=g++
DEBUG=
FLAGS=-Wall
LIB=#-L/usr/
INC=-I/home/winnen/Dropbox/BrianChen/code -I/home/winnen/Dropbox/BrianChen/code/determinant -I/home/winnen/Dropbox/BrianChen/surfaceExtractorSource/06.29.2012/source
OBJ=vertex.o Geom.o main.o mathfunc.o
#triIntersect.o mathlib.o funclib.o


Alpha: $(OBJ)
	$(CC) $(DEBUG) $(FLAGS) -o Alpha $(OBJ)

main.o: main.cpp Geom.h vertex.h
	$(CC) -c $(FLAGS) $(INC) $(LIB) main.cpp

Geom.o: Geom.h Geom.cpp vertex.h# Atom.h triIntersect.h
	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) Geom.cpp

#Atom.o: Atom.h Atom.cpp vertex.h
#	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) Atom.cpp

vertex.o: vertex.h vertex.cpp
	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) vertex.cpp

#triIntersect.o: triIntersect.h triIntersect.cpp
#	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) triIntersect.cpp

#mathlib.o: mathlib.h mathlib.cpp
#	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) mathlib.cpp

mathfunc.o: mathfunc.h mathfunc.cpp
	$(CC) -c $(DEBUG) $(FLAGS) $(INC) $(LIB) mathfunc.cpp

clean:
	rm -rf *.o Alpha ~* *~

#valgrind is on corona --leakcheck=yes