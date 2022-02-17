CC = g++
CFLAGS = -Wall -O3 --std=c++11

VMC.exe : VMC.o random.o
	$(CC) random.o VMC.o -o VMC.exe
VMC.o : VMC.cpp
	$(CC) -c VMC.cpp -o VMC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o VMC.exe *.dat
