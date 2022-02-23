CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o TSP.o random.o
	$(CC) random.o TSP.o main.o -o main.exe
main.o : main.cpp TSP.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TSP.o: TSP.cpp TSP.h random.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)

clean:
	rm *.o main.exe *.dat
