CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o position.o pdf.o metropolis.o statistic.o
	$(CC) random.o main.o position.o pdf.o metropolis.o statistic.o -o main.exe
main.o : main.cpp random.h measure.h position.h pdf.h metropolis.h statistic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
position.o : position.cpp position.h
	$(CC) -c position.cpp -o position.o $(CFLAGS)
pdf.o : pdf.cpp pdf.h position.h
	$(CC) -c pdf.cpp -o pdf.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h position.h pdf.h measure.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h random.h measure.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)

clean:
	rm *.o main.exe
