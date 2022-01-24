CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o functions.o integrate.o statistic.o
	$(CC) random.o main.o functions.o integrate.o statistic.o -o main.exe
main.o : main.cpp random.h measure.h functions.h integrate.h statistic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)
integrate.o : integrate.cpp integrate.h
	$(CC) -c integrate.cpp -o integrate.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h random.h measure.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)

clean:
	rm *.o main.exe
