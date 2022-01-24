CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o thismeasure.o statistic.o
	$(CC) random.o main.o thismeasure.o statistic.o -o main.exe
main.o : main.cpp random.h measure.h thismeasure.h statistic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
thismeasure.o : thismeasure.cpp thismeasure.h measure.h
	$(CC) -c thismeasure.cpp -o thismeasure.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h random.h measure.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)

clean:
		rm *.o main.exe
