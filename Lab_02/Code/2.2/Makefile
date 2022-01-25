CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o randomwalk.o statistic.o
	$(CC) random.o main.o randomwalk.o statistic.o -o main.exe
main.o : main.cpp random.h measure.h randomwalk.h statistic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
randomwalk.o: randomwalk.cpp randomwalk.h random.h
	$(CC) -c randomwalk.cpp -o randomwalk.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h random.h measure.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)

clean:
	rm *.o main.exe
