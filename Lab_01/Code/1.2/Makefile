CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o statistic.o
	$(CC) random.o main.o statistic.o -o main.exe
main.o : main.cpp random.h statistic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
statistic.o : statistic.cpp statistic.h
	$(CC) -c statistic.cpp -o statistic.o $(CFLAGS)

clean:
		rm *.o main.exe
