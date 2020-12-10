CC=gcc
CXX=g++
CFLAGS=-I/Users/daniele/local/include -I/usr/include -O3
LDFLAGS=-L/Users/daniele/local/lib
LIBS=-lsundials_cvode -lsundials_nvecserial -lm -lpthread
OBJS=dynasty.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dynasty : dynasty.c
	$(CC) -c -o dynasty.o $(CFLAGS) dynasty.c
	$(CC) -o dynasty $(OBJS) $(LDFLAGS) $(LIBS)
	
lib : dynasty.c
	$(CC) -g -c -o dynasty.o -fPIC -DLIB $(CFLAGS) dynasty.c
	$(CXX) -g -o dynasty.dylib -fPIC -shared $(LDFLAGS) dynasty.o $(LIBS)

clean :
	rm -f dynasty *.o *.dylib *.so
