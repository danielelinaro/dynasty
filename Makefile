CC=gcc
CFLAGS=-I/Users/daniele/local/include -I/usr/include -Wall -O3
LDFLAGS=-L/Users/daniele/local/lib
LIBS=-lsundials_cvode -lsundials_nvecserial -lm
OBJS=dynasty.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dynasty : $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)
	
clean :
	rm -f dynasty *.o
