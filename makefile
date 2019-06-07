CC = gcc
DEPS = neuro.h statebuf.h
OBJ = statebuf.o sim.o
TARGET = sim

#%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $<

#statebuf.o: statebuf.h statebuf.c
#	$(CC) -c statebuf.h statebuf.c

#$(TARGET): $(OBJ) neuro.h
#	$(CC) -o $@ $^ -lm

sim: sim.c neuro.h neuro.c statebuf.c statebuf.h interpolation.c interpolation.h
	$(CC) -o $@ $^ -lm

clean: 
	rm -f $(OBJ) $(TARGET) *~
