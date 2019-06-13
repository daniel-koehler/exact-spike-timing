CC = gcc
CFLAGS = -Wall
# Libraries: math.h
LDFLAGS = -lm

OBJDIR = obj
SRCDIR = src
RESPATH = results

TARGET = sim
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

all: dir $(OBJDIR)/$(TARGET)
	
dir:
	mkdir -p $(OBJDIR)
	mkdir -p $(RESPATH)

$(OBJDIR)/$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $< $(LDFLAGS)

# make graphics using gnuplot
raster: plots/raster_plot
	gnuplot $^

voltage: plots/voltage_plot
	gnuplot $^

.PHONY: clean
clean: 
	rm -f $(OBJS) $(TARGET) *~
	@echo "Cleaned up obj directory."