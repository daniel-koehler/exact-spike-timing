CC = gcc
CFLAGS = -Wall
# Libraries: math.h
LDFLAGS = -lm

OBJDIR = obj
SRCDIR = src

TARGET = sim
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))

all: dir $(OBJDIR)/$(TARGET)
	
dir:
	mkdir -p $(OBJDIR)

$(OBJDIR)/$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $< $(LDFLAGS)

.PHONY: clean
clean: 
	rm -f $(OBJS) $(TARGET) *~
	@echo "Cleaned up build directory."

