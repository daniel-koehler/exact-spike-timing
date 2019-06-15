CC = gcc
CFLAGS = -Wall -g
# Libraries: math.h
LDFLAGS = -lm

OBJDIR = obj
SRCDIR = src
RESDIR = results
PLOTDIR = plots

TARGET = sim
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))
PLOTS_ = $(wildcard $(PLOTDIR)/*)
PLOTS = $(patsubst $(PLOTDIR)/%,%,$(PLOTS_))

all: dir $(OBJDIR)/$(TARGET)
	
test:
	@echo $(PLOTS_)
dir:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(RESDIR)

run:
	$(OBJDIR)/$(TARGET)

debug:
	valgrind -v --leak-check=full $(OBJDIR)/$(TARGET)

$(OBJDIR)/$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $< $(LDFLAGS)

# make graphics using gnuplot
.PHONY: plots
plots:
	@echo "Available plots:\n"$(PLOTS)

$(PLOTS):
	gnuplot $(PLOTDIR)/$@

.PHONY: clean
clean: 
	rm -f $(OBJS) $(TARGET) *~
	@echo "Cleaned up obj directory."