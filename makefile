CC = gcc
CFLAGS = -Wall -g
# Libraries: math.h
LDFLAGS = -lm

# Directories
OBJDIR = obj
SRCDIR = src
RESDIR = results
PLOTDIR = plots

# Files
TARGET = sim
SRCS := $(wildcard $(SRCDIR)/*.c)
SRCS := $(filter-out $(SRCDIR)/testing.c,$(SRCS))

# Canonical or prescient implementation?
# Use 'make c=1 <rule>' for canonical implementation.
ifdef c
	CFLAGS += -DCANONICAL
	SRCS := $(filter-out $(SRCDIR)/prescient.c,$(SRCS))
else
	SRCS := $(filter-out $(SRCDIR)/canonical.c,$(SRCS))
endif

OBJS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))
PLOTS := $(wildcard $(PLOTDIR)/*)
PLOTS := $(patsubst $(PLOTDIR)/%,%,$(PLOTS))

all: dir clean $(OBJDIR)/$(TARGET)

dir:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(RESDIR)

run:
	@$(OBJDIR)/$(TARGET)

debug:
	valgrind -v --leak-check=full --track-origins=yes --show-leak-kinds=all $(OBJDIR)/$(TARGET)

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
	@rm -f $(OBJS) $(TARGET) *~
	@echo "Cleaned up obj directory."