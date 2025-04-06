# programs
TARGET := akhal
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

GXX := gcc
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) -o $@ $^ -lm
	rm $(OBJS)

%.o: %.c
	$(GXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
