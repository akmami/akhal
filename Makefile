# akhal - assembly graph analysis toolkit
#
# Two build products, mirroring the htslib + samtools split:
#   libakhal.a   the reusable library  (src/lib/*.c, public headers in include/)
#   akhal        the command-line tool (src/cli/*.c, links against libakhal.a)
#

CC      := gcc
CFLAGS  := -O3 -Wall -Wextra -Wpedantic -Iinclude -Isrc/cli
LDLIBS  := -lm -lz
AR      := ar

TARGET  := akhal
LIB     := libakhal.a

# library
LIB_SRC := $(wildcard src/lib/*.c)
LIB_OBJ := $(LIB_SRC:.c=.o)

# CLI (all commands + dispatcher)
CLI_SRC := src/cli/main.c \
		   src/cli/cmd_stats.c \
		   src/cli/cmd_extract.c \
		   src/cli/cmd_parse.c \
		   src/cli/cmd_sort.c \
		   src/cli/cmd_rank.c \
		   src/cli/cmd_vg2gfa.c \
		   src/cli/cmd_gaf2sam.c \
		   src/cli/cmd_sampoke.c \
		   src/cli/cmd_annotate.c \
		   src/cli/cmd_annotget.c
CLI_OBJ := $(CLI_SRC:.c=.o)

ALL_OBJ := $(CLI_OBJ)

.PHONY: all clean
all: $(TARGET)

$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^
	@mkdir -p lib && mv $(LIB) lib/

$(TARGET): $(ALL_OBJ) $(LIB)
	$(CC) $(CFLAGS) -o $@ $(ALL_OBJ) lib/$(LIB) $(LDLIBS)
	@rm -f $(LIB_OBJ) $(CLI_OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) lib/$(LIB) $(LIB_OBJ) $(CLI_OBJ)
