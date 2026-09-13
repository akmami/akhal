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

PREFIX      ?= /usr/local
BINDIR      := $(PREFIX)/bin
LIBDIR      := $(PREFIX)/lib
INCLUDEDIR  := $(PREFIX)/include
DATADIR     := $(PREFIX)/share
BASHCOMPDIR := $(DATADIR)/bash-completion/completions
ZSHCOMPDIR  := $(DATADIR)/zsh/site-functions

INSTALL         := install
INSTALL_PROGRAM := $(INSTALL) -m 755
INSTALL_DATA    := $(INSTALL) -m 644

# library
LIB_SRC := $(wildcard src/lib/*.c)
LIB_OBJ := $(LIB_SRC:.c=.o)

# CLI (all commands + dispatcher)
CLI_SRC := src/cli/main.c \
		   src/cli/cmd_stats.c \
		   src/cli/cmd_extract.c \
		   src/cli/cmd_compare.c \
		   src/cli/cmd_parse.c \
		   src/cli/cmd_compact.c \
		   src/cli/cmd_sort.c \
		   src/cli/cmd_rank.c \
		   src/cli/cmd_vg2gfa.c \
		   src/cli/cmd_gfa2rgfa.c \
		   src/cli/cmd_gfa2dot.c \
		   src/cli/cmd_gaf2sam.c \
		   src/cli/cmd_sampoke.c \
		   src/cli/cmd_annotate.c \
		   src/cli/cmd_annotget.c
CLI_OBJ := $(CLI_SRC:.c=.o)

ALL_OBJ := $(CLI_OBJ)

.PHONY: all clean install uninstall
all: $(TARGET)

$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^
	@mkdir -p lib && mv $(LIB) lib/

$(TARGET): $(ALL_OBJ) $(LIB)
	$(CC) $(CFLAGS) -o $@ $(ALL_OBJ) lib/$(LIB) $(LDLIBS)
	@rm -f $(LIB_OBJ) $(CLI_OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

install: all
	$(INSTALL) -d $(DESTDIR)$(BINDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(INCLUDEDIR)/akhal $(DESTDIR)$(BASHCOMPDIR) $(DESTDIR)$(ZSHCOMPDIR)
	$(INSTALL_PROGRAM) $(TARGET) $(DESTDIR)$(BINDIR)/
	$(INSTALL_DATA) completions/akhal.bash $(DESTDIR)$(BASHCOMPDIR)/$(TARGET)
	$(INSTALL_DATA) completions/_akhal $(DESTDIR)$(ZSHCOMPDIR)/_$(TARGET)
	@[ "$(DESTDIR)$(LIBDIR)" -ef lib ] || $(INSTALL_DATA) lib/$(LIB) $(DESTDIR)$(LIBDIR)/
	@[ "$(DESTDIR)$(INCLUDEDIR)/akhal" -ef include/akhal ] || $(INSTALL_DATA) include/akhal/*.h $(DESTDIR)$(INCLUDEDIR)/akhal/
	@# zsh scans this directory, so completion is in place - but compinit caches
	@# what it found, and oh-my-zsh reuses that cache for a day
	@if [ -z "$(DESTDIR)" ] && command -v zsh >/dev/null 2>&1 && \
	   zsh -c 'print -l $$fpath' 2>/dev/null | grep -qx "$(ZSHCOMPDIR)"; then \
		echo "if a new tab does not complete: rm -f ~/.zcompdump* && exec zsh"; \
	fi

uninstall:
	rm -f  $(DESTDIR)$(BINDIR)/$(TARGET)
	rm -f  $(DESTDIR)$(BASHCOMPDIR)/$(TARGET)
	rm -f  $(DESTDIR)$(ZSHCOMPDIR)/_$(TARGET)
	[ "$(DESTDIR)$(LIBDIR)" -ef lib ] || rm -f $(DESTDIR)$(LIBDIR)/$(LIB)
	[ "$(DESTDIR)$(INCLUDEDIR)/akhal" -ef include/akhal ] || rm -rf $(DESTDIR)$(INCLUDEDIR)/akhal

clean:
	rm -f $(TARGET) lib/$(LIB) $(LIB_OBJ) $(CLI_OBJ)
