VARIANT =

ifneq ($(VARIANT), )
	SUFFIX = -$(VARIANT)
else
	SUFFIX =
endif

CC ?= gcc
DEFAULTCFLAGS = -std=gnu99 -g
RELEASECFLAGS = -O3 -DNDEBUG
DEBUGCFLAGS = -O0 -fstack-protector-all -fsanitize=address -Wall
BUILDCFLAGS =  $(DEFAULTCFLAGS) $(RELEASECFLAGS) $(CFLAGS)

RINMINER_BIN = rinminer$(SUFFIX)

all : clean $(RINMINER_BIN)

$(RINMINER_BIN):
	$(CC) $(BUILDCFLAGS) -o $(RINMINER_BIN) -pthread rinminer.c

clean:
	rm -f $(RINMINER_BIN)

debug : BUILDCFLAGS =  $(DEFAULTCFLAGS) $(DEBUGCFLAGS) $(CFLAGS)
debug : all
