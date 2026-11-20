CC     ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L

all: libregion_timer.a

libregion_timer.a: region_timer.o
	ar rcs $@ $^

test: test_region_timer
	./test_region_timer

test_region_timer: test_region_timer.c libregion_timer.a
	$(CC) $(CFLAGS) $< -L. -lregion_timer -lm -o $@

clean:
	rm -f *.o *.a test_region_timer
