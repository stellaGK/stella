#include "region_timer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double now(void)
{
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC, &ts);
   return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#define NBUCKETS 256

typedef struct region {
   char          *name;
   double         total;   /* accumulated seconds */
   double         t_start; /* start stamp of the outermost open interval */
   int            depth;   /* open (nested) starts of this region */
   struct region *next;    /* next region in the same hash bucket */
} region;

static region *buckets[NBUCKETS];

static unsigned hash(const char *s)
{
   unsigned h = 5381u; /* djb2 */
   while (*s) h = h * 33u + (unsigned char)*s++;
   return h % NBUCKETS;
}

static region *find(const char *name)
{
   region *r;
   for (r = buckets[hash(name)]; r; r = r->next)
      if (strcmp(r->name, name) == 0) return r;
   return NULL;
}

static region *find_or_create(const char *name)
{
   region *r = find(name);
   unsigned h;
   if (r) return r;
   r = calloc(1, sizeof *r);
   if (r) r->name = malloc(strlen(name) + 1);
   if (!r || !r->name) { fprintf(stderr, "region_timer: out of memory\n"); abort(); }
   strcpy(r->name, name);
   h = hash(name);
   r->next = buckets[h];
   buckets[h] = r;
   return r;
}

void region_start(const char *name)
{
   region *r = find_or_create(name);
   if (r->depth++ == 0) r->t_start = now();
}

void region_end(const char *name)
{
   region *r = find(name);
   if (!r || r->depth == 0) {
      fprintf(stderr, "region_timer: region_end(\"%s\") without matching region_start\n", name);
      return;
   }
   if (--r->depth == 0) r->total += now() - r->t_start;
}

double region_total_time(const char *name)
{
   region *r = find(name);
   return r ? r->total : 0.0;
}
