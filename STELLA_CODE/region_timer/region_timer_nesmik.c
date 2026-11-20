/*
 * region_timer implementation backed by the nesmik profiling library
 * (https://gitlab.pm.bsc.es/beppp/nesmik). Every region_start/region_end
 * call is forwarded to nesmik_region_start/nesmik_region_stop so that
 * nesmik's own backends (Extrae, DLB-TALP, Nsys, ...) see the regions.
 *
 * region_total_time keeps working exactly as in the plain region_timer.c
 * implementation (same hash table bookkeeping) since nesmik has no public
 * API to query back accumulated time, and stella's own end-of-run timing
 * report relies on it.
 *
 * region_end("total") additionally calls nesmik_finalize(), since "total"
 * marks the end of the whole run.
 */
#include "region_timer.h"

#include <nesmik/nesmik.h>

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
static int nesmik_initialised = 0;

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
   if (!r || !r->name) { fprintf(stderr, "region_timer_nesmik: out of memory\n"); abort(); }
   strcpy(r->name, name);
   h = hash(name);
   r->next = buckets[h];
   buckets[h] = r;
   return r;
}

void region_start(const char *name)
{
   region *r;

   r = find_or_create(name);
   if (r->depth++ == 0) r->t_start = now();

   nesmik_region_start(name);
}

void region_end(const char *name)
{
   region *r = find(name);
   if (!r || r->depth == 0) {
      fprintf(stderr, "region_timer_nesmik: region_end(\"%s\") without matching region_start\n", name);
      return;
   }
   if (--r->depth == 0) r->total += now() - r->t_start;

   nesmik_region_stop(name);

   if (strcmp(name, "total") == 0) nesmik_finalize();
}

double region_total_time(const char *name)
{
   region *r = find(name);
   return r ? r->total : 0.0;
}
