#include "region_timer.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>

static void sleep_ms(int ms)
{
   struct timespec ts = {0, ms * 1000000L};
   nanosleep(&ts, NULL);
}

/* Only the outermost start/end pair of a recursive region is timed */
static void recurse(int n)
{
   region_start("recursive");
   if (n > 0) recurse(n - 1);
   else sleep_ms(20);
   region_end("recursive");
}

int main(void)
{
   int i;
   region_start("main");
   for (i = 0; i < 5; i++) {
      region_start("fields");
      sleep_ms(20);
      region_end("fields");
   }
   region_start("gke");
   sleep_ms(50);
   region_end("gke");
   recurse(3);
   region_end("main");

   assert(region_total_time("fields") >= 0.1 && region_total_time("fields") < 0.2);
   assert(region_total_time("gke") >= 0.05 && region_total_time("gke") < 0.1);
   assert(region_total_time("recursive") >= 0.02 && region_total_time("recursive") < 0.04);
   assert(region_total_time("main") >= region_total_time("fields") + region_total_time("gke") + region_total_time("recursive"));
   assert(region_total_time("never_entered") == 0.0);

   region_end("never_started");        /* prints a warning, is otherwise ignored */
   puts("OK");
   return 0;
}
