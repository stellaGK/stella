/*
 * region_timer: named wall-clock timers.
 *
 *   region_start("fields");
 *   ...
 *   region_end("fields");
 *   t = region_total_time("fields");   // accumulated seconds
 *
 * Semantics
 *  - Regions are created on first region_start. Names are case-sensitive.
 *  - Regions may be nested/overlap freely (each has its own clock).
 *  - Recursive re-entry of the SAME region is supported: only the outermost
 *    start/end pair is timed.
 *  - region_end without a matching start prints a warning to stderr and is ignored.
 *  - NOT thread-safe: call from one thread (e.g. outside OpenMP parallel
 *    regions). Every MPI rank keeps its own timers.
 *  - Clock: CLOCK_MONOTONIC.
 */
#ifndef REGION_TIMER_H
#define REGION_TIMER_H

#ifdef __cplusplus
extern "C" {
#endif

void   region_start(const char *name);
void   region_end(const char *name);

/* Accumulated time in seconds; 0.0 for a region that was never entered. */
double region_total_time(const char *name);

#ifdef __cplusplus
}
#endif
#endif
