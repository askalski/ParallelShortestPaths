#ifndef __tests_h__
#define __tests_h__

void noop();

#ifdef MYDEBUG
#define deb(...) fprintf (stderr, __VA_ARGS__)
#else
#define deb(...) noop()
#endif

#endif
