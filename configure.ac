#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.71])
AC_INIT
AC_CONFIG_SRCDIR([lfe])
AC_CONFIG_SRCDIR([src/lfe.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_AUX_DIR([build/autoconf])
AC_CONFIG_FILES([src/Makevars])
AC_CONFIG_MACRO_DIRS([tools])
dnl AC_LANG(C)
AC_ARG_ENABLE([threads],
        [AS_HELP_STRING([--disable-threads],[do not use threads])],
        [do_threads=$enableval],[do_threads=yes])

AC_ARG_ENABLE([huge],
	[AS_HELP_STRING([--enable-huge],[use 64 bit integers])],
	[AC_DEFINE([HUGE_INT], [1], [64bit])],
	[])

dnl local debug 
AC_ARG_ENABLE([pedantry],,
	[DFLAGS='-Wall -Wextra -Wpedantic'],
	[])

AC_SUBST(DFLAGS)
# pick up R's compilers, we need to use that one
CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
LDFLAGS=`${R_HOME}/bin/R CMD config LDFLAGS`

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

if test "$do_threads" = "yes"; then
  AX_PTHREAD([], [
  AC_MSG_WARN([No threads on this architecture])
  AC_DEFINE([NOTHREADS], [1], [no threads])
  do_threads="no"])
else
  AC_DEFINE([NOTHREADS], [1], [no threads])
fi
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)
# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

if test "$do_threads" = "yes"; then
# bugger fails on macosx
AC_MSG_CHECKING(for timed semaphore wait)
# we may use these when checking for semaphores
LDFLAGS="$LDFLAGS $PTHREAD_LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
AC_RUN_IFELSE([AC_LANG_SOURCE([[#include <semaphore.h>
 #include <stdlib.h>
 #include <stdio.h>
 #include <errno.h>
 #include <time.h>
  int main(int argc, char *argv[]) {
    sem_t foo;
    int stat;
    struct timespec tmo = {time(NULL),1000000};
    if(sem_init(&foo,0,1) != 0) {
      perror("sem_init");
      exit(1);
    }
    stat = sem_timedwait(&foo,&tmo);
    if(stat != 0 && stat != ETIMEDOUT) {
      perror("sem_timedwait");
      exit(1);
    }
    exit(0);
  }
]])],[AC_MSG_RESULT(yes)
AC_DEFINE(HAVE_SEM, 1, timed sem wait)],[AC_MSG_RESULT(no)
],[])

# check for capability of setting thread name
AC_SEARCH_LIBS([pthread_setname_np], [], AC_DEFINE([HAVE_THREADNAME],1, [chg thread name]), [], [$PTHREAD_CFLAGS])
AC_SUBST(HAVE_THREADNAME)
fi

AC_OUTPUT
