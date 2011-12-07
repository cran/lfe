#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.59])
AC_INIT(lfe, {version} , frisch@gaure.no)
AC_CONFIG_SRCDIR([src/lfe.c])

AC_ARG_ENABLE([threads],
        [AS_HELP_STRING([--disable-threads],[do not use threads])],
        [do_threads=$enableval],[do_threads=yes])

# pick up R's compilers, we need to use that one
CC=`${R_HOME}/bin/R CMD config CC`

# Checks for programs.
AC_PROG_CC
AC_PROG_AWK

# Checks for libraries.

if test "$do_threads" = "yes"; then
  ACX_PTHREAD([],[PTHREAD_CFLAGS=-DNOTHREADS])
else
  PTHREAD_CFLAGS=-DNOTHREADS
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# bugger fails on macosx
AC_MSG_CHECKING(for timed semaphore wait)
# we may use these when checking for semaphores
CFLAGS=$PTHREAD_CFLAGS
LDFLAGS=$PTHREAD_LIBS

AC_TRY_RUN(dnl
[#include <semaphore.h>
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
],
AC_MSG_RESULT(yes)
SEM_FLAGS=-DHAVE_SEM,
AC_MSG_RESULT(no)
SEM_FLAGS=
)

AC_SUBST(SEM_FLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT


