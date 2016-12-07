/*
 * * This program calculates the time required to
 * * execute the program specified as its first argument.
 * * The time is printed in seconds, on standard out.
 * */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

#define BILLION  1000000000L;
double tic()
{
    struct timespec start;
    double accum;

    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) {
      perror( "clock gettime" );
      return EXIT_FAILURE;
    }

    accum =   start.tv_sec + (double) start.tv_nsec/(double)BILLION;
    return accum;
}
double toc()
{
   struct timespec stop;
   double accum;

   if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) 
    {
      perror( "clock gettime" );
      return EXIT_FAILURE;
    }
    accum =   stop.tv_sec + (double) stop.tv_nsec/(double)BILLION;
    return accum;

}
int main( int argc, char** argv )
  {
    struct timespec start, stop;
    double accum;
    int i;

    if( clock_gettime( CLOCK_REALTIME, &start) == -1 ) {
      perror( "clock gettime" );
      return EXIT_FAILURE;
    }

    for (i=0;i<1000000;i++);
    if( clock_gettime( CLOCK_REALTIME, &stop) == -1 ) {
      perror( "clock gettime" );
      return EXIT_FAILURE;
    }

    accum = ( stop.tv_sec - start.tv_sec )
            + (double)( stop.tv_nsec - start.tv_nsec )
              / (double)BILLION;
    printf( "%lf\n", accum );
    return EXIT_SUCCESS;
  }
