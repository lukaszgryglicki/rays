#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#define REAL long double
int __times;
time_t __t_start, __t_cur;

void catch_signal(int signo)
{
 time(&__t_cur);
 __t_cur -= __t_start;
 printf("Executed %d times in %d seconds, %f/s\n", 
	 __times, __t_cur, (float)__times/(float)__t_cur);
 exit(0);
}

void init()
{
 signal(SIGINT, catch_signal);
}

void test_engine()
{
 time(&__t_start);
 __times = 0;
 while (1)
   {
    __times ++;
   }
}

int main()
{
 init();
 test_engine();
 return 0;
}

