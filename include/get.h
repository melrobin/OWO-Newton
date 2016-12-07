/* get.h - prototypes for get.c functions and common macros */

char *get_string(char *);
int get_int(char *,int,int);
double get_float(char *,double,double);

/* MIN, MAX, ROUND macros */

#define MAX(a,b)    (((a) > (b)) ? (a) : (b))
#define MIN(a,b)    (((a) < (b)) ? (a) : (b))
#define ROUND(a)  (((a) < 0) ? (int)((a)-0.5) : (int)((a)+0.5))
