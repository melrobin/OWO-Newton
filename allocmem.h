/* allocmem.h - prototypes for allocmem.c functions  */

extern FILE *OpenFile(char *, char *);
extern double  *FarAllocateMemory(int );
extern double  *AllocateMemory(int );
extern double  **FarAllocateMatrixMemory(int, int);
extern double  **FarAllocateDMatrixMemory(int, int);

