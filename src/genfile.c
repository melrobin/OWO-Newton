#include <stdio.h>

int main(void)
{
   int i;
   FILE *fp;

   fp = fopen("splittest.txt","w");
   for (i=0;i<7000;i++)
      fprintf(fp,"%d %f\n",i,2.*i);
   fclose(fp);
   return(0);
}
   
