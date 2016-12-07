/*****************************************************************************/
/*random.c
 author : Hema Chandrasekaran*/
/*****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "random.h"

#define PI      3.141592654

/* seed # 1*/


#define SEED1		29
#define SEED2		3191
#define SEED3		1987

/* seed # 2*/

/*
#define SEED1		163
#define SEED2		4057
#define SEED3		2179


/* seed # 3*/

/*
#define SEED1		11
#define SEED2		4799
#define SEED3		2621

/* seed # 4*/

/*
#define SEED1		463
#define SEED2		17351
#define SEED3		20269


/* seed # 5*/
/*
#define SEED1		2141
#define SEED2		24001
#define SEED3		6257

/* seed # 6*/

/*
#define SEED1		8297
#define SEED2		19
#define SEED3		1153


/* seed # 7*/

/*
#define SEED1		27241
#define SEED2		26119
#define SEED3		3


/* seed # 8*/
/*
#define SEED1		15461
#define SEED2		31
#define SEED3		79

/* seed # 9*/
/*
#define SEED1		13757
#define SEED2		10639
#define SEED3		16433

/* seed # 10*/
/*
#define SEED1		18401
#define SEED2		5903
#define SEED3		1627
*/



struct RandomNum
{
	int ix, iy, iz;
};

/*static struct RandomNum  random, *randptr = &random;*/
static struct RandomNum  random1 = { SEED1,SEED2,SEED3},*randptr = &random1;

float rand1(void)
{
	int ixx, iyy, izz, itemp;
	float temp, uniform;

	ixx = randptr->ix/177;
	randptr->ix  = 171*(randptr->ix % 177) - 2*ixx;

	if ( randptr->ix < 0)
	{
		randptr->ix +=  30269;
	}
	iyy = randptr->iy/176;
	randptr->iy  = 172*(randptr->iy % 176) - 35*iyy;

	if ( randptr->iy < 0)
	{
		randptr->iy += 30307;
	}
	izz = randptr->iz/178;
	randptr->iz  = 170*(randptr->iz % 178) - 63*izz;

	if ( randptr->iz < 0)
	{
		randptr->iz += 30323;
	}

	temp = (float)((float)randptr->ix/30269.0 + (float)randptr->iy/30307.0 +
				(float)randptr->iz/30323.0);
	itemp = (int)temp;

	uniform = temp - itemp;
	return(uniform);
}
/******************************************************************************/
	/* function slete() generates a gaussian random process by invoking
		rand1() and manipulating the sample returned by it. */
/******************************************************************************/
float slete(float std, float mean)
{
	float rand1(void);
	float gaussian, temp1, temp2;

	temp1 =  (float)cos(2*PI*rand1());
	temp2 =  (float)sqrt( -2.0 * log(rand1()));
	gaussian = mean +  std * temp1 * temp2;

	return(gaussian);
}
/***************************************************************************/


