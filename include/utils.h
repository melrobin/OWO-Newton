#include <stdio.h>
#define	 MAXFEATURES		100
#define  MAXHIDDENLAYERS	2
#define	 SIGDIGITS			9
#define	 FILENAMELEN		50
#define	 TRUE				1
#define	 FALSE				0

typedef struct  UserEntry
{
	int 		Features;
	int 		DesiredOutputs;
	int 		MaxPatterns;
	char		SaveOutput;
	FILE		*fpwt, *fpfeat, *fplog, *fplist, *fpout;
	char		*WtFile, *ListFile, *FeatFile, *LogFile, *OutFile;
} USERINPUT;

struct TopologyInfo
{
	FILE		*fpTopology;
	char		*TopologyFile;
	int			oldTopology;
	int			HiddenLayers;
	int			HiddenLayerUnits[MAXHIDDENLAYERS];
	int 		LayerUnits[MAXHIDDENLAYERS+2];
};

struct MiscArrays
{
	double 		*Inputs, *Outputs;
	double  	****Wts;
	double  	**Thresholds;
	double  	**Nets;
	double  	**LayerOutputs;
	double		MSError[MAXFEATURES];
	double		TotalError;
};

struct Origin
{
	struct UserEntry 	*userptr;
	struct TopologyInfo  *tptr;
	struct MiscArrays    *arrayptr;
};

typedef struct _INPUTINFO
	{
	  int num;
	  double *mean;
	  double *std;
        }
        XINFO;

