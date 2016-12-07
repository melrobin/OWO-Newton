#include <valarray>
struct MiscArrays
{
double 	*ColumnWts,   *p, *g,  *ColumnWts2;
double 	*ColumnCrossCor;
double 	*Inputs, *Outputs, *HidInputs;
       std::valarray<double> OutputMS;
double	****Wts;
double	**Thresholds;
double	**Nets;
double	**Mean, **Var, *InputMean, *InputVar;
double	**LayerOutputs;
double	MinMse, ErrorCG;
int 	HwoOrOwo, HwoLayer, FirstTime;
double 	chisq, alamda;
int 	*ia, ma, ErrorIncrease, *PickHiddenUnits;
double	**covar, **alpha, *a;
};
 
extern void ConjugateGradient(int);
void ConjugateGradientLM(double**covar, int mfit, double*oneda);
