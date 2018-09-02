#include <Rcpp.h>
using namespace Rcpp;


// Declarations of all the little functions used by the main functions

bool InCylinder(double PointX, double PointY, double PointZ, double Radius, double Height, double CtrX, double CtrY, double CtrZ);
double VerticalDistance(double Height, double CtrZ, double PointZ);
double VerticalMask(double Height, double CtrZ, double PointZ);
double EpanechnikovFunction(double Height, double CtrZ, double PointZ);
double GaussFunction(double Width, double CtrX, double CtrY, double PointX, double PointY);













