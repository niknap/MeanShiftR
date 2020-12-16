// Copyright (C) 2018 Dr. Nikolai Knapp, UFZ
//
// This file is part of the MeanShiftR R package.
//
// The MeanShiftR R package is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// MeanShiftR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MeanShiftR If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// Collection of all the little functions used by the main functions

// Function to check whether a point [PointX, PointY, PointZ] is within a cylider of a given radius
// and height from the center point of the top circle [TopX, TopY, TopZ]
bool InCylinder(double PointX, double PointY, double PointZ, double Radius, double Height, double CtrX, double CtrY, double CtrZ){
  if ((pow((PointX - CtrX), 2.0) + pow((PointY - CtrY), 2.0) <= pow(Radius, 2.0)) && (PointZ >= (CtrZ - (0.5*Height))) && (PointZ <= (CtrZ + (0.5*Height))) == true) {
    return true;
  }	else {
    return false;
  }
}

// Help functions for vertical filter
double VerticalDistance(double Height, double CtrZ, double PointZ){
  double BottomDistance = (double) std::abs((CtrZ-Height/4-PointZ)/(3*Height/8));
  double TopDistance = (double) std::abs((CtrZ+Height/2-PointZ)/(3*Height/8));
  double MinDistance = std::min(BottomDistance, TopDistance);
  return MinDistance;
}

//Equivalent R code
//distx <- function(h, CtrZ, PointZ){
//  bottomdist <- abs((CtrZ-h/4-PointZ)/(3*h/8))
//  topdist <- abs((CtrZ+h/2-PointZ)/(3*h/8))
//  mindist <- pmin(bottomdist, topdist)
//  return(mindist)
//}

double VerticalMask(double Height, double CtrZ, double PointZ){
  if((PointZ >= CtrZ-Height/4) && (PointZ <= CtrZ+Height/2)){
    return 1;
  }
  else
  {
    return 0;
  }
}

//Equivalent R code
//maskx <- function(h, CtrZ, PointZ){
//  maskvec <- ifelse(PointZ >= CtrZ-h/4 & PointZ <= CtrZ+h/2, 1, 0)
//  return(maskvec)
//}

// Epanechnikov function for vertical filter
double EpanechnikovFunction(double Height, double CtrZ, double PointZ){
  double Result = VerticalMask(Height, CtrZ, PointZ)*(1-(pow((1-VerticalDistance(Height, CtrZ, PointZ)), 2.0)));
  return Result;
}

//Equivalent R code
//Epanechnikov <- function(h, CtrZ, PointZ){
//  output <- maskx(h, CtrZ, PointZ)*(1-(1-distx(h, CtrZ, PointZ))^2)
//  return(output)
//}

// Gauss function for horizontal filter
double GaussFunction(double Width, double CtrX, double CtrY, double PointX, double PointY){
  double Distance = pow((pow((PointX-CtrX), 2.0)+pow((PointY-CtrY), 2.0)), 0.5);
  double NormDistance = Distance/Width;
  double Result = std::exp(-5.0*pow(NormDistance, 2.0));
  return Result;
}

//Equivalent R code
//gauss <- function(w, CtrX, CtrY, PointX, PointY){
//  distance <- ((PointX-CtrX)^2+(PointY-CtrY)^2)^0.5
//  norm.distance <- distance/w
//  output <- exp(-5*norm.distance^2)
//  return(output)
//}
