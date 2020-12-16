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
using namespace Rcpp;


// Declarations of all the little functions used by the main functions

bool InCylinder(double PointX, double PointY, double PointZ, double Radius, double Height, double CtrX, double CtrY, double CtrZ);
double VerticalDistance(double Height, double CtrZ, double PointZ);
double VerticalMask(double Height, double CtrZ, double PointZ);
double EpanechnikovFunction(double Height, double CtrZ, double PointZ);
double GaussFunction(double Width, double CtrX, double CtrY, double PointX, double PointY);













