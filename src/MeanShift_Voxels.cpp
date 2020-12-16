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
#include "LittleFunctionsCollection.h"
using namespace Rcpp;
using namespace std;


//' Mean shift clustering using a discrete voxel space
//'
//' @title Mean shift clustering using a discrete voxel space
//' @description
//' Adaptive mean shift clustering to delineate tree crowns from lidar point clouds. This is a version using 1-m³ voxels instead of exact point coordinates, to speed up processing.
//' @param pc Point cloud has to be in matrix format with 3-columns representing X, Y and Z and each row representing one point
//' @param H2CW_fac Factor for the ratio of height to crown width. Determines kernel diameter based on its height above ground.
//' @param H2CL_fac Factor for the ratio of height to crown length. Determines kernel height based on its height above ground.
//' @param UniformKernel Boolean to enable the application of a simple uniform kernel without distance weighting (Default False)
//' @param MaxIter Maximum number of iterations, i.e. steps that the kernel can move for each point. If centroid is not found after all iteration, the last position is assigned as centroid and the processing jumps to the next point
//' @param maxx Maximum X-coordinate
//' @param maxy Maximum Y-coordinate
//' @param maxz Maximum Z-coordinate
//' @return data.frame with X, Y and Z coordinates of each point in the point cloud and  X, Y and Z coordinates of the centroid to which the point belongs
//' @export
// [[Rcpp::export]]
DataFrame MeanShift_Voxels(NumericMatrix pc, double H2CW_fac, double H2CL_fac, bool UniformKernel=false, int MaxIter=20, int maxx=100, int maxy=100, int maxz=60){

  int nrows = pc.nrow();
  int minx = 0;
  int miny = 0;
  int minz = 0;

  NumericVector centroidx(nrows);
  NumericVector centroidy(nrows);
  NumericVector centroidz(nrows);

  std::vector<vector<vector<double> > >array3D;

  double myx;
  double myy;
  double myz;
  int vx;
  int vy;
  int vz;
  int np;
  double r;
  double d;
  double h;
  int boxminx;
  int boxmaxx;
  int boxminy;
  int boxmaxy;
  int boxminz;
  int boxmaxz;

  // Adjust the size of the array3D
  array3D.resize((int) maxx+1);
  for (int xi=0; xi<maxx+1; xi++) {
    array3D[xi].resize((int) maxy+1);
    for (int yi=0; yi<maxy+1; yi++) {
      array3D[xi][yi].resize((int) maxz+1);
    }
  }

  // Fill the array3D with zeros
  for (int xi=0; xi<maxx+1; xi++) {
    for (int yi=0; yi<maxy+1; yi++) {
      for (int zi=0; zi<maxz+1; zi++) {
        array3D[xi][yi][zi] = 0;
      }
    }
  }

  // Loop through all points
  for(int i=0; i<nrows; i++){

     myx = (double) pc(i,0);
     myy = (double) pc(i,1);
     myz = (double) pc(i,2);

     vx = (int) floor(myx);
     vy = (int) floor(myy);
     vz = (int) floor(myz);

     array3D[vx][vy][vz] = array3D[vx][vy][vz] + 1;

  }

  // Loop again over each point and iteratively shift the centroid of
  // neighbor points
  for(int i=0; i<nrows; i++){

    double meanx = (double) pc(i,0);
    double meany = (double) pc(i,1);
    double meanz = (double) pc(i,2);

    double oldx = -100.0;
    double oldy = -100.0;
    double oldz = -100.0;

    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;
    double sump = 0.0;

    double verticalweight = 0.0;
    double horizontalweight = 0.0;
    double weight = 0.0;

    // Keep iterating as long as the centroid is not reached
    int IterCounter = 0;
    do {

      sumx = 0.0;
      sumy = 0.0;
      sumz = 0.0;
      sump = 0.0;

      IterCounter = IterCounter + 1;
      // Calculate cylinder dimensions based on point height
      d = H2CW_fac * meanz;
      r = 0.5 * d;
      h = H2CL_fac * meanz;

      oldx = meanx;
      oldy = meany;
      oldz = meanz;

      // Calculate neighborhood margins to consider for this point
      boxminx = (int)floor(meanx - r);
      boxmaxx = (int)ceil(meanx + r);
      boxminy = (int)floor(meany - r);
      boxmaxy = (int)ceil(meany + r);
      boxminz = (int)floor(meanz - h/2);
      boxmaxz = (int)ceil(meanz + h/2);

      // Make sure the boxes don't exceed the global margins
      if (boxminx < minx) {boxminx = minx;}
      if (boxmaxx > maxx+1) {boxmaxx = maxx;}
      if (boxminy < miny) {boxminy = miny;}
      if (boxmaxy > maxy+1) {boxmaxy = maxy;}
      if (boxminz < minz) {boxminz = minz;}
      if (boxmaxz > maxz+1) {boxmaxz = maxz;}

      // Loop through all voxels in the neighborhood of the point
      for(int xi=boxminx; xi<boxmaxx; xi++){
        for(int yi=boxminy; yi<boxmaxy; yi++){
          for(int zi=boxminz; zi<boxmaxz; zi++){
            if(array3D[xi][yi][zi] > 0){
              if(InCylinder(xi, yi, zi, r, h, meanx, meany, meanz) == true){

                // If the option of a uniform kernel is set to false calculate the centroid
                // by multiplying all coodinates by their weights, depending on their relative
                // position within the cylinder, summing up the products and dividing by the
                // sum of all weights
                if(UniformKernel == false){
                  np = array3D[xi][yi][zi];
                  verticalweight = EpanechnikovFunction(h, meanz, zi);
                  horizontalweight = GaussFunction(d, meanx, meany, xi, yi);
                  weight = verticalweight * horizontalweight;
                  sumx = sumx + weight * np * xi;
                  sumy = sumy + weight * np * yi;
                  sumz = sumz + weight * np * zi;
                  sump = sump + weight * np;
                }
                // If the option of a uniform kernel is set to true calculate the centroid
                // by summing up all coodinates and dividing by the number of points
                else
                {
                  np = array3D[xi][yi][zi];
                  sumx = sumx + np * xi;
                  sumy = sumy + np * yi;
                  sumz = sumz + np * zi;
                  sump = sump + np * 1.0;
                }
              }
            }
          }
        }
      }

      meanx = sumx / sump;
      meany = sumy / sump;
      meanz = sumz / sump;

    }while(meanx != oldx && meany != oldy && meanz != oldz && IterCounter < MaxIter);

    centroidx[i] = meanx;
    centroidy[i] = meany;
    centroidz[i] = meanz;
  }

  return DataFrame::create(_["X"]= pc(_,0),_["Y"]= pc(_,1),_["Z"]= pc(_,2),_["CtrX"]= centroidx,_["CtrY"]= centroidy,_["CtrZ"]= centroidz);
}




