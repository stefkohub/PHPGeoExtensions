#include <iostream>
#include "clipper.hpp"  
#include <math.h>

using namespace ClipperLib;
#define PRECISION 1000000000000.0

/*
var paths = [
  new google.maps.LatLng(40.7969256,8.5729409),
  new google.maps.LatLng(40.804107,8.571548),
  new google.maps.LatLng(40.8368187,8.40098),
  new google.maps.LatLng(40.835903,8.887521)
];


void grow() {
  KMdataArray dp = (this->dataPts)->getPts();

  vector <clusterPoint> cPoints(vectorPoint((this->dataPts)->getNPts()));

  vector <clusterPoint> CH;

  Path subj;
  Paths solution;

  for (int i = 0; i < (this->dataPts)->getNPts(); i++) {
    subj << IntPoint(dp[i][0]*PRECISION, dp[i][1]*PRECISION);
  }
  ClipperOffset co;
  co.AddPath(subj, jtRound, etClosedPolygon);
  co.Execute(solution, 5.0);
  for (int i = solution.begin(); i != solution.end(); i++) {
    clusterPoint cpoint;
    cpoint.lat=floorf(solution[i][0])/PRECISION;
    cpoint.lng=floorf(solution[i][1])/PRECISION;
    CH.push_back(cpoint);
  }
}
*/

int main()
{
  Path subj;
  Paths solution;

  float coord[4][2];

  coord[0][0]=40.7969256;
  coord[0][1]=8.5729409;
  coord[1][0]=40.804107;
  coord[1][1]=8.571548;
  coord[2][0]=40.8368187;
  coord[2][1]=8.40098;
  coord[3][0]=40.835903;
  coord[3][1]=8.887521;


  for (int i=0;i<4;i++) {
    subj << IntPoint(coord[i][0]*PRECISION, coord[i][1]*PRECISION);
  }
/*
  subj <<
    IntPoint(348,257) << IntPoint(364,148) << IntPoint(362,148) <<
    IntPoint(326,241) << IntPoint(295,219) << IntPoint(258,88) <<
    IntPoint(440,129) << IntPoint(370,196) << IntPoint(372,275);*/
  ClipperOffset co;
  co.AddPath(subj, jtMiter, etClosedPolygon);
  co.Execute(solution, -7.0);

  std::cout << solution <<std::endl;
  //for(std::vector<int>::size_type i = 0; i != solution.size(); i++) {
  int i=0;
    std::cout << solution[i] << " - ";
    for(std::vector<int>::size_type j = 0; j != solution[i].size(); j++) {
      std::cout << ">>> " << floorf(solution[i][j].X)/PRECISION << "--" << floorf(solution[i][j].Y)/PRECISION <<std::endl;
    }
  //}
  std::cout << "Fine" <<std::endl;
   
  //draw solution ...
  // DrawPolygons(solution, 0x4000FF00, 0xFF009900);
} 
