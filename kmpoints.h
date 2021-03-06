#ifndef KMLOCAL_POINTS_H
#define KMLOCAL_POINTS_H

#include "KMlocal.h"                    // k-means algorithms
#include <vector>


#define PHP_CLUSTERPOINTS_RES_NAME "kmpoints::clusterPoints"

typedef struct _clusterPoint{
  int idp;
  KMcoord lat;
  KMcoord lng;
  bool operator <(const _clusterPoint &p) const {
    return lat < p.lat || (lat == p.lat && lng < p.lng);
  }
} clusterPoint;

typedef struct _circlePoint {
  KMcoord lat;
  KMcoord lng;
  double radius;
} circlePoint;

typedef std::vector<clusterPoint> vectorPoint;
typedef std::vector<double> vectorDistance;


class kmpoints {
public:
    kmpoints(int maxPoints);
    void newPoint(double lat, double lng);
    void newPoint(double lat, double lng, long id);
    KMcoord getLat(int);
    KMcoord getLng(int);
    vector < vectorPoint > getPolygons(int);
    vector < vectorPoint > getClusters(int);
    vector <clusterPoint> getPolygon();
    vector <clusterPoint> getPolygon(float);
    circlePoint getCircle();
    int getNumPts();
    int getHullNum();
    int getNumIntersects(double, double, double);
    int getNumIntersectsv2(double, double, double);
    vector <long> getIdIntersects(double, double, double);
    vector <clusterPoint> getOffendingPts();
    KMdataArray getDataPts();
    vector <vectorDistance> dMatrix;

private:
    //KMcoord lat;
    //KMcoord lng;
    int nPts;
    KMdata *dataPts;
    long *dataPtsID;
    int theDim;
    int maxPoints;
    int hullNum;
    vector <clusterPoint> _expandPolygon(vector <clusterPoint>, float);
    void createDistanceMatrix();
};

#endif /* KMLOCAL_POINTS_H */
