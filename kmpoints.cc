#include "kmpoints.h"
#include <algorithm>
#include <vector>
#include "Miniball.hpp"

#define RESOLUTION 0.1
#define R 6371.0
#define R2 40589641.0
#define PI (3.141592653589793)

clusterPoint center;

KMterm  term(100, 0, 0, 0,              // run for 100 stages
                0.10,                   // min consec RDL
                0.10,                   // min accum RDL
                3,                      // max run stages
                0.50,                   // init. prob. of acceptance
                10,                     // temp. run length
                0.95);                  // temp. reduction factor

// Private helper functions

clusterPoint getCenterPoint(KMdataArray pts, int pointNum) {
  clusterPoint thecenter;
  // int pointNum = npts;
  for (int i=0;i<pointNum; i++) {
    thecenter.lat+=pts[i][0];
    thecenter.lng+=pts[i][1];
  }
  thecenter.lat/=pointNum;
  thecenter.lng/=pointNum;
  return thecenter;
}

bool coordinateOrder(const clusterPoint &a, const clusterPoint &b)
{
    if (a.lat - center.lat >= 0 && b.lat - center.lat < 0)
        return true;
    if (a.lat - center.lat < 0 && b.lat - center.lat >= 0)
        return false;
    if (a.lat - center.lat == 0 && b.lat - center.lat == 0) {
        if (a.lng - center.lng >= 0 || b.lng - center.lng >= 0)
            return a.lng > b.lng;
        return b.lng > a.lng;
    }

    // compute the cross product of vectors (center -> a) x (center -> b)
    int det = (a.lat - center.lat) * (b.lng - center.lng) - (b.lat - center.lat) * (a.lng - center.lng);
    if (det < 0)
        return true;
    if (det > 0)
        return false;

    // points a and b are on the same line from the center
    // check which point is closer to the center
    int d1 = (a.lat - center.lat) * (a.lat - center.lat) + (a.lng - center.lng) * (a.lng - center.lng);
    int d2 = (b.lat - center.lat) * (b.lat - center.lat) + (b.lng - center.lng) * (b.lng - center.lng);
    return d1 > d2;
}

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
double cross(const clusterPoint &O, const clusterPoint &A, const clusterPoint &B)
{
	return (A.lat - O.lat) * (B.lng - O.lng) - (A.lng - O.lng) * (B.lat - O.lat);
}

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
vector<clusterPoint> convex_hull(vector<clusterPoint> P)
{
	int n = P.size(), k = 0;
	vector<clusterPoint> H(2*n);

	// Sort points lexicographically
	sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k);
	return H;
}

static inline double equirectangularApproxDist(KMpoint pa, KMpoint pb)
{
    KMpoint a = kmAllocPt(2);
    KMpoint b = kmAllocPt(2);
    a[0]=pa[0]/(180.0*PI);
    a[1]=pa[1]/(180.0*PI);
    b[0]=pb[0]/(180.0*PI);
    b[1]=pb[1]/(180.0*PI);

    double x = b[1]-a[1]*cos((a[0]+b[0])/2);
    double y = b[0]-a[0];
    // returning square distance
    double dist=R2*((x*x)+(y*y));
    if (dist<RESOLUTION)
      return 0.0;
    else
      return dist;

}

// LINEAR VERSION
static inline double cosinesLaw(double alat, double alng, double blat, double blng)
{
    double d1, d2;
    d1=(alat - blat)/180.0*PI*R;
    if (d1<0.0)
      d1=-d1;
    d2=(alng - blng)/180.0*PI*R; 
    if (d2<0.0)
      d2=-d2;
    return d1 > d2 ? d1 : d2;
}

static inline double cosinesLaw(KMpoint a, KMpoint b)
{
    double d1, d2;
    d1=(a[0] - b[0])/180.0*PI*R;
    if (d1<0.0)
      d1=-d1;
    d2=(a[1] - b[1])/180.0*PI*R; //bisognerebbe correggere il raggio terrestre alla data[0]itudine
    if (d2<0.0)
      d2=-d2;
    //printf("d1 e d2 tra %f,%f e %f,%f è: %f\n",a[0],a[1],b[0],b[1],d1,d2);
    return d1 > d2 ? d1 : d2;
}

// LINEAR VERSION
static inline double cosinesLaw(clusterPoint a, clusterPoint b)
{
    double d1, d2;
    d1=(a.lat - b.lat)/180.0*PI*R;
    if (d1<0.0)
      d1=-d1;
    d2=(a.lng - b.lng)/180.0*PI*R; //bisognerebbe correggere il raggio terrestre alla data latitudine
    if (d2<0.0)
      d2=-d2;
    //printf("d1 e d2 tra %f,%f e %f,%f è: %f\n",a.lat,a.lng,b.lat,b.lng,d1,d2);
    return d1 > d2 ? d1 : d2;
}

// Class functions

kmpoints::kmpoints(int maxPoints) {
  // KMdata dataPts(dim, nPts);
  this->theDim=2;
  this->maxPoints=maxPoints;
  term.setAbsMaxTotStage(200);
  // KMdata dp(this->theDim, this->maxPoints);
  this->dataPts = new KMdata(this->theDim, this->maxPoints);
  this->dataPtsID = (int*)calloc(this->maxPoints, sizeof(int));
  //this->dataPtsPtr=&dp;
  this->nPts=0;
  this->hullNum=0;
}

void kmpoints::newPoint(double lat, double lng, int id) {
  KMpointArray pa = (this->dataPts)->getPts();
  pa[this->nPts][0] = (KMcoord) lat;
  pa[this->nPts][1] = (KMcoord) lng;
  this->dataPtsID[this->nPts]=id;
  this->nPts++;
  (this->dataPts)->setNPts(this->nPts);
}

void kmpoints::newPoint(double lat, double lng) {
  KMpointArray pa = (this->dataPts)->getPts();
  pa[this->nPts][0] = (KMcoord) lat;
  pa[this->nPts][1] = (KMcoord) lng;
  this->nPts++;
  (this->dataPts)->setNPts(this->nPts);
}

KMcoord kmpoints::getLat(int index) {
  // KMdata dp=*(this->dataPts);
  KMdataArray dp = (this->dataPts)->getPts();
  KMcoord c = dp[index][0];
  return c;
}

KMcoord kmpoints::getLng(int index) {
  // KMdata dp=*(this->dataPts);
  KMdataArray dp = (this->dataPts)->getPts();
  KMcoord c = dp[index][1];
  return c;
}

vector < vectorPoint > kmpoints::getPolygons(int k) {
  KMdataArray dp = (this->dataPts)->getPts();

  (this->dataPts)->buildKcTree();
  KMfilterCenters ctrs(k, *(this->dataPts));

  // KMlocalLloyds       kmAlg(ctrs, term);
  // ctrs = kmAlg.execute();
  KMlocalHybrid Local_Hybrid(ctrs, term);   // EZ-Hybrid heuristic
  ctrs = Local_Hybrid.execute();

  KMctrIdxArray closeCtr = new KMctrIdx[(this->dataPts)->getNPts()];
  double* sqDist = new double[(this->dataPts)->getNPts()];
  ctrs.getAssignments(closeCtr, sqDist);

  vector < vectorPoint > cPoints(ctrs.getK(), vectorPoint((this->dataPts)->getNPts()));

  for (int c=0;c<ctrs.getK(); c++) {
    int p=0;
    vectorPoint innerPoints; //((this->dataPts)->getNPts());
    for (int i = 0; i < (this->dataPts)->getNPts(); i++) {
      if (c==closeCtr[i]) {
        clusterPoint cpoint;
        cpoint.lat=dp[i][0];
        cpoint.lng=dp[i][1];
        innerPoints.push_back(cpoint);
        p++;
      }
    }
    center.lat=ctrs[c][0];
    center.lng=ctrs[c][1];
    std::sort(innerPoints.begin(), innerPoints.end(), coordinateOrder);
    cPoints[c]=innerPoints;
    innerPoints.clear();
  }

  vector <clusterPoint> CH;
  for (int i = 0; i < ctrs.getK(); i++) {
        CH=convex_hull(cPoints[i]);
	cPoints[i]=CH;
  }
  this->hullNum=ctrs.getK();

  delete [] closeCtr;
  delete [] sqDist;

  return cPoints;

}

vector < vectorPoint > kmpoints::getClusters(int k) {
  KMdataArray dp = (this->dataPts)->getPts();

  (this->dataPts)->buildKcTree();
  KMfilterCenters ctrs(k, *(this->dataPts));

  // KMlocalLloyds       kmAlg(ctrs, term);
  // ctrs = kmAlg.execute();
  KMlocalHybrid Local_Hybrid(ctrs, term);   // EZ-Hybrid heuristic
  ctrs = Local_Hybrid.execute();

  KMctrIdxArray closeCtr = new KMctrIdx[(this->dataPts)->getNPts()];
  double* sqDist = new double[(this->dataPts)->getNPts()];
  ctrs.getAssignments(closeCtr, sqDist);

  vector < vectorPoint > cPoints(ctrs.getK(), vectorPoint((this->dataPts)->getNPts()));

  for (int c=0;c<ctrs.getK(); c++) {
    int p=0;
    vectorPoint innerPoints; //((this->dataPts)->getNPts());
    for (int i = 0; i < (this->dataPts)->getNPts(); i++) {
      if (c==closeCtr[i]) {
        clusterPoint cpoint;
        cpoint.lat=dp[i][0];
        cpoint.lng=dp[i][1];
        innerPoints.push_back(cpoint);
        p++;
      }
    }
    cPoints[c]=innerPoints;
    innerPoints.clear();
  }

  delete [] closeCtr;
  delete [] sqDist;

  return cPoints;

}


vector <clusterPoint> kmpoints::getPolygon() {
  KMdataArray dp = (this->dataPts)->getPts();

  vector <clusterPoint> cPoints(vectorPoint((this->dataPts)->getNPts()));

  vectorPoint innerPoints;
  for (int i = 0; i < (this->dataPts)->getNPts(); i++) {
    clusterPoint cpoint;
    cpoint.lat=dp[i][0];
    cpoint.lng=dp[i][1];
    innerPoints.push_back(cpoint);
  }
  center = getCenterPoint(dp, (this->dataPts)->getNPts());
  std::sort(innerPoints.begin(), innerPoints.end(), coordinateOrder);
  cPoints=innerPoints;

  vector <clusterPoint> CH;
  CH=convex_hull(cPoints);
  cPoints=CH;
  this->hullNum=1;

  return cPoints;

}

circlePoint kmpoints::getCircle() {
  //vector <clusterPoint> cPoints;

  list<vector<KMcoord> > lp;
  // define the types of iterators through the points and their coordinates
  // ----------------------------------------------------------------------
  typedef std::list<std::vector<KMcoord> >::const_iterator PointIterator;
  typedef std::vector<KMcoord>::const_iterator CoordIterator;

  // create an instance of Miniball
  // ------------------------------
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> >
    MB;

  KMdataArray dp = (this->dataPts)->getPts();
  for (int i=0;i<this->getNumPts();i++) {
    vector <KMcoord> p(2);
    p[0]=dp[i][0];
    p[1]=dp[i][1];
    lp.push_back(p);
  }

  MB mb (2, lp.begin(), lp.end());

  const KMcoord *Ccenter = mb.center();
  circlePoint retVal;
  retVal.lat=Ccenter[0];
  retVal.lng=Ccenter[1];
  retVal.radius=sqrt(mb.squared_radius());
  if (retVal.radius<RESOLUTION) retVal.radius=0;

  return retVal;
}

int kmpoints::getNumIntersects(double lat, double lng, double criteria) {
  KMdataArray dp = (this->dataPts)->getPts();
  int nIntersects=0;

  criteria*=criteria;

  for (int i=0;i<this->getNumPts();i++) {
    if (cosinesLaw(dp[i][0],dp[i][1],lat,lng)<=criteria)
      nIntersects++;
  }
  return nIntersects;
}

vector <int> kmpoints::getIdIntersects(double lat, double lng, double criteria) {
  KMdataArray dp = (this->dataPts)->getPts();
  KMpoint thePoint = kmAllocPt(2);
  vector <int> theIntersects;

  thePoint[0]=(KMcoord)lat;
  thePoint[1]=(KMcoord)lng;

  criteria*=criteria;

  for (int i=0;i<this->getNumPts();i++) {
    if (cosinesLaw(dp[i],thePoint)<=criteria)
      theIntersects.push_back(this->dataPtsID[i]);
  }
  return theIntersects;
}


int kmpoints::getNumPts() {
  return (this->dataPts)->getNPts();
}

int kmpoints::getHullNum() {
  return this->hullNum;
}
