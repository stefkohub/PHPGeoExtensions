#include "kmpoints.h"
#include <algorithm>
#include <vector>
#include "Miniball.hpp"
#include "mercator.h"
#include "clipper.hpp"

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

static inline double mercatorDistance(double alat, double alng, double blat, double blng)
{
  double d1, d2;

  d1=(lat2y_m(alat)-lat2y_m(blat))*(lat2y_m(alat)-lat2y_m(blat));
  d2=(lon2x_m(alng)-lon2x_m(blng))*(lon2x_m(alng)-lon2x_m(blng));
  // printf("Distance squared: %f\n",d1+d2);
  return d1+d2;

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
  this->dataPtsID = (long*)calloc(this->maxPoints, sizeof(long));
  //this->dataPtsPtr=&dp;
  this->nPts=0;
  this->hullNum=0;
}

void kmpoints::newPoint(double lat, double lng, long id) {
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

#ifdef NOCLIPPER
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
#else 

vector <clusterPoint> kmpoints::getPolygon() {
  return getPolygon(0.0);
}

vector <clusterPoint> kmpoints::_expandPolygon(vector <clusterPoint> CH, float delta) {
    ClipperLib::Path subj;
    ClipperLib::Paths solution;

    for (int i = 0; i < CH.size(); i++) {
      subj << ClipperLib::IntPoint(floorf(lat2y_m(CH[i].lat)), floorf(lon2x_m(CH[i].lng)));
    }
    ClipperLib::ClipperOffset co;
    co.AddPath(subj, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);
    co.Execute(solution, delta);
    vector <clusterPoint> oCH;

    if (!solution.empty()) {
      for(std::vector<int>::size_type j = 0; j != solution[0].size(); j++) {
        clusterPoint cpoint;
        cpoint.lat=(y2lat_m(solution[0][j].X));
        cpoint.lng=(x2lon_m(solution[0][j].Y));
        oCH.push_back(cpoint);
      }  
    } 
    return oCH;
}

vector <clusterPoint> kmpoints::getPolygon(float delta) {

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
  if (delta) {
    cPoints=_expandPolygon(CH, delta);
  } else {
    cPoints=CH;
  }

  this->hullNum=1;

  return cPoints;

}
#endif

#ifdef VECTOR
circlePoint kmpoints::getCircle() {

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
    p[0]=lat2y_m(dp[i][0]);
    p[1]=lon2x_m(dp[i][1]);
    lp.push_back(p);
  }

  MB mb (2, lp.begin(), lp.end());

  const KMcoord *Ccenter = mb.center();
  circlePoint retVal;
  retVal.lat=y2lat_m(Ccenter[0]);
  retVal.lng=x2lon_m(Ccenter[1]);
  retVal.radius=sqrt(mb.squared_radius());
  // retVal.radius=mb.squared_radius();

  if (retVal.radius<RESOLUTION) retVal.radius=0;

  // number of support points
  *kmOut << "Number of support points:\n  ";
  *kmOut << mb.nr_support_points() << endl;

  // support points on the boundary determine the smallest enclosing ball
  *kmOut << "Support point indices (numbers refer to the input order):\n  ";
  MB::SupportPointIterator it = mb.support_points_begin();
  PointIterator first = lp.begin();
  for (; it != mb.support_points_end(); ++it) {
    *kmOut << distance(first, *it) << " "; // 0 = first point
  }
  *kmOut << endl;
  
  // relative error: by how much does the ball fail to contain all points? 
  //                 tiny positive numbers come from roundoff and are ok
  *kmOut << "Relative error:\n  ";
  KMcoord suboptimality;
  *kmOut << mb.relative_error (suboptimality) <<  endl;
  
  // suboptimality: by how much does the ball fail to be the smallest
  //                enclosing ball of its support points? should be 0 
  //                in most cases, but tiny positive numbers are again ok
  *kmOut << "Suboptimality:\n  ";
  *kmOut << suboptimality <<  endl;

  // validity: the ball is considered valid if the relative error is tiny
  //           (<= 10 times the machine epsilon) and the suboptimality is zero
  *kmOut << "Validity:\n  ";
  *kmOut << (mb.is_valid() ? "ok" : "possibly invalid") << endl;

  // computation time
  *kmOut << "Computation time was "<< mb.get_time() << " seconds\n";



  return retVal;
}
#else
circlePoint kmpoints::getCircle() {
  
  // *kmOut << "ALGORITMO SENZA ITERATOR" <<endl;

  // generate random points and store them in a 2-d array
  // ----------------------------------------------------
  int n = this->getNumPts();
  KMcoord** ap = new KMcoord*[n];
  
  KMdataArray dp = (this->dataPts)->getPts();

  for (int i=0; i<n; ++i) {
    KMcoord* p = new KMcoord[2];
    p[0]=lat2y_m(dp[i][0]);
    p[1]=lon2x_m(dp[i][1]);
    ap[i]=p;
  }

  // define the types of iterators through the points and their coordinates
  // ----------------------------------------------------------------------
  typedef KMcoord* const* PointIterator; 
  typedef const KMcoord* CoordIterator;

  // create an instance of Miniball
  // ------------------------------
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
    MB;
  MB mb (2, ap, ap+n);
  
  // output results
  // --------------
  // center
  const KMcoord* center = mb.center(); 

  circlePoint retVal;
  retVal.lat=y2lat_m(center[0]);
  retVal.lng=x2lon_m(center[1]);
  retVal.radius=sqrt(mb.squared_radius());
  return retVal;
}

#endif

int kmpoints::getNumIntersects(double lat, double lng, double criteria) {
  KMdataArray dp = (this->dataPts)->getPts();
  int nIntersects=0;

  for (int i=0;i<this->getNumPts();i++) {
    if (cosinesLaw(dp[i][0],dp[i][1],lat,lng)<=criteria)
      nIntersects++;
  }
  return nIntersects;
}

vector <long> kmpoints::getIdIntersects(double lat, double lng, double criteria) {
  KMdataArray dp = (this->dataPts)->getPts();
  KMpoint thePoint = kmAllocPt(2);
  vector <long> theIntersects;

  thePoint[0]=(KMcoord)lat;
  thePoint[1]=(KMcoord)lng;

  for (int i=0;i<this->getNumPts();i++) {
    if (cosinesLaw(dp[i][0],dp[i][1],lat,lng)<=criteria)
      theIntersects.push_back(this->dataPtsID[i]);
  }
  return theIntersects;
}

vector <clusterPoint> kmpoints::getOffendingPts() {
  vector <clusterPoint> oPts;
  return oPts;
  // KMdataArray dp = (this->dataPts)->getPts();
  // return (this->dataPts)->getoPts();
}


int kmpoints::getNumPts() {
  return (this->dataPts)->getNPts();
}

int kmpoints::getHullNum() {
  return this->hullNum;
}

KMdataArray kmpoints::getDataPts() {
  return (this->dataPts)->getPts();
}

void kmpoints::createDistanceMatrix() {
  KMdataArray dp = this->getDataPts();
  // vector< vectorDistance > myvec(this->getNumPts(), vectorDistance(this->getNumPts()));
  this->dMatrix.resize(this->getNumPts(), vectorDistance(this->getNumPts()));

  for (int i=0;i<this->getNumPts();i++) {
    for (int j=i+1;j<this->getNumPts();j++) {
      this->dMatrix[i][j]=cosinesLaw(dp[i][0],dp[i][1],dp[j][0],dp[j][1]);
    }
  }
}

