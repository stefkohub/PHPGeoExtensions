PHPGeoExtensions
================

A PHP estension to add some fast algorithm for spatial clustering, partitioning and shapes (i.e. convex hull, minimum circle, k-means clustering, etc)

**WARNING**
The kmlocal subdirectory is the 1.7.2 version of kmlocal from: http://www.cs.umd.edu/~mount/Projects/KMeans/

I renamed kmlocal-1.7.2 directory to kmlocal to overcome some parsing issues in config.m4 (filename is defined as the first token before the first point)

Example Usage in PHP
====================

Initialize kmpoints class that will use maximum 1000 points:

```PHP
$k = new kmpoints(1000);
```

Add some datapoints with latitude, longitude, and an optional numerical identifier:

```PHP
$k->newPoint(41.9100711,12.5359979,1);
$k->newPoint(41.7872566,12.243663,2);
```

Get the minimum circle containing those points:

```PHP
print_r($k->getCircle());
```

Get the polygon (convex hull):

```PHP
print_r($k->getPolygon());
```

Get an array of convex hull obtained from a first kmeans clustering step:

```PHP
$theKparameter=1;
print_r($k->getPolygons($theKparameter));
```

Credits
================

|Package|URL|
|------------|----------------------------------------------------------------------------|
|Miniball | http://www.inf.ethz.ch/personal/gaertner/miniball.html|
|KMLocal | http://www.cs.umd.edu/~mount/Projects/KMeans/|
|Clipperlib | http://www.angusj.com/delphi/clipper.php|
|Mercator Transformations | http://wiki.openstreetmap.org/wiki/Mercator#Spherical_Mercator|

