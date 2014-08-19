PHPGeoExtensions
================

A PHP estension to add some fast algorithm for spatial clustering, partitioning and shapes (i.e. convex hull, minimum circle, k-means clustering, etc)

**WARNING**
The kmlocal subdirectory is the 1.7.2 version of kmlocal from: http://www.cs.umd.edu/~mount/Projects/KMeans/

I renamed kmlocal-1.7.2 directory to kmlocal to overcome some parsing issues in config.m4 (filename is defined as the first token before the first point)

Example Usage in PHP
====================

Initialize kmponts class that will use maximum 1000 points:

```
$k = new kmpoints(1000);
```

Add some datapoints with latitude, longitude, and an optional numerical identifier:

```
$k->newPoint(41.9100711,12.5359979,1);
$k->newPoint(41.7872566,12.243663,2);
```

Get the minimum circle containing those points:

```
print_r($k->getCircle());
```

Get the polygon (convex hull):

```
print_r($k->getPolygon());
```

Get an array of convex hull obtained from a first kmeans clustering pass:

```
$theKparameter=1;
print_r($k->getPolygons($theKparameter));
```

