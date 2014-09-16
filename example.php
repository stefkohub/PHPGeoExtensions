<?php

include "minbound-altro.inc";

$MAXP=count($OrigPoints);
$CENTERS=(int)sqrt($MAXP/2);

$k = new kmpoints($MAXP);

for ($n=0;$n<$MAXP;$n++) {
  // 		latitude, longitude, point id 
  $k->newPoint($OrigPoints[$n][0],$OrigPoints[$n][1],$n);
}

print "The minimum circle containing all the points is: \n";
print_r($k->getCircle());

$lat= 45.43812;
$lng = 12.31816;
print "Looking for intersections from center: $lat, $lng with radius 5Km\n";

print "Number of intersections: ".$k->getNumIntersects($lat, $lng, 5)."\n";
print "Array of IDs intersecating: \n";
print_r($k->getIdIntersects($lat, $lng, 5));

print "The polygon (convex hull) containing all the points:\n";
$hullRes = $k->getPolygon();
print_r($hullRes);

print "The array of polygons (convex hull) containing all the points:\n";
print "Using kmeans with $CENTERS centers (approx formula: square root of half total points)\n";
$hullRes = $k->getPolygons($CENTERS);
print_r($hullRes);

print "The array of clusters of all the points:\n";
print "Using kmeans with $CENTERS centers (approx formula: square root of half total points)\n";
$theClusters = $k->getClusters($CENTERS);
print_r($theClusters);

