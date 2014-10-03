PHP_ARG_ENABLE(phpgeoextensions,
    [Whether to enable the "PHPGeoExtensions" extension],
    [  --enable-phpgeoextensions      Enable "PHPGeoExtensions" PHP extension support])

if test $PHP_PHPGEOEXTENSIONS != "no"; then
    PHP_REQUIRE_CXX()
    KMLSOURCES="kmlocal/src/KM_ANN.cpp kmlocal/src/KMeans.cpp kmlocal/src/KMterm.cpp \
		kmlocal/src/KMrand.cpp kmlocal/src/KCutil.cpp kmlocal/src/KCtree.cpp \
		kmlocal/src/KMdata.cpp kmlocal/src/KMcenters.cpp \
		kmlocal/src/KMfilterCenters.cpp kmlocal/src/KMlocal.cpp \
		clipperlib/clipper.cpp"
    PHP_ADD_INCLUDE("kmlocal/src/")
    PHP_ADD_INCLUDE("clipperlib/")
    PHP_ADD_LIBRARY(stdc++, 1, PHPGEOEXTENSIONS_SHARED_LIBADD)
    ALLSOURCES="php_geoextensions.cc kmpoints.cc $KMLSOURCES"
    PHP_NEW_EXTENSION(phpgeoextensions, $ALLSOURCES, $ext_shared)
    PHP_SUBST(PHPGEOEXTENSIONS_SHARED_LIBADD)
fi
