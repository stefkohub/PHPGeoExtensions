#include "php_geoextensions.h"
// #include "KMlocal.h"
#include "kmpoints.h"

#include <map>

zend_object_handlers kmpoints_object_handlers;
zend_class_entry *kmpoints_ce;
int le_clusterPoints;

struct kmpoints_object {
    zend_object std;
    kmpoints *kmpoint;
};

void kmpoints_free_storage(void *object TSRMLS_DC)
{
    kmpoints_object *obj = (kmpoints_object *)object;
    delete obj->kmpoint; 

    zend_hash_destroy(obj->std.properties);
    FREE_HASHTABLE(obj->std.properties);

    efree(obj);
}

zend_object_value kmpoints_create_handler(zend_class_entry *type TSRMLS_DC)
{
    zval *tmp;
    zend_object_value retval;

    kmpoints_object *obj = (kmpoints_object *)emalloc(sizeof(kmpoints_object));
    memset(obj, 0, sizeof(kmpoints_object));
    obj->std.ce = type;

    ALLOC_HASHTABLE(obj->std.properties);
    zend_hash_init(obj->std.properties, 0, NULL, ZVAL_PTR_DTOR, 0);
    zend_hash_copy(obj->std.properties, &type->default_properties,
        (copy_ctor_func_t)zval_add_ref, (void *)&tmp, sizeof(zval *));

    retval.handle = zend_objects_store_put(obj, NULL, kmpoints_free_storage, NULL TSRMLS_CC);
    retval.handlers = &kmpoints_object_handlers;

    return retval;
}

PHP_METHOD(kmpoints, __construct)
{
    long maxp;
    kmpoints *kmp = NULL;
    zval *object = getThis();

    if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "l", &maxp) == FAILURE) {
        RETURN_NULL();
    }

    kmp = new kmpoints((int)maxp);
    kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(object TSRMLS_CC);
    obj->kmpoint = kmp;
}

PHP_METHOD(kmpoints, getLng)
{
    kmpoints *kmp;
    long theIndex;

    kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
    kmp = obj->kmpoint;
    if (kmp != NULL) {
      if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "l", &theIndex) == FAILURE) {
        RETURN_NULL();
      }
      // php_printf("Index: %d\n",(int)kmp->getLng((int)theIndex));
      RETURN_DOUBLE((kmp->getLng((int)theIndex)));
    }

}

PHP_METHOD(kmpoints, getLat)
{
    kmpoints *kmp;
    long theIndex;

    kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
    kmp = obj->kmpoint;
    if (kmp != NULL) {
      if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "l", &theIndex) == FAILURE) {
        RETURN_NULL();
      }
      RETURN_DOUBLE((kmp->getLat((int)theIndex)));
    }

}

PHP_METHOD(kmpoints, newPoint)
{
  double lat;
  double lng;
  double uid=-1;
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  
  if (kmp != NULL) {
    if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "dd|d", &lat, &lng, &uid) == FAILURE) {
      RETURN_NULL();
    }
    //php_printf("lat, lng: %f, %f\n",lat,lng);
    if (uid==-1) {
      kmp->newPoint(lat, lng);
    } else {
      kmp->newPoint(lat, lng, uid);
    }
    RETURN_NULL();
  }
  
}

PHP_METHOD(kmpoints, getPolygons)
{
  long nCenter;
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  vector < vectorPoint > retval;
  char hullName[100];
  zval	*theHull;
  zval  *thePoint;

  if (kmp != NULL) {

      if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "l", &nCenter) == FAILURE) {
          RETURN_NULL();
      }

      // php_printf("nPts=%d\n",(int)nCenter);

      retval = kmp->getPolygons((int)nCenter);

      array_init(return_value);

      for(int c=0;c<kmp->getHullNum();c++) {
        ALLOC_INIT_ZVAL(theHull);
	array_init(theHull);
	sprintf(hullName, "%d", c);
        for (int i=0;i<retval[c].size();i++) {
          ALLOC_INIT_ZVAL(thePoint);
          // php_printf("Hull %d-%d: (%f,%f)\n",c,i,retval[c][i].lat,retval[c][i].lng);
	  array_init(thePoint);
	  add_next_index_double(thePoint, retval[c][i].lat);
	  add_next_index_double(thePoint, retval[c][i].lng);
	  add_next_index_zval(theHull, thePoint);
        }
	add_assoc_zval(return_value, hullName, theHull);
      }

  }
  
}

PHP_METHOD(kmpoints, getPolygon)
{
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  vector < clusterPoint > retval;
  zval	*theHull;
  zval  *thePoint;

  if (kmp != NULL) {

      retval = kmp->getPolygon();
      array_init(return_value);

        ALLOC_INIT_ZVAL(theHull);
	array_init(theHull);
        for (int i=0;i<retval.size();i++) {
          ALLOC_INIT_ZVAL(thePoint);
          // php_printf("Hull (%d): (%f,%f)\n",i,retval[i].lat,retval[i].lng);
	  array_init(thePoint);
	  add_next_index_double(thePoint, retval[i].lat);
	  add_next_index_double(thePoint, retval[i].lng);
	  // add_next_index_zval(theHull, thePoint);
	  add_next_index_zval(return_value, thePoint);
        }
	//add_assoc_zval(return_value, "hull", theHull);

  }
  
}

PHP_METHOD(kmpoints, getCircle)
{
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  circlePoint retval;
  zval	*theCenter;
  zval  *theRadius;

  if (kmp != NULL) {

      retval = kmp->getCircle();
      array_init(return_value);

        ALLOC_INIT_ZVAL(theCenter);
	array_init(theCenter);
        add_next_index_double(theCenter, retval.lat);
        add_next_index_double(theCenter, retval.lng);
	add_assoc_zval(return_value, "center", theCenter);
        ALLOC_INIT_ZVAL(theRadius);
        array_init(theRadius);
	add_next_index_double(theRadius, retval.radius);
	add_assoc_zval(return_value, "radius", theRadius);

  }
  
}

PHP_METHOD(kmpoints, getNumPts)
{
    kmpoints *kmp;

    kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
    kmp = obj->kmpoint;
    if (kmp != NULL) {
      RETURN_DOUBLE((kmp->getNumPts()));
    }

}

PHP_METHOD(kmpoints, getHullNum)
{
    kmpoints *kmp;

    kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
    kmp = obj->kmpoint;
    if (kmp != NULL) {
      RETURN_LONG((kmp->getHullNum()));
    }

}

PHP_METHOD(kmpoints, getNumIntersects)
{
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  double lat;
  double lng;
  double radius;
  long retval;

  if (kmp != NULL) {

      if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddd", &lat, &lng, &radius) == FAILURE) {
          RETURN_NULL();
      }

      retval = kmp->getNumIntersects(lat, lng, radius);
   
      RETURN_LONG(retval);

  }
  
}

PHP_METHOD(kmpoints, getIdIntersects)
{
  kmpoints *kmp;
  kmpoints_object *obj = (kmpoints_object *)zend_object_store_get_object(getThis() TSRMLS_CC);
  kmp = obj->kmpoint;
  double lat;
  double lng;
  double radius;
  vector <int> retval;

  if (kmp != NULL) {

      if (zend_parse_parameters(ZEND_NUM_ARGS() TSRMLS_CC, "ddd", &lat, &lng, &radius) == FAILURE) {
          RETURN_NULL();
      }

      retval = kmp->getIdIntersects(lat, lng, radius);
      array_init(return_value);

      for (int i=0;i<retval.size();i++) {
	  add_next_index_long(return_value, retval[i]);
      }

  }
  
}

function_entry kmpoints_methods[] = {
    PHP_ME(kmpoints,  __construct,     		NULL, ZEND_ACC_PUBLIC | ZEND_ACC_CTOR)
    PHP_ME(kmpoints,  getLng,      		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getLat,           	NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  newPoint,           	NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getPolygons, 		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getPolygon, 		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getCircle, 		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getNumPts, 		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getHullNum, 		NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getNumIntersects, 	NULL, ZEND_ACC_PUBLIC)
    PHP_ME(kmpoints,  getIdIntersects,	 	NULL, ZEND_ACC_PUBLIC)
    {NULL, NULL, NULL}
};

PHP_MINIT_FUNCTION(phpgeoextensions)
{
    zend_class_entry ce;
    INIT_CLASS_ENTRY(ce, "kmpoints", kmpoints_methods);
    kmpoints_ce = zend_register_internal_class(&ce TSRMLS_CC);
    kmpoints_ce->create_object = kmpoints_create_handler;
    memcpy(&kmpoints_object_handlers,
        zend_get_std_object_handlers(), sizeof(zend_object_handlers));
    kmpoints_object_handlers.clone_obj = NULL;

    le_clusterPoints = zend_register_list_destructors_ex(NULL, NULL, PHP_CLUSTERPOINTS_RES_NAME, module_number);

    return SUCCESS;
}


zend_module_entry phpgeoextensions_module_entry = {
#if ZEND_MODULE_API_NO >= 20010901
    STANDARD_MODULE_HEADER,
#endif
    PHP_PHPGEOEXTENSIONS_EXTNAME,
    NULL,        /* Functions */
    PHP_MINIT(phpgeoextensions),        /* MINIT */
    NULL,        /* MSHUTDOWN */
    NULL,        /* RINIT */
    NULL,        /* RSHUTDOWN */
    NULL,        /* MINFO */
#if ZEND_MODULE_API_NO >= 20010901
    PHP_PHPGEOEXTENSIONS_EXTVER,
#endif
    STANDARD_MODULE_PROPERTIES
};

#ifdef COMPILE_DL_PHPGEOEXTENSIONS
extern "C" {
ZEND_GET_MODULE(phpgeoextensions)
}
#endif

