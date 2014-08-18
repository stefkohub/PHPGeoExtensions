#ifndef PHP_GEOEXTENSIONS_H
#define PHP_GEOEXTENSIONS_H

#define PHP_PHPGEOEXTENSIONS_EXTNAME  "geoextensions"
#define PHP_PHPGEOEXTENSIONS_EXTVER   "0.5"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

extern "C" {
#include "php.h"
}

extern zend_module_entry geoextensions_module_entry;
#define phpext_geoextensions_ptr &geoextensions_module_entry;

#endif /* PHP_GEOEXTENSIONS_H */
