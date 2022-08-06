#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <inttypes.h>

const int32_t TRITET_TRUE = 1;
const int32_t TRITET_FALSE = 0;

const int32_t TRITET_SUCCESS = 0;

const int32_t TRITET_ERROR_NULL_DATA = 10;
const int32_t TRITET_ERROR_STRING_CONCAT = 20;

const int32_t TRITET_ERROR_NULL_POINT_LIST = 100;
const int32_t TRITET_ERROR_NULL_SEGMENT_LIST = 200;
const int32_t TRITET_ERROR_NULL_FACET_LIST = 300;
const int32_t TRITET_ERROR_NULL_FACET_POLYGON_LIST = 400;
const int32_t TRITET_ERROR_NULL_REGION_LIST = 500;
const int32_t TRITET_ERROR_NULL_HOLE_LIST = 600;

const int32_t TRITET_ERROR_INVALID_POINT_INDEX = 1000;
const int32_t TRITET_ERROR_INVALID_SEGMENT_INDEX = 2000;
const int32_t TRITET_ERROR_INVALID_SEGMENT_POINT_ID = 3000;
const int32_t TRITET_ERROR_INVALID_FACET_INDEX = 4000;
const int32_t TRITET_ERROR_INVALID_FACET_NUM_POLYGON = 5000;
const int32_t TRITET_ERROR_INVALID_FACET_POINT_INDEX = 6000;
const int32_t TRITET_ERROR_INVALID_FACET_POINT_ID = 7000;
const int32_t TRITET_ERROR_INVALID_REGION_INDEX = 8000;
const int32_t TRITET_ERROR_INVALID_HOLE_INDEX = 9000;

#endif  // CONSTANTS_H
