#pragma once

#include <cuda_fp16.h>

#define IND_TYPE    uint32_t
#define CON_TYPE    __half
#define TIMER_TYPE  uint8_t
#define INF_TYPE    uint32_t
#define MASK_TYPE   uint8_t
#define ATTR_TYPE   uint16_t

#define STATUS_SUSCEPTIBLE (0)
#define STATUS_EXPOSED     (1)
#define STATUS_INFECTIOUS  (2)
#define STATUS_RECOVERED   (3)

struct TxEvent {
    int32_t    source;
    int32_t    time_exposed;
    int32_t    time_infectious;
    int32_t    time_recovered;
};
