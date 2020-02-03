#pragma once

#include <cuda_fp16.h>

#define IND_TYPE    uint32_t
#define CON_TYPE    __half
#define TIMER_TYPE  uint8_t
#define INF_TYPE    uint32_t
#define MASK_TYPE   uint8_t

#define STATUS_SUSCEPTIBLE (0)
#define STATUS_EXPOSED     (1)
#define STATUS_INFECTIOUS  (2)
#define STATUS_RECOVERED   (3)

struct TxEvent {
    IND_TYPE    source;
    IND_TYPE    target;
    uint32_t    time_step;
};
