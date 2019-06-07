#ifndef statebuf_h
#define statebuf_h

#include "neuro.h"

typedef struct state_buf_t{
    state_t **states;
    int slots;
    int slot_size;
    int size;
    int curr_slot;
} state_buf_t;

state_buf_t *create_buffer(int slots, int slot_size);
void state_buf_read(state_buf_t *buf, state_t *res);
void state_buf_write(state_buf_t *buf, state_t *states, int state_len, int rel_slot);
void state_buf_add(state_buf_t *buf, state_t *state, int index, int rel_slot);

#endif