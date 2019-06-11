#ifndef statebuf_h
#define statebuf_h

#include "va_model.h"

typedef struct state_buf_t{
    state_t **states;
    int slots;
    int slot_size;
    int size;
    int curr_slot;
} state_buf_t;

state_buf_t *create_buffer(int slots, int slot_size);

void buf_read(state_buf_t *buf, state_t *res, int index);
void buf_write(state_buf_t *buf, state_t *state, int index, int rel_slot);
void buf_add(state_buf_t *buf, state_t *state, int index, int rel_slot);

void buf_read_all(state_buf_t *buf, state_t *res);
void buf_write_all(state_buf_t *buf, state_t *states, int rel_slot);
void buf_add_all(state_buf_t *buf, state_t *states, int rel_slot);
#endif