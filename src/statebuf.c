#include <stdlib.h>
#include <stdio.h>
#include "statebuf.h"


state_buf_t *create_buffer(int slots, int slot_size){
    /*
    Creates a state buffer. 'states' is a multidimensional array of size slots x slot_size.
    */
    state_buf_t *buf = (state_buf_t *) malloc(sizeof(state_buf_t));
    buf->slots     = slots;
    buf->slot_size = slot_size;
    buf->size      = slots * slot_size;
    buf->states = (state_t **) malloc(sizeof(state_t*) * slots);
    for(int i = 0; i < slots; i++){
        buf->states[i] = malloc(sizeof(state_t) * slot_size);
        for(int j = 0; j < slot_size; j++){
            buf->states[i][j] = (state_t) {
                .t_ela  = 0.0,
                .V_m    = 0.0,
                .g_ex   = 0.0,
                .g_in   = 0.0};
        }
    }
    buf->curr_slot = 0;
    return buf;
}

void buf_read(state_buf_t *buf, state_t *res, int index){
    /*
    Reads the state variables of neuron 'index' buffered for current time slot and writes them to the address pointed by *res.
    */
    res = &buf->states[buf->curr_slot][index];
}

void buf_write(state_buf_t *buf, state_t *state, int index, int rel_slot){
    /*
    */
    if (rel_slot > buf->slots || rel_slot < 0){
        printf("Invalid 'rel_slot': %i.\n", rel_slot);
        return;
    }
    int abs_slot = buf->curr_slot + rel_slot;
    if(abs_slot >= buf->slots) abs_slot -= buf->slots;
    buf->states[abs_slot][index] = *state;
}

void buf_add(state_buf_t *buf, state_t *state, int index, int rel_slot){
    /*
    */
   if (index >= buf->slot_size || index < 0){
        printf("Invalid 'index': %i.\n", index);
        return;
    }
    else if (rel_slot > buf->slots || rel_slot < 0){
        printf("Invalid 'rel_slot': %i.\n", rel_slot);
        return;
    }
    int abs_slot = buf->curr_slot + rel_slot;
    if(abs_slot >= buf->slots) abs_slot -= buf->slots;
    buf->states[abs_slot][index] = *state;      // TODO: is that right?
}

void buf_read_all(state_buf_t *buf, state_t *res){
    for(int i = 0; i < buf->slot_size; i++){
        buf_read(buf, &res[i], i);
    }
    buf->curr_slot += 1;
    if(buf->curr_slot >= buf->slots) buf->curr_slot = 0;
};

void buf_write_all(state_buf_t *buf, state_t *states, int rel_slot){  
    for(int i = 0; i < buf->slot_size; i++){
        buf_write(buf, &states[i], i, rel_slot);
    }
}
void buf_add_all(state_buf_t *buf, state_t *states, int rel_slot){
    for(int i = 0; i < buf->slot_size; i++){
        buf_add(buf, &states[i], i, rel_slot);
    }
}