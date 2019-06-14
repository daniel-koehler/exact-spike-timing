/*
Defines 'state_buf_t' which is a buffer that can hold the state variables of multiple 
neurons for multiple time slots. E.g. create_buffer(3, 5) returns a buffer that holds
the state variables of 5 neurons for 3 different time slots.

Also provides functions to read, write and add state variables to the buffer.
*/

#include <stdlib.h>
#include <stdio.h>
#include "statebuf.h"


state_buf_t *create_buffer(int slots, int slot_size){
    /*
    Creates a buffer that stores state variables in 'states' (a multidimensional array of size slots x slot_size).
    */
    state_buf_t *buf = (state_buf_t *) malloc(sizeof(state_buf_t));
    buf->slots     = slots;
    buf->slot_size = slot_size;
    buf->states = (state_t **) malloc(sizeof(state_t*) * slots);
    for(int i = 0; i < slots; i++){
        buf->states[i] = malloc(sizeof(state_t) * slot_size);
        for(int j = 0; j < slot_size; j++){
            buf->states[i][j] = ZERO_STATE;
        }
    }
    buf->curr_slot = 0;
    return buf;
}

void free_buffer(state_buf_t * buf){
    /*
    Free allocated memory of buf.
    */
    for(int i = 0; i < buf->slots; i++){
        free(buf->states[i]);
    }
    free(buf->states);
    free(buf);
}

void buf_read(state_buf_t *buf, state_t *res, int index){
    /*
    Reads the state variables of neuron 'index' buffered for current time slot and writes them to the address pointed by *res.
    */
    *res = buf->states[buf->curr_slot][index];
}

void buf_write(state_buf_t *buf, state_t *state, int index, int rel_slot){
    /*
    Writes state variables in 'state' into the buffer for neuron 'index'.
    'rel_slot' is the slot index relative to 'curr_slot' of the buffer.
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
    Increments state variables for neuron 'index' by the ones in 'state'.
    'rel_slot' is the slot index relative to 'curr_slot' of the buffer.
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
    state_add(&buf->states[abs_slot][index], state);
}

void buf_read_all(state_buf_t *buf, state_t *res){
    /*
    Reads buffered state variables for all neurons, clears the current slot of the buffer and increments 'curr_slot'.
    */   
    for(int i = 0; i < buf->slot_size; i++){
        buf_read(buf, &res[i], i);
        buf_write(buf, &ZERO_STATE, i, 0);
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