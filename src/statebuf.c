#include <stdlib.h>
#include "statebuf.h"
#include <stdio.h>

state_buf_t *create_buffer(int slots, int slot_size){
    /*
    Creates a state buffer. 'states' is a multidimensional array of size slots x slot_size.
    */
    state_buf_t *buffer = (state_buf_t *) malloc(sizeof(state_buf_t));
    buffer->slots     = slots;
    buffer->slot_size = slot_size;
    buffer->size      = slots * slot_size;
    buffer->states = (state_t **) malloc(sizeof(state_t*) * slots);
    for(int i = 0; i < slots; i++){
        buffer->states[i] = malloc(sizeof(state_t) * slot_size);
        for(int j = 0; j < slot_size; j++){
            buffer->states[i][j] = (state_t) {
                .t_ela  = 0.0,
                .V_m    = 0.0,
                .g_ex   = 0.0,
                .g_in   = 0.0};
        }
    }
    buffer->curr_slot = 0;
    return buffer;
}

void state_buf_read(state_buf_t *buf, state_t *res){
    /*
    Reads the state variables buffered for current time slot and writes them to the address pointed by *res.
    */
    for(int i = 0; i < buf->slot_size; i++){
        res[i] = buf->states[buf->curr_slot][i];
    }
    buf->curr_slot += 1;
    if(buf->curr_slot >= buf->slots) buf->curr_slot = 0;
}

void state_buf_write(state_buf_t *buf, state_t *states, int state_len, int rel_slot){
    /*
    */
    if (buf->slot_size != state_len){
        printf("Passed argument 'states' has wrong size.\n");
        return;
    }
    else if (rel_slot > buf->slots || rel_slot < 0){
        printf("Invalid 'rel_slot': %i.\n", rel_slot);
        return;
    }
    int abs_slot = buf->curr_slot + rel_slot;
    if(abs_slot >= buf->slots) abs_slot -= buf->slots;
    for(int i = 0; i < buf->slot_size; i++){
        buf->states[abs_slot][i] = states[i];
    }
}

void state_buf_add(state_buf_t *buf, state_t *state, int index, int rel_slot){
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