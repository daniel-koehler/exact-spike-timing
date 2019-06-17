#include <stdlib.h>
#include <stdio.h>
#include "spike.h"


int compare_spikes(const void *p1, const void *p2){
    /*
    Compare function for qsort, that determines the interspike interval between the given spikes
    */
    spike_t *spike1 = *(spike_t **) p1;
    spike_t *spike2 = *(spike_t **) p2;

    if(spike1->t < spike2->t) return -1;
    if(spike1->t > spike2->t) return  1;
    else return 0;
}

spike_t *sort_spikes(spike_t *head, int spike_cnt){
    /* 
    Sorts spikes in a linked list based on their time t returns the new head of the list.
    */
    spike_t *curr_spike = head;
    spike_t **spike_array = (spike_t **) malloc(sizeof(spike_t *)*spike_cnt);
    int i = 0;
    if(spike_cnt > 0){
        while(curr_spike){
            spike_array[i++] = curr_spike;
            curr_spike = curr_spike->next;
        }
        qsort(spike_array, spike_cnt, sizeof(spike_t *), compare_spikes);
        for(i = 0; i < spike_cnt - 1; i++){
            spike_array[i]->next = spike_array[i+1];
        }
        spike_array[spike_cnt - 1]->next = NULL;
        curr_spike = spike_array[0];
        free(spike_array);
        return curr_spike;
    }
    else return NULL;
}

void print_spike(spike_t *spike){
    printf("Neuron %i at %f ms\n", spike->index, spike->t);
}

void print_spikes(spike_t *head){
    /*
    Print all spikes contained in a linked list.
    */
    spike_t *curr_spike = head;
    while(curr_spike){
        print_spike(curr_spike);
        curr_spike = curr_spike->next;
    }
}

spike_t *append_spike(spike_t *spike, float t, int idx){
    /*
    Appends a spike to node 'spike' and returns pointer to it.
    */
    spike->next  = (spike_t *) malloc(sizeof(spike_t));
    spike        = spike->next;
    spike->next  = NULL;
    spike->index = idx;
    spike->t     = t;
    return spike;
}

void free_spikes(spike_t *top_input){
    spike_t *tmp;
    while(top_input){
        tmp = top_input;
        top_input = top_input->next;
        free(tmp);
    }
}