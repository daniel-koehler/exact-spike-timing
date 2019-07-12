#include <stdlib.h>
#include <stdio.h>
#include "spike.h"

spike_t *new_spike(int index, float t, float weight){
    /*
    Allocate and initialize new spike.
    */
    spike_t *spike = (spike_t *) malloc(sizeof(spike_t));
    spike->index  = index;
    spike->t      = t;
    spike->weight = weight;
    spike->next   = NULL;
    return spike;
}

int compare_spikes(const void *p1, const void *p2){
    /*
    Compare function for qsort, that determines the interspike interval between the given spikes
    */
    spike_t *spike1 = *(spike_t **) p1;
    spike_t *spike2 = *(spike_t **) p2;

    if (spike1->t < spike2->t) return -1;
    if (spike1->t > spike2->t) return  1;
    else return 0;
}

spike_t *sort_spikes(spike_t *head, int spike_cnt){
    /* 
    Sorts spikes in a linked list based on their time t returns the new head of the list.
    */
    spike_t *curr_spike = head;
    spike_t **spike_array = (spike_t **) malloc(sizeof(spike_t *)*spike_cnt);
    int i = 0;
    if (spike_cnt > 0){
        while (curr_spike){
            spike_array[i++] = curr_spike;
            curr_spike = curr_spike->next;
        }
        qsort(spike_array, spike_cnt, sizeof(spike_t *), compare_spikes);
        for (i = 0; i < spike_cnt - 1; i++){
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
    printf("Neuron %i at %f ms with weight %f.\n", spike->index, spike->t, spike->weight);
}

void print_spikes(spike_t *head){
    /*
    Print all spikes contained in a linked list.
    */
    spike_t *curr_spike = head;
    while (curr_spike){
        print_spike(curr_spike);
        curr_spike = curr_spike->next;
    }
}

void append_spike(spike_t **head, spike_t *spike){
    /*
    Appends 'spike' to a linked list, whose start is given by 'head'.
    */
    if (!head){  // Special case: empty list
        *head = spike;
        return;
    }
    spike_t *curr = *head;
    while (curr){
        if (curr->next){
            curr = curr->next;
        }
        else{
            curr->next = spike;
            break;
        }
    }

}

void sortin_spike(spike_t **head, spike_t *spike){
    /*
    Sorts 'spike' w.r.t its time t into a linked list, whose start is given by 'head'.
    */
    if (!(*head) || (*head)->t >= spike->t){
        spike->next = *head;
        *head = spike;
        return;
    }
    spike_t *curr = *head;
    while (curr->next && curr->next->t < spike->t){
        curr = curr->next;
    }
    spike->next = curr->next;
    curr->next = spike;
}

void free_spikes(spike_t *top_input){
    spike_t *tmp;
    while (top_input){
        tmp = top_input;
        top_input = top_input->next;
        free(tmp);
    }
}