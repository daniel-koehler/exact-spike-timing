#include <stdlib.h>
#include <math.h>
#include "canonical.h"

void setup_sim(sim_t *sim){
    /*
    Implementation-specific setup of simulation structure.
    */
}

void clear_sim(sim_t *sim){
}

void create_events(sim_t *sim){
    /*
    Creates a poisson distributed input spike train for each neuron.
    */
    int n = sim->n;
    int i;
    float dt;
    int     spike_cnt    = 0;
    spike_t *new_input   = NULL;
    spike_t *top_input   = NULL;
    spike_t *next_input  = NULL;
    spike_t **spikes = (spike_t **) malloc(sizeof(spike_t*) * n);
    //sim->spikes = spikes 

    for(i = 0; i < n; i++){
        spikes[i] = NULL;
        next_input = NULL;
        for(float t = sim->t_start; t <= sim->t_input; ){
            /* Add up exponentially distributed intervals to generate Poisson spike train */
            dt  = -log(1.0 - frand()) * sim->t_avg;
            t  += dt;
            if (t > sim->t_input) break;

            /* Create a new input spike */
            new_input        = (spike_t *) malloc(sizeof(spike_t));
            new_input->t     = t;
            new_input->index = i;
            new_input->next  = NULL;

            if (spikes[i] == NULL){
                spikes[i] = next_input = new_input;
            }
            else{
                next_input->next = new_input;
            }
            next_input = new_input;
            spike_cnt++;
        }
    }
    sim->spikes = spikes;
    //sim->top_input  = sort_spikes(top_input, spike_cnt);
}

void process_input(sim_t *sim, float t){}
void subthreshold_dynamics(sim_t *sim, int i){}
void generate_spike(sim_t *sim, int i, float t, FILE *fd_raster){}