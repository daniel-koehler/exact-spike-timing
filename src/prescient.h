#ifndef prescient_h
#define prescient_h

#include <stdio.h>
#include "spike.h"
#include "va_model.h"
#include "statebuf.h"

#define frand() ((double) rand() / (RAND_MAX + 1.0))

typedef int bool;
#define true 1
#define false 0

typedef struct sim {
    /* simulation parameters */
    float h;                // simulation time step
    int   n;                // total number of neurons
    int   n_ex;
    int   n_in;
    int   n_syn;
    float r_ei;             // ratio of exc. to inh. neurons
    float p_conn;           // connection probability
    float min_delay;        // must fulfill h <= min_delay
    float max_delay;
    bool  rand_delays;
    bool  rand_states;      // initialize state variables randomly
    factor_sel_t  calc_factors;     // 'Lookup': use lookup table, 'Calculate': exact factors calculated for each time interval

    float t_start;
    float t_end;
    float t_input;
    float t_avg;            // average ISI for external stimulation

    /* Factors and for exponential integration */
    float *factors_h;
    float *factors_dt;
    
    interpolation_t interpolation;

    /* external stimulus / spikes */
    spike_t *next_input;
    spike_t *top_input;

    /* internal spikes (used to plot ISI distribution) */
    spike_t *next_spike;
    spike_t *top_spike;
    int spike_cnt;
    state_buf_t *state_buf;
    state_t *state_mem;             // stores state variables for each neuron
    state_t *tmp_mem;               // holds copy of state_mem
    
    state_t *buffered_state;
    state_t *update;

    synapse_t **synapses;           // stores synaptic connections between neurons
} sim_t;

void setup_sim(sim_t *sim);
void clear_sim(sim_t *sim);
void create_events(sim_t *sim);
void process_input(sim_t *sim, float t);
void subthreshold_dynamics(sim_t *sim, int i);
void generate_spike(sim_t *sim, int i, float t, FILE *fd_raster);

#endif