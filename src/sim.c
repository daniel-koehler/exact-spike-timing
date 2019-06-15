/*
TODO:
- calculate c3 in va_model.c in dependence of I_inj
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "statebuf.h"
#include "va_model.h"
#include "spike.h"

#define true 1
#define false 0
#define frand() ((double) rand() / (RAND_MAX + 1.0))
#define N 1000

typedef enum interpolation_t {
    NONE,
    LINEAR,
    QUADRATIC,
    CUBIC
} interpolation_t;

typedef struct synapse_t{
    int target;
    int type;       // 0: inhibitory, 1: excitatory
    float weight;
    float delay;
} synapse_t;

typedef struct sim_t {
    /* simulation parameters */
    float h;                        // simulation time step
    int n;                          // total number of neurons
    int n_ex;
    int n_in;
    int n_syn;
    float ratio_ex_in;              // ratio of exc. to inh. neurons
    float p_conn;                   // connection probability
    float min_delay;
    float max_delay;
    int rand_delays;

    float t_start;
    float t_end;
    float t_input;

    float *factors_h;               // integration factors for time step h
    float *factors_dt;   
    
    interpolation_t interpolation;  // order of interpolation for exact spike timing

    /* external stimulus / spikes */
    spike_t *next_input;
    spike_t *top_input;

    /* internal spikes (used to plot ISI distribution) */
    spike_t *next_spike;
    spike_t *top_spike;
    int spike_cnt;

    state_t *state_mem;             // stores state variables for each neuron
    state_t *tmp_mem;               // holds copy of state_mem
    state_buf_t *state_buf;
    state_t *buffered_state;
    state_t *update;

    synapse_t **synapses;           // stores synaptic connections between neurons
  
    /* FILE objects used to save simulation results */
    FILE *fd_raster;
    FILE *fd_voltage;
    FILE *fd_conductance;
    FILE *fd_isi;
    FILE *fd_stats;
} sim_t;

void create_events(sim_t *sim, float t_start, float t_end, float t_avg){
    /*
    Creates a poisson distributed input spike train for each neuron.
    */
    int n = sim->n;
    int i;
    float dt;
    int     spike_cnt    = 0;
    spike_t *new_input   = NULL;
    spike_t *top_input = NULL;
    spike_t *next_input  = NULL;

    for(i = 0; i < n; i++){
        for(float t = t_start; t <= t_end; ){
            dt  = -log(1.0 - frand()) * t_avg;
            t  += dt;

            if (t > t_end){
                break;
            }

            /* Create a new input spike */
            new_input        = (spike_t *) malloc(sizeof(spike_t));
            new_input->t     = t;
            new_input->index = i;
            new_input->next  = NULL;

            if (top_input == NULL){
                sim->top_input = sim->next_input = top_input = new_input;
            }
            else{
                next_input->next = new_input;
            }
            next_input = new_input;
            spike_cnt++;
        }
    }
    //sim->top_input = sim->next_input = sort_spikes(top_input, spike_cnt);
    sim->top_input  = sort_spikes(top_input, spike_cnt);
}

void initialize_state_mem(sim_t *sim){
    /*
    Initialize state variables of all neurons
    */
    int i;
    state_t *state_mem;
	float h          = sim->h;
	int   n          = sim->n;
	float max_delay  = sim->max_delay;
	float min_delay  = sim->min_delay;
	int   slots      = (max_delay-min_delay)/h + 1;

	int   slot_size  = n;
    sim->state_mem = state_mem = (state_t *) malloc(sizeof(state_t)*n);
    sim->state_buf = create_buffer(slots, slot_size);
    sim->tmp_mem        = (state_t *) malloc(sizeof(state_t) * n);
    sim->buffered_state = (state_t *) malloc(sizeof(state_t) * n);
    sim->update         = (state_t *) malloc(sizeof(state_t));


    for(i = 0; i < n; i++){
        state_mem[i] = (state_t) {
            .t_ela  = tau_ref,
            //.V_m    = frand() * (V_th - E_rest) + E_rest,
            .V_m    = 0.0,
            .g_ex   = 0.0,
            .g_in   = 0.0};
    }
}

void create_network(sim_t *sim){
    /*
    Create synapses that represent the neural network
    */
    int i, j;
    int   n_syn      = sim->n_syn;
    int   n          = sim->n;
    int   n_ex       = sim->n_ex;
    float max_delay  = sim->max_delay;
    float min_delay  = sim->min_delay;
    float delay, weight;
    synapse_t **synapses = (synapse_t **) malloc(sizeof(synapse_t *) * n);
    for(i = 0; i < n; i++){
        synapses[i] = (synapse_t *) malloc(sizeof(synapse_t) * n_syn);
        for(j = 0; j < n_syn; j++){
            /* Target index */
            synapses[i][j].target = floor(frand()*n);

            /* Synaptic propagation delays */
            if(sim->rand_delays){               
                // TODO: let delays be multiple of h
                delay = frand()*(max_delay - min_delay) + min_delay;
            }
            else     delay      = min_delay;
            synapses[i][j].delay = delay;

            /* Synaptic weights */
            if(i < n_ex){
                weight = dg_ex;
                synapses[i][j].type = 1;              
            }
            else{
                weight = dg_in;
                synapses[i][j].type = 0;
            }
            synapses[i][j].weight = weight;            
        }
    }
    sim->synapses = synapses;
}

sim_t *setup_sim(void){
    /*
    Set simulation parameters, initalize memory, generate stimulus
    */

    sim_t *sim = (sim_t *) malloc(sizeof(sim_t));
    float h;
    int   n;
    int   n_ex;
    float ratio;
    float p_conn;
    float t_start;
    float t_end;
    float t_input;
    float t_avg   = 4.0;

    /* simulation parameters */
    sim->h       = h       = 0.1;
    sim->t_start = t_start = 0.0;
    sim->t_end   = t_end   = 100.0;
    sim->t_input = t_input = 50.0;

    sim->interpolation = NONE; 

    sim->n           = n      = N;
    sim->ratio_ex_in = ratio  = 4;
    sim->p_conn      = p_conn = 0.02;
    sim->min_delay   = sim->h;
    sim->n_ex        = n_ex   = floor(n * ratio / (ratio + 1));
    sim->n_in        = n - n_ex;
    sim->n_syn       = (n - 1) * p_conn;
    
    sim->min_delay   = h;
    sim->max_delay   = 10*h;
    sim->rand_delays = false;
    
    sim->top_spike   = NULL;
    sim->next_spike  = NULL;
    sim->spike_cnt   = 0;

    create_events(sim, t_start, t_input, t_avg);

    create_network(sim);

    initialize_state_mem(sim);

    sim->factors_h  = (float *) malloc(sizeof(float) * 7);
    sim->factors_dt = (float *) malloc(sizeof(float) * 7);
    calc_factors(h, sim->factors_h);
    
    /* Open output files */
    sim->fd_raster      = fopen("results/raster", "w+");
    sim->fd_voltage     = fopen("results/voltage", "w+");
    sim->fd_conductance = fopen("results/conductance", "w+");
    sim->fd_isi         = fopen("results/isi", "w+");
    sim->fd_stats       = fopen("results/stats.txt", "w+");
    if ( (!sim->fd_raster) || (!sim->fd_voltage) || (!sim->fd_isi) || (!sim->fd_stats) || (!sim->fd_conductance)){
        printf("Could not open files to save results. Exit.\n");
        exit(1);
    }
    return(sim);

}

void clear_sim(sim_t *sim){
    /*
    Free allocated memory
    */
    if (sim){      
        if (sim->synapses){
            for(int i = 0; i < sim->n; i++){
                free(sim->synapses[i]);
            }
            free(sim->synapses);
        }
        if (sim->state_mem)         free(sim->state_mem);
        if (sim->factors_h)         free(sim->factors_h);
        if (sim->factors_dt)        free(sim->factors_dt);
        if (sim->state_buf)         free_buffer(sim->state_buf);
        if (sim->tmp_mem)           free(sim->tmp_mem);
        if (sim->buffered_state)    free(sim->buffered_state);
        if (sim->update)            free(sim->update);
        if (sim->top_input)         free_spikes(sim->top_input);
        if (sim->top_spike)         free_spikes(sim->top_spike);

        if (sim->fd_raster)         fclose(sim->fd_raster);
        if (sim->fd_voltage)        fclose(sim->fd_voltage);
        if (sim->fd_conductance)    fclose(sim->fd_conductance);
        if (sim->fd_isi)            fclose(sim->fd_isi);
        if (sim->fd_stats)          fclose(sim->fd_stats);
    }
    free(sim);
    
}

void print_state_mem(sim_t *sim){
    state_t *state_mem = sim->state_mem;
    int   n          = sim->n;
    for(int i=0; i<n; i++){
        printf("Index %i: t_ela = %.2f, V_m = %.2f, g_ex = %.2f, g_in = %.2f\n", i, state_mem[i].t_ela,\
        state_mem[i].V_m, state_mem[i].g_ex, state_mem[i].g_in);
    }
}

void calc_update(state_t *update, float *factors, float g_ex, float g_in){
    update->V_m = 0;
    update->g_ex = g_ex;
    update->g_in = g_in;
    solve_analytic(update, factors);
    update->t_ela = 0;
}

float calc_isi(spike_t *head, int n, int spike_cnt, FILE *fd){
    /* 
    Calculates interspike intervals (ISIs) and writes them into a file. Returns avergage ISI.
    */
    head = sort_spikes(head, spike_cnt);
    if(head){
        spike_t *curr_spike = head->next;
        float t0[n];
        int i, j = 0;
        float isi;
        float sum = 0.0;
        for(i = 0; i < n; i++){
            t0[i] = 0;
        }
        while(curr_spike){
            i = curr_spike->index;
            if (t0[i] == 0){    // first spike of neuron i
                t0[i] = curr_spike->t;   
            }
            else{
                isi  = curr_spike->t - t0[i];
                sum += isi;
                j++;
                fprintf(fd, "%f\n", curr_spike->t - t0[i]);
                t0[i] = curr_spike->t;         
            }
            curr_spike = curr_spike->next;        
        }
        return sum/j;
    }
    else return 0.0;
}

void statistics(sim_t *sim){
    /*
    Calculate and save various statistics for the current simulation run.
    */
    state_t *state_mem = sim->state_mem;

    float mean_isi = calc_isi(sim->top_spike, sim->n, sim->spike_cnt, sim->fd_isi);
    float mean_V_m = 0.0, mean_g_ex = 0.0, mean_g_in = 0.0;
    int spike_cnt  = sim->spike_cnt;
    for (int i = 0; i < sim->n; i++){
        mean_V_m  += state_mem[i].V_m;
        mean_g_ex += state_mem[i].g_ex;
        mean_g_in += state_mem[i].g_in;
    }
    mean_V_m  /= sim->n;
    mean_g_ex /= sim->n;
    mean_g_in /= sim->n;

    fprintf(sim->fd_stats, "Number of spikes: %i\n", spike_cnt);
    fprintf(sim->fd_stats, "Average membrane voltage: %.2f mV\n", mean_V_m);
    fprintf(sim->fd_stats, "Average exc. conductance: %.2f nS\n", mean_g_ex);
    fprintf(sim->fd_stats, "Average inh. conductance: %.2f nS\n", mean_g_in);
    fprintf(sim->fd_stats, "Average interspike interval: %.2f ms\n", mean_isi);
}

void process_input(sim_t *sim, float t){
    /*
    Calculate influence of input spikes and write it to corresponding entry in state_buf.
    */
    spike_t *spike = sim->next_input;
    float t_s;
    while(spike){
        if (spike->t > t + sim->min_delay) break;
        t_s = spike->t - t;
        calc_factors(sim->h - t_s, sim->factors_dt);
        calc_update(sim->update, sim->factors_dt, dg_stim, 0);
        buf_add(sim->state_buf, sim->update, spike->index, 0);
        spike = spike->next;
    }
    sim->next_input = spike;
}

void subthreshold_dynamics(sim_t *sim, int i){
    /*
    Calculate subthreshold dynamics for neuron i.
    */
    float t_em;
    float h = sim->h;
    float *factors_dt = sim->factors_dt;
    state_t     *state_mem      = sim->state_mem;
    state_t     *buffered_state = sim->buffered_state;

    /* Standard timing */
    if (sim->interpolation == NONE){   
        solve_analytic(&state_mem[i], sim->factors_h);
        add_state(&state_mem[i], &buffered_state[i]);
    }

    /* Exact timing */
    else{ 
        memcpy(sim->tmp_mem, state_mem, sizeof(state_t) * sim->n);    // copy of state variables at time t
        /* Neuron emerges from refractory period this update interval */
        if (state_mem[i].t_ela > tau_ref - h && state_mem[i].t_ela <= tau_ref){
            t_em = state_mem[i].t_ela + h - tau_ref;
            calc_factors(t_em, factors_dt);
            solve_analytic(&state_mem[i], factors_dt);
            state_mem[i].V_m = 0;
            calc_factors(h - t_em, factors_dt);
            solve_analytic(&state_mem[i], factors_dt);
            add_state(&state_mem[i], &buffered_state[i]);
            state_mem[i].V_m -= t_em / h * buffered_state[i].V_m;
        }
        /* Neuron does not emerge from refractory period */
        else{
            solve_analytic(&state_mem[i], sim->factors_h);
            add_state(&state_mem[i], &buffered_state[i]);
        }
    }
}

void process_spike(sim_t *sim, int i, float t){
    state_t     *state_mem      = sim->state_mem;
    state_t     *tmp_mem        = sim->tmp_mem;
    state_t     *update         = sim->update;
    state_buf_t *state_buf      = sim->state_buf;
    float       *factors_dt     = sim->factors_dt;
    synapse_t   **synapses      = sim->synapses;
    float h = sim->h;
    float t_s, delay;
    float y0_dot, yh_dot;
    int target;

    /* Calculate spike time */
    switch (sim->interpolation){
        case NONE:
            t_s = h;
            break;
        
        case LINEAR:
            t_s = linear_int(tmp_mem[i].V_m, state_mem[i].V_m, V_th, h);
            break;
        
        case QUADRATIC:
            y0_dot = voltage_deriv(0.0, tmp_mem[i].V_m, tmp_mem[i].g_ex, tmp_mem[i].g_in);
            t_s = quadratic_int(tmp_mem[i].V_m, y0_dot, state_mem[i].V_m, V_th, h);
            break;

        case CUBIC: // not implemented yet
            y0_dot = voltage_deriv(0.0, tmp_mem[i].V_m, tmp_mem[i].g_ex, tmp_mem[i].g_in);
            yh_dot = voltage_deriv(h, tmp_mem[i].V_m, tmp_mem[i].g_ex, tmp_mem[i].g_in);
            t_s = cubic_int(tmp_mem[i].V_m, y0_dot, state_mem[i].V_m, yh_dot, V_th, h);
            break;
    }
    fprintf(sim->fd_raster, "%f %i\n", t + t_s, i);

    /* Calculate update */
    // Note: if a single neuron can have both excitatory and inhibitory synapses
    // the calculation of the update has to be done in the loop below  

    calc_factors(h - t_s, factors_dt);            
    if (synapses[i][0].type == 0){  // inhibitory synapse
        calc_update(update, factors_dt, 0, synapses[i][0].weight);
    }
    else{                           // excitatory synapse
        calc_update(update, factors_dt, synapses[i][0].weight, 0);
    }
    
    /* Propagate spike to each target neuron */
    for(int j = 0; j < sim->n_syn; j++){
        target = synapses[i][j].target;
        delay  = synapses[i][j].delay / h;
        buf_add(state_buf, update, target, delay);
    }

    state_mem[i].t_ela = 0;

    /* Save spikes in linked list */
    if (!sim->top_spike){
        sim->top_spike = sim->next_spike = (spike_t *) malloc(sizeof(spike_t));
    }
    else{
        sim->next_spike->next = (spike_t *) malloc(sizeof(spike_t));
        sim->next_spike = sim->next_spike->next;
    }
    sim->next_spike->next   = NULL;
    sim->next_spike->index  = i;
    sim->next_spike->t      = t + t_s;
    sim->spike_cnt         += 1;
}

void simulation_loop(sim_t *sim){
    /*
    Main loop of the time based simulation
    */

    float t_start = sim->t_start;
    float t_end   = sim->t_end;
    float t;
    float t_output = 0.1 * (t_end - t_start);
    float h = sim->h;
    int   n = sim->n;

    state_t     *state_mem = sim->state_mem;

    for(t=t_start; t <= t_end; t+=h){
        if (t >= t_output){
            printf("t = %.0f ms\n", t);
            t_output += 0.1 * (t_end - t_start);
        }

        /* Update interval (t, t+h] */
        process_input(sim, t);
        buf_read_all(sim->state_buf, sim->buffered_state);
        for(int i = 0; i < n; i++){
            subthreshold_dynamics(sim, i);
            
            /* Neuron spikes */
            if (state_mem[i].V_m >= V_th){
                process_spike(sim, i, t);
            }

            /* Neuron is in refractory period - clamp to resting potential */
            if (state_mem[i].t_ela < tau_ref){
                state_mem[i].V_m = E_rest;
            }

            /* Write state variables to output files */
            fprintf(sim->fd_voltage, "%f %f\n", t, state_mem[0].V_m);
            fprintf(sim->fd_conductance, "%f %f %f\n", t, state_mem[0].g_ex, -state_mem[0].g_in);
        }
    }
}

int main(void){
    sim_t *sim = setup_sim();
    simulation_loop(sim);
    statistics(sim);
    clear_sim(sim);   
}
