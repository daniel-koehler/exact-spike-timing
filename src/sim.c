#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "statebuf.h"
#include "va_model.h"

#define true 1
#define false 0
#define frand() ((double) rand() / (RAND_MAX + 1.0))

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

typedef struct spike_t {
    float   t;
    int     index;
    struct spike_t *next;
} spike_t;

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
    
    interpolation_t interpolation;  // order of interpolation for exact spike timing

    spike_t *next_spike;
    spike_t *top_spike;
    int spike_cnt;

    state_t *state_mem;             // stores state variables for each neuron
    state_buf_t *state_buf;
    synapse_t **synapses;            // stores synaptic connections between neurons
} sim_t;

int compare_spikes(const void *p1, const void *p2){
    /*
    Compare function for qsort, that determines the interspike interval between the given spikes
    */
    spike_t *spike1 = *(spike_t **) p1;
    spike_t *spike2 = *(spike_t **) p2;

    if(spike1->t < spike2->t) return -1;
    if(spike1->t > spike2->t) return 1;
    else return 0;
}

void create_events(sim_t *sim, float t_start, float t_end, float t_avg){
    /*
    Creates a poisson distributed input spike train for each neuron.
    */
    int n = sim->n;
    int i;
    float dt;
    int     spike_cnt    = 0;
    spike_t *new_spike   = NULL;
    spike_t *first_spike = NULL;
    spike_t *next_spike  = NULL;

    for(i = 0; i < n; i++){
        for(float t = t_start; t <= t_end; ){
            dt  = -log(1.0 - frand()) * t_avg;
            t  += dt;

            if (t > t_end){
                break;
            }

            /* Create a new input spike */
            new_spike        = (spike_t *) malloc(sizeof(spike_t));
            new_spike->t     = t;
            new_spike->index = i;
            new_spike->next  = NULL;

            if (first_spike == NULL){
                sim->top_spike = sim->next_spike = first_spike = new_spike;
            }
            else{
                next_spike->next = new_spike;
            }
            next_spike = new_spike;
            spike_cnt++;
        }
    }
    sim->spike_cnt = spike_cnt;

    /* Sort spikes */
    i             = 0;
    next_spike    = first_spike;
    spike_t **spike_array = (spike_t **) malloc(sizeof(spike_t *)*spike_cnt);
    while(next_spike){
        spike_array[i++] = next_spike;
        next_spike = next_spike->next;
    }
    qsort(spike_array, spike_cnt, sizeof(spike_t *), compare_spikes);
    for(i = 0; i < spike_cnt - 1; i++){
        spike_array[i]->next = spike_array[i+1];
    }
    spike_array[spike_cnt - 1]->next = NULL;
    sim->top_spike = spike_array[0];
    free(spike_array);
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
    
    for(i = 0; i < n; i++){
        state_mem[i] = (state_t) {
            .t_ela  = tau_ref,
            //.V_m    = frand() * (V_th - E_rest) + E_rest,
            .V_m = 0.0,
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

    sim->n           = n      = 10;
    sim->ratio_ex_in = ratio  = 4;
    sim->p_conn      = p_conn = 0.02;
    sim->min_delay   = sim->h;
    sim->n_ex  = n_ex = floor(n * ratio / (ratio + 1));
    sim->n_in  = n - n_ex;
    sim->n_syn = (n - 1) * p_conn;
    
    sim->min_delay   = h;
    sim->max_delay   = 2*h;
    sim->rand_delays = false;

    /* create stimulus */
    sim->top_spike = (spike_t *) malloc(sizeof(spike_t));
    sim->next_spike  = (spike_t *) malloc(sizeof(spike_t));
    create_events(sim, t_start, t_input, t_avg);

    create_network(sim);

    initialize_state_mem(sim);

    sim->factors_h = (float *) malloc(sizeof(float) * 7);
    calc_factors(h, sim->factors_h);
    
    return(sim);
}

void print_state_mem(sim_t *sim){
    state_t *state_mem = sim->state_mem;
    int   n          = sim->n;
    for(int i=0; i<n; i++){
        printf("Index %i: t_ela = %.2f, V_m = %.2f, g_ex = %.2f, g_in = %.2f\n", i, state_mem[i].t_ela,\
        state_mem[i].V_m, state_mem[i].g_ex, state_mem[i].g_in);
    }
}

void print_spikes(sim_t *sim){
    spike_t *top_spike = sim->top_spike;
    while(top_spike){
        printf("Neuron %i at %f ms\n", top_spike->index, top_spike->t);
        top_spike = top_spike->next;
    }
}

void calc_update(state_t *update, float *factors, float g_ex, float g_in){
    update->V_m = 0;
    update->g_ex = g_ex;
    update->g_in = g_in;
    solve_analytic(update, factors);
    update->t_ela = 0;
}

void simulation_loop(sim_t *sim){
    /*
    Main loop of the time based simulation
    */

    float t_start = sim->t_start;
    float t_end = sim->t_end;
    float t;
    float h = sim->h;
    int n = sim->n;
    int n_syn = sim->n_syn;

    int i,j;
    float t_em, t_s;
    int target;
    float delay;
    float *factors_t1 = (float *) malloc(sizeof(float) * 7);
    float *factors_t2 = (float *) malloc(sizeof(float) * 7);

    spike_t *spike = sim->top_spike;
    spike_t *tmp;

    state_buf_t *state_buf      = sim->state_buf;
    synapse_t   **synapses      = sim->synapses;
    state_t     *state_mem      = sim->state_mem;
    state_t     *tmp_mem        = (state_t *) malloc(sizeof(state_t) * n);  // holds copy of state_mem
    state_t     *buffered_state = (state_t *) malloc(sizeof(state_t) * n);
    state_t     *update         = (state_t *) malloc(sizeof(state_t));

    for(t=t_start; t <= t_end; t+=h){
        /*if (fmod(t, 5) < h){
            printf("%f\n", t);
        }*/

        /* Update interval (t, t+h] */
        memcpy(tmp_mem, state_mem, sizeof(state_t) * n);
        buf_read_all(state_buf, buffered_state);
        
        /* Process input spikes */
        while (spike){
            if (spike->t > t + sim->min_delay){
                sim->top_spike = spike;
                break;
            }
            target = spike->index;
            t_s = spike->t - t;
            calc_factors(h - t_s, factors_t1);
            calc_update(update, factors_t1, dg_stim, 0);
            buf_add(sim->state_buf, update, target, 0);
            printf("Process spike for neuron %i at %.2f ms\n", target, spike->t);
            tmp = spike;
            spike = spike->next;
            free(tmp);
        }

        for(i = 0; i < n; i++){
            /************************* Subthreshold dynamics ***********************/
            /* Neuron emerges from refractory period this update interval */
            if (state_mem[i].t_ela > tau_ref - h || state_mem[i].t_ela <= tau_ref){
                t_em = state_mem[i].t_ela + h - tau_ref;
                calc_factors(t_em,     factors_t1);
                calc_factors(h - t_em, factors_t2);
                solve_analytic(&state_mem[i], factors_t1);
                state_mem[i].V_m = 0;
                
                solve_analytic(&state_mem[i], factors_t2);
                state_add(&state_mem[i], &buffered_state[i], &state_mem[i]);
                state_mem[i].V_m -= t_em / h * buffered_state[i].V_m;
            }
            /* Neuron does not emerge from refractory period */
            else{
                solve_analytic(&state_mem[i], sim->factors_h);
                state_add(&state_mem[i], &buffered_state[i], &state_mem[i]);
            }
            /************************* Collect spikes ***********************/
            /* Neuron spikes */
            if (state_mem[i].V_m >= V_th){
                t_s = linear_int(tmp_mem[i].V_m, state_mem[i].V_m, V_th, h);
                /* Calculate update */
                calc_factors(h - t_s, factors_t1);
                // Note: if a single neuron can have both excitatory and inhibitory synapses
                // the calculation of the update has to be done in the loop below
                if (synapses[i][j].type == 0){
                    calc_update(update, factors_t1, 0, synapses[i][j].weight);
                }
                else{
                    calc_update(update, factors_t1, synapses[i][j].weight, 0);
                }
                for(j = 0; j < n_syn; j++){
                    target = synapses[i][j].target;
                    delay  = synapses[i][j].delay / h;  // TODO check if the result is right
                    buf_add(state_buf, update, target, delay);
                }

                state_mem[i].t_ela = 0;
            }

            /* Neuron is in refractory period - clamp to resting potential */
            if (state_mem[i].t_ela < tau_ref){
                state_mem[i].V_m = E_rest;
            }
        }
    }
}

int main(void){
    sim_t *sim = setup_sim();
    state_t *buffered_state = (state_t *) malloc(sizeof(state_t) * sim->n);
    buf_read_all(sim->state_buf, buffered_state);
    state_t *states = (state_t *) malloc(sizeof(state_t) * sim->n);
    buf_write(sim->state_buf, states, sim->n, 0);
    state_t *state = (state_t *) malloc(sizeof(state_t));
    buf_add(sim->state_buf, state, 3, 3);
    simulation_loop(sim);
    //printf("%f, %f\n", sim->state_mem[0], sim->model.tau_ref);
    print_state_mem(sim);
    //print_spikes(sim);
}

