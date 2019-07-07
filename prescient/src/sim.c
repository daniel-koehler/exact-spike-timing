#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "statebuf.h"
#include "va_model.h"
#include "spike.h"

clock_t M1, M2;
float P1;

#define frand() ((double) rand() / (RAND_MAX + 1.0))

typedef int bool;
#define true 1
#define false 0

typedef enum {Excitatory, Inhibitory} synapse_type;

typedef struct synapse{
    int target;
    synapse_type type;
    float weight;
    float delay;
} synapse_t;

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

/* FILE objects used to save simulation results */
FILE *fd_raster;
FILE *fd_voltage;
FILE *fd_conductance;

void print_state_mem(sim_t *sim){
    state_t *state_mem = sim->state_mem;
    int   n          = sim->n;
    for(int i=0; i<n; i++){
        printf("Index %i: t_ela = %.2f, V_m = %.2f, g_ex = %.2f, g_in = %.2f\n", i, state_mem[i].t_ela,\
        state_mem[i].V_m, state_mem[i].g_ex, state_mem[i].g_in);
    }
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

    for(i = 0; i < n; i++){
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
    sim->top_input  = sort_spikes(top_input, spike_cnt);
}

void setup_parameters(sim_t *sim){
    /* 
    Setup all important parameters of the simulation.
    */
    sim->n       = 4000;
    sim->h       = 0.1;
    sim->t_start = 0.0;
    sim->t_end   = 1000.0;
    sim->t_input = 50.0;
    sim->t_avg   = 4.0;
    sim->r_ei    = 4;
    sim->p_conn  = 0.02;
    sim->interpolation = Linear;
    sim->min_delay     = sim->h;
    sim->max_delay     = 10 * sim->h;
    sim->rand_delays   = false;
    sim->rand_states   = true;
    sim->calc_factors  = false;
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
        if (sim->rand_states){
            state_mem[i] = (state_t) {
            .t_ela = tau_ref + h,
            .V_m    = frand() * (V_th - E_rest) + E_rest,
            .g_ex   = frand() * dg_ex * 2,
            .g_in   = frand() * dg_in * 2};
        }
        else
        {
            state_mem[i] = (state_t) {
            .t_ela = tau_ref + h,
            .V_m    = 0.0,
            .g_ex   = 0.0,
            .g_in   = 0.0};
        }        
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
                delay = frand()*(max_delay - min_delay) + min_delay;
                delay = ((int)(delay / sim->h) * sim->h);   // let delay be multiple of h
                printf("Delay: %f\n", delay);
            }
            else     delay      = min_delay;
            synapses[i][j].delay = delay;

            /* Synaptic weights */
            if(i < n_ex){
                weight = dg_ex;
                synapses[i][j].type = Excitatory;              
            }
            else{
                weight = dg_in;
                synapses[i][j].type = Inhibitory;
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

    setup_parameters(sim);

    sim->n_ex      = floor(sim->n * sim->r_ei / (sim->r_ei + 1));
    sim->n_in      = sim->n - sim->n_ex;
    sim->n_syn     = (sim->n - 1) * sim->p_conn;
    
    sim->top_spike   = NULL;
    sim->next_spike  = NULL;
    sim->spike_cnt   = 0;

    create_events(sim);

    create_network(sim);

    initialize_state_mem(sim);

    generate_lut(sim->h, 100);      // generates lookup table with 100+1 entries for time interval [0 h]

    sim->factors_h  = (float *) malloc(sizeof(float) * NUM_FACTORS);
    sim->factors_dt = (float *) malloc(sizeof(float) * NUM_FACTORS);

    calc_factors(sim->h, sim->factors_h);
    
    /* Open output files */
    fd_raster      = fopen("results/raster", "w+");
    fd_voltage     = fopen("results/voltage", "w+");
    fd_conductance = fopen("results/conductance", "w+");
    if ( !fd_raster || !fd_voltage || !fd_conductance){
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
        //free_lut();
        free_buffer(sim->state_buf);
        if (sim->tmp_mem)           free(sim->tmp_mem);
        if (sim->buffered_state)    free(sim->buffered_state);
        if (sim->update)            free(sim->update);
        if (sim->top_input)         free_spikes(sim->top_input);
        if (sim->top_spike)         free_spikes(sim->top_spike);
        if (fd_raster)              fclose(fd_raster);
        if (fd_voltage)             fclose(fd_voltage);
        if (fd_conductance)         fclose(fd_conductance);
    }
    free(sim);
    
}

float calc_avgvoltage(state_t *state_mem, int n){
    FILE *fd = fopen("results/avgvoltage", "w+");
    int sum = 0;
    for (int i = 0; i < n; i++){
        fprintf(fd, "%f\n", state_mem[i].V_m);
        sum += state_mem[i].V_m;
    }
    fclose(fd);
    return sum/n;
}

float calc_firingrate(spike_t *head, int n, float t_span){
    /*
    Calculates the average firing rate of each neuron and returns average firing rate of all neurons.
    */
    if(head){
        FILE *fd = fopen("results/firingrate", "w+");
        spike_t *curr_spike = head;
        int count[n];
        int   i   = 0;
        float sum = 0.0;
        for(i = 0; i < n; i++){
            count[i] = 0;
        }
        while(curr_spike){
            i = curr_spike->index;
            count[i] += 1;
            curr_spike = curr_spike->next;        
        }
        for(i = 0; i < n; i++){
            sum += count[i];
            fprintf(fd, "%d\n", count[i]);
        }
        fclose(fd);
        return sum/n;
    }
    else return 0.0;
}

float calc_isi(spike_t *head, int n, int spike_cnt){
    /* 
    Calculates interspike intervals (ISIs) and writes them into a file. Returns avergage ISI.
    */
    head = sort_spikes(head, spike_cnt);
    if(head){
        FILE *fd = fopen("results/isi", "w+");
        spike_t *spike = head;
        float t0[n];
        int i, j = 0;
        float isi;
        float sum = 0.0;
        for(i = 0; i < n; i++){
            t0[i] = 0;
        }
        while(spike){
            i = spike->index;
            if (t0[i] == 0){    // first spike of neuron i
                t0[i] = spike->t;   
            }
            else{
                isi  = spike->t - t0[i];
                sum += isi;
                j++;
                fprintf(fd, "%f\n", spike->t - t0[i]);
                t0[i] = spike->t;         
            }
            spike = spike->next;        
        }
        fclose(fd);
        return sum/j;
    }
    else return 0.0;
}

void statistics(sim_t *sim){
    /*
    Calculate and save various statistics for the current simulation run.
    */
    state_t *state_mem = sim->state_mem;
    FILE *fd = fopen("results/stats.txt", "w+");
    float mean_isi = calc_isi(sim->top_spike, sim->n, sim->spike_cnt);
    float mean_firingrate = calc_firingrate(sim->top_spike, sim->n, sim->t_end - sim->t_start);
    float mean_V_m = calc_avgvoltage(state_mem, sim->n);
    float mean_g_ex = 0.0, mean_g_in = 0.0;
    int spike_cnt  = sim->spike_cnt;
    for (int i = 0; i < sim->n; i++){
        mean_V_m  += state_mem[i].V_m;
        mean_g_ex += state_mem[i].g_ex;
        mean_g_in += state_mem[i].g_in;
    }
    mean_V_m  /= sim->n;
    mean_g_ex /= sim->n;
    mean_g_in /= sim->n;

    fprintf(fd, "Number of spikes: %i\n", spike_cnt);
    fprintf(fd, "Average membrane voltage: %.2f mV\n", mean_V_m);
    fprintf(fd, "Average exc. conductance: %.2f nS\n", mean_g_ex);
    fprintf(fd, "Average inh. conductance: %.2f nS\n", mean_g_in);
    fprintf(fd, "Average interspike interval: %.2f ms\n", mean_isi);
    fprintf(fd, "Average firing rate : %.2f ms\n", mean_firingrate);
}

void calc_update(state_t *update, float *factors, float g_ex, float g_in){
    update->V_m  = 0;
    update->g_ex = g_ex;
    update->g_in = g_in;
    solve_analytic(update, factors);
    update->t_ela = 0;
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
        get_factors(sim->h - t_s, sim->factors_dt, sim->calc_factors);
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
    float   *factors_dt     = sim->factors_dt;
    state_t *state_mem      = sim->state_mem;
    state_t *buffered_state = sim->buffered_state;

    /* Standard timing */
    if (sim->interpolation == None){   
        solve_analytic(&state_mem[i], sim->factors_h);
        add_state(&state_mem[i], &buffered_state[i]);
    }

    /* Exact timing */
    else{ 
        /* Neuron emerges from refractory period this update interval */  
        if (state_mem[i].t_ela > tau_ref - h && state_mem[i].t_ela <= tau_ref){
            t_em = state_mem[i].t_ela + h - tau_ref;
            get_factors(t_em, factors_dt, sim->calc_factors);
            solve_analytic(&state_mem[i], factors_dt);            
            state_mem[i].V_m = 0;

            get_factors(h - t_em, factors_dt, sim->calc_factors);
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

void generate_spike(sim_t *sim, int i, float t){
    state_t     *state_mem      = sim->state_mem;
    state_t     *tmp_mem        = sim->tmp_mem;
    state_t     *update         = sim->update;
    state_buf_t *state_buf      = sim->state_buf;
    float       *factors_dt     = sim->factors_dt;
    synapse_t   **synapses      = sim->synapses;
    float h = sim->h;
    float t_s, delay;
    float y0 = tmp_mem[i].V_m;
    float yh = state_mem[i].V_m;
    float y0_dot, yh_dot;
    int target;

    /* Calculate spike time */
    switch (sim->interpolation){
        case None:
            t_s = h;
            break;
        
        case Linear:
            t_s = linear_int(y0, yh, V_th, h);
            break;
        
        case Quadratic:
            y0_dot = voltage_deriv(y0, tmp_mem[i].g_ex, tmp_mem[i].g_in);
            t_s = quadratic_int(y0, y0_dot, yh, V_th, h);
            break;

        case Cubic: // not implemented yet
            y0_dot = voltage_deriv(y0, tmp_mem[i].g_ex, tmp_mem[i].g_in);
            yh_dot = voltage_deriv(yh, state_mem[i].g_ex, state_mem[i].g_in);
            t_s = cubic_int(y0, y0_dot, yh, yh_dot, V_th, h);
            break;
    }
    fprintf(fd_raster, "%f %i\n", t + t_s, i);

    /* Calculate update */
    // Note: if a single neuron can have both excitatory and inhibitory synapses
    // the calculation of the update has to be done in the loop below
    get_factors(h - t_s, factors_dt, sim->calc_factors);    
    if (synapses[i][0].type == Inhibitory){
        calc_update(update, factors_dt, 0, synapses[i][0].weight);
    }
    else{
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
        memcpy(sim->tmp_mem, state_mem, sizeof(state_t) * sim->n);    // copy of state variables at time t
        for(int i = 0; i < n; i++){
            subthreshold_dynamics(sim, i);

            /* Neuron spikes */
            if (state_mem[i].V_m >= V_th){
                generate_spike(sim, i, t);
            }
            
            /* Neuron is in refractory period - clamp to resting potential */
            if (state_mem[i].t_ela < tau_ref){
                state_mem[i].V_m = E_rest;
            }
        }
        /* Write state variables to output files */
        fprintf(fd_voltage, "%f %f\n", t, state_mem[990].V_m);
        fprintf(fd_conductance, "%f %f %f\n", t, state_mem[990].g_ex, -state_mem[990].g_in);
    }
}

int main(void){
    sim_t *sim = setup_sim();
    clock_t start, end;
    float simulation_time;
    start = clock();
    simulation_loop(sim);
    end = clock();
    
    simulation_time = ((float) (end - start) / CLOCKS_PER_SEC);
    printf("Simulation time: %f s\n", simulation_time);
    statistics(sim);
    clear_sim(sim);
    
    // For measurement of execution time
    /*M1 = clock();
    M2 = clock();
    P1 = ((float) (M2 - M1) / CLOCKS_PER_SEC);
    printf("Performance time: %f s\n", P1);*/    
}