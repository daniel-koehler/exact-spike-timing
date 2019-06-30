#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "statebuf.h"
#include "va_model.h"
#include "spike.h"

#ifdef CANONICAL

#else // Prescient implementation
#include "prescient.h"
#endif


clock_t M1, M2;
float P1;

#define frand() ((double) rand() / (RAND_MAX + 1.0))



/* FILE objects used to save simulation results */
FILE *fd_raster;
FILE *fd_voltage;
FILE *fd_conductance;
FILE *fd_isi;
FILE *fd_firingrate;
FILE *fd_avgvoltage;
FILE *fd_stats;

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
    fd_isi         = fopen("results/isi", "w+");
    fd_firingrate  = fopen("results/firingrate", "w+");
    fd_avgvoltage  = fopen("results/avgvoltage", "w+");
    fd_stats       = fopen("results/stats.txt", "w+");
    if ( !fd_raster || !fd_voltage || !fd_isi || !fd_stats || !fd_conductance || !fd_firingrate || !fd_avgvoltage){
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
        if (fd_isi)                 fclose(fd_isi);
        if (fd_firingrate)          fclose(fd_firingrate);
        if (fd_avgvoltage)          fclose(fd_avgvoltage);
        if (fd_stats)               fclose(fd_stats);
    }
    free(sim);
    
}

float calc_avgvoltage(state_t *state_mem, int n, FILE *fd){
    int sum = 0;
    for (int i = 0; i < n; i++){
        fprintf(fd, "%f\n", state_mem[i].V_m);
        sum += state_mem[i].V_m;
    }
    return sum/n;
}

float calc_firingrate(spike_t *head, int n, float t_span, FILE *fd){
    /*
    Calculates the average firing rate of each neuron and returns average firing rate of all neurons.
    */
    if(head){
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
        return sum/n;
    }
    else return 0.0;
}

float calc_isi(spike_t *head, int n, int spike_cnt, FILE *fd){
    /* 
    Calculates interspike intervals (ISIs) and writes them into a file. Returns avergage ISI.
    */
    head = sort_spikes(head, spike_cnt);
    if(head){
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
        return sum/j;
    }
    else return 0.0;
}

void statistics(sim_t *sim){
    /*
    Calculate and save various statistics for the current simulation run.
    */
    state_t *state_mem = sim->state_mem;

    float mean_isi = calc_isi(sim->top_spike, sim->n, sim->spike_cnt, fd_isi);
    float mean_firingrate = calc_firingrate(sim->top_spike, sim->n, sim->t_end - sim->t_start, fd_firingrate);
    float mean_V_m = calc_avgvoltage(state_mem, sim->n, fd_avgvoltage);
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

    fprintf(fd_stats, "Number of spikes: %i\n", spike_cnt);
    fprintf(fd_stats, "Average membrane voltage: %.2f mV\n", mean_V_m);
    fprintf(fd_stats, "Average exc. conductance: %.2f nS\n", mean_g_ex);
    fprintf(fd_stats, "Average inh. conductance: %.2f nS\n", mean_g_in);
    fprintf(fd_stats, "Average interspike interval: %.2f ms\n", mean_isi);
    fprintf(fd_stats, "Average firing rate : %.2f ms\n", mean_firingrate);
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
                process_spike(sim, i, t, fd_raster);
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
