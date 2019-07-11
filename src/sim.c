#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "statebuf.h"
#include "va_model.h"
#include "spike.h"

#ifdef CANONICAL
#include "canonical.h"
#else
#include "prescient.h"
#endif

clock_t M1, M2;
float P1;

void print_state_mem(sim_t *sim){
    state_t *state_mem = sim->state_mem;
    int   n          = sim->n;
    for(int i=0; i<n; i++){
        printf("Index %i: t_ela = %.2f, V_m = %.2f, g_ex = %.2f, g_in = %.2f\n", i, state_mem[i].t_ela,\
        state_mem[i].V_m, state_mem[i].g_ex, state_mem[i].g_in);
    }
}

void setup_parameters(sim_t *sim){
    /* 
    Setup all important parameters of the simulation.
    */
    sim->n       = 4000;
    sim->h       = 0.05;
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
    sim->calc_factors  = Lookup;

    sim->n_ex      = floor(sim->n * sim->r_ei / (sim->r_ei + 1));
    sim->n_in      = sim->n - sim->n_ex;
    sim->n_syn     = (sim->n - 1) * sim->p_conn;
}

void initialize_state_mem(sim_t *sim){
    /*
    Initialize state variables of all neurons
    */
    int i;
    state_t *state_mem;
	float h          = sim->h;
	int   n          = sim->n;

    sim->state_mem = state_mem = (state_t *) malloc(sizeof(state_t)*n);

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

sim_t *setup(void){
    /*
    Set simulation parameters, initalize memory, generate stimulus
    */

    sim_t *sim = (sim_t *) malloc(sizeof(sim_t));

    setup_parameters(sim);
    setup_sim(sim);
    create_events(sim);

    create_network(sim);

    initialize_state_mem(sim);

    if(sim->calc_factors == Lookup){
        generate_lut(sim->h, 100);      // generates lookup table with 100+1 entries for time interval [0 h]
    }
    sim->factors_h  = (float *) malloc(sizeof(float) * NUM_FACTORS);
    sim->factors_dt = (float *) malloc(sizeof(float) * NUM_FACTORS);
    sim->spike_cnt  = 0;
    sim->top_spike  = NULL;
    calc_factors(sim->h, sim->factors_h);
    
    return(sim);

}

void clean_up(sim_t *sim){
    /*
    Free allocated memory
    */
    if (sim){
        clear_sim(sim);     
        if (sim->synapses){
            for(int i = 0; i < sim->n; i++){
                free(sim->synapses[i]);
            }
            free(sim->synapses);
        }
        if (sim->state_mem)         free(sim->state_mem);
        if (sim->factors_h)         free(sim->factors_h);
        if (sim->factors_dt)        free(sim->factors_dt);
        free_lut();
        if (sim->top_spike)         free_spikes(sim->top_spike);   
        
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

    fclose(fd);
}



int main(void){
    sim_t *sim = setup();
    clock_t start, end;
    float simulation_time;
    start = clock();
    simulation_loop(sim);
    end = clock();
    
    simulation_time = ((float) (end - start) / CLOCKS_PER_SEC);
    printf("Simulation time: %f s\n", simulation_time);
    statistics(sim);
    clean_up(sim);
    
    
    // For measurement of execution time
    /*M1 = clock();
    M2 = clock();
    P1 = ((float) (M2 - M1) / CLOCKS_PER_SEC);
    printf("Performance time: %f s\n", P1);*/    
}
