#include <stdlib.h>
#include <math.h>
#include "canonical.h"

/* FILE objects used to save simulation results */
FILE *fd_raster;
FILE *fd_voltage;
FILE *fd_conductance;

void setup_sim(sim_t *sim){
    /*
    Implementation-specific setup of simulation structure.
    */

    sim->tmp_state = (state_t *) malloc(sizeof(state_t));
    sim->spike_influence = (state_t *) malloc(sizeof(state_t));
    /* Open output files */
    fd_raster      = fopen("results/raster", "w+");
    fd_voltage     = fopen("results/voltage", "w+");
    fd_conductance = fopen("results/conductance", "w+");
    if ( !fd_raster || !fd_voltage || !fd_conductance){
        printf("Could not open files to save results. Exit.\n");
        exit(1);
    }
}

void clear_sim(sim_t *sim){
    if (sim->spikes){
        for(int i = 0; i < sim->n; i++){
            free_spikes(sim->spikes[i]);
        }
        free(sim->spikes);
    }
    if (sim->tmp_state)         free(sim->tmp_state);
    if (sim->spike_influence)   free(sim->spike_influence);
    if (fd_raster)              fclose(fd_raster);
    if (fd_voltage)             fclose(fd_voltage);
    if (fd_conductance)         fclose(fd_conductance);
}

void create_events(sim_t *sim){
    /*
    Creates a poisson distributed input spike train for each neuron i and stores it to sim->spike[i].
    */
    int n = sim->n;
    int i;
    float dt;
    int     spike_cnt    = 0;
    spike_t *new_input   = NULL;
    spike_t *next_input  = NULL;
    spike_t **spikes = (spike_t **) malloc(sizeof(spike_t*) * n);

    for(i = 0; i < n; i++){
        spikes[i] = NULL;
        next_input = NULL;
        for(float t = sim->t_start; t <= sim->t_input; ){
            /* Add up exponentially distributed intervals to generate Poisson spike train */
            dt  = -log(1.0 - frand()) * sim->t_avg;
            t  += dt;
            if (t > sim->t_input) break;

            /* Create a new input spike */
            new_input         = (spike_t *) malloc(sizeof(spike_t));
            new_input->t      = t;
            new_input->index  = i;
            new_input->weight = dg_stim;
            new_input->next   = NULL;

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
}

void calc_update(sim_t *sim, int i, float t, float dt){
    state_t *state_mem = sim->state_mem;
    get_factors(dt, sim->factors_dt, sim->calc_factors);
    solve_analytic(&state_mem[i], sim->factors_dt);
    /* Neuron in refractory period - clamp to resting potential */
    
    if (state_mem[i].t_ela < tau_ref){
        state_mem[i].V_m = E_rest;
    }
    /* Neuron spikes */
    else if (state_mem[i].V_m >= V_th){
        generate_spike(sim, i, t, dt, sim->tmp_state);
    }
}


void neuron_dynamics(sim_t *sim, int i, float t){
    spike_t **spikes   = sim->spikes;
    spike_t *next_spike;
    state_t *state_mem = sim->state_mem;

    state_t *tmp       = sim->tmp_state;
    state_t *spike_influence = sim->spike_influence;
    float h     = sim->h;
    float t_s   = 0.0;      // Offset of current spike
    float t_pre = 0.0;      // Offset of previous spike
    float interval;         // Time interval between subsequent spikes
    float t_em;         
    bool emerging = false;

    /* Does neuron i emerge from refractory period within (t, t+h] ? */
    if (state_mem[i].t_ela > tau_ref - h && state_mem[i].t_ela <= tau_ref){
                t_em = state_mem[i].t_ela + h - tau_ref;
                emerging = true;
    }

    /* Go through all spikes reaching neuron i in current update interval (t, t+h] */
    while(spikes[i] && spikes[i]->t <= t+h){
        t_s  = spikes[i]->t - t;
        
        *tmp = state_mem[i];        // store state variables at beginning of interval, in case it is needed for interpolation
        if (emerging && t_s > t_em){
            /* Interval (t_pre, t_em]: neuron is still in refractory period */
            emerging = false;
            interval = t_em - t_pre;
            get_factors(interval, sim->factors_dt, sim->calc_factors);
            solve_analytic(&state_mem[i], sim->factors_dt);
            state_mem[i].V_m = E_rest;
            /* Interval (t_em, t_s]: neuron has left refractory period, Vm is not clamped to resting potential anymore */
            t_pre = t_em;
        }

        /* Integrate to next incoming spike */
        interval = t_s - t_pre;
        calc_update(sim, i, t, interval);

        /* Add influence of incoming spike */
        spike_influence->t_ela = 0.0;
        spike_influence->V_m   = E_rest;
        if(spikes[i]->weight >= 0.0){
            spike_influence->g_in = 0.0;
            spike_influence->g_ex = spikes[i]->weight;
        }
        else{
            spike_influence->g_ex = 0.0;
            spike_influence->g_in = -spikes[i]->weight;
        }

        add_state(&state_mem[i], spike_influence);

        /* Go to next spike */
        t_pre = t_s;
        next_spike = spikes[i]->next;
        free(spikes[i]);
        spikes[i] = next_spike;
    }

    /* Integrate to end of the time step */
    interval = h - t_s;
    calc_update(sim, i, t, interval);
    
}

void generate_spike(sim_t *sim, int i, float t, float dt, state_t *state0){
    /*
    Exact spike time of neuron i is interpolated using points (0, V_m(0)) and (dt, V_m(dt)) and 
    derivatives of the voltage for higher-order interpolations respectively. The state variables for time t
    are passed in 'state0'.
    */
    state_t     *state_mem      = sim->state_mem;
    synapse_t   **synapses      = sim->synapses;
    float t_s;
    float y0 = state0->V_m;
    float yh = state_mem[i].V_m;
    float y0_dot, yh_dot;
    
    spike_t *spike;
    float weight;
    int target;
    float tspike;

    /* Calculate spike time */
    switch (sim->interpolation){
        case None:
            t_s = t + dt;
            break;
        
        case Linear:
            t_s = t + linear_int(y0, yh, V_th, dt);
            break;
        
        case Quadratic:
            y0_dot = voltage_deriv(y0, state0->g_ex, state0->g_in);
            t_s = t + quadratic_int(y0, y0_dot, yh, V_th, dt);
            break;

        case Cubic: // not implemented yet
            y0_dot = voltage_deriv(y0, state0->g_ex, state0->g_in);
            yh_dot = voltage_deriv(yh, state_mem[i].g_ex, state_mem[i].g_in);
            t_s = t + cubic_int(y0, y0_dot, yh, yh_dot, V_th, dt);
            break;
    }
    fprintf(fd_raster, "%f %i\n", t_s, i);

    /* For each neuron j connected to spiking neuron i add a new spike */
    for(int j = 0; j < sim->n_syn; j++){
        target = synapses[i][j].target;        
        tspike     = t_s + synapses[i][j].delay;
        if (synapses[i][j].type == Inhibitory){
            weight = -synapses[i][j].weight;
        }
        else{
            weight = synapses[i][j].weight;
        }
        spike = new_spike(target, tspike, weight);
        sortin_spike(&sim->spikes[target], spike);
    }

    state_mem[i].t_ela = 0;
    state_mem[i].V_m   = E_rest;

    /* Save spikes off all neuron together in linked list (for evaluation of statistics) */
    if (!sim->top_spike){
        sim->top_spike = sim->next_spike = (spike_t *) malloc(sizeof(spike_t));
    }
    else{
        sim->next_spike->next = (spike_t *) malloc(sizeof(spike_t));
        sim->next_spike = sim->next_spike->next;
    }
    sim->next_spike->next   = NULL;
    sim->next_spike->index  = i;
    sim->next_spike->t      = t_s;
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
    state_t *state_mem = sim->state_mem;
    for(t=t_start; t <= t_end; t+=h){
        if (t >= t_output){
            printf("t = %.0f ms\n", t);
            t_output += 0.1 * (t_end - t_start);
        }
        /* Update interval (t, t+h] */   
        for(int i = 0; i < n; i++){
            neuron_dynamics(sim, i, t);
        }
        /* Write state variables to output files */
        fprintf(fd_voltage, "%f %f\n", t, state_mem[0].V_m);
        fprintf(fd_conductance, "%f %f %f\n", t, state_mem[0].g_ex, -state_mem[0].g_in);
    }
}