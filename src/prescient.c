#include <stdlib.h>
#include <math.h>
#include "prescient.h"

void setup_sim(sim_t *sim){
    /*
    Implementation-specific setup of simulation structure.
    */
    int   slots      = (sim->max_delay - sim->min_delay)/sim->h + 1;
    sim->top_spike   = NULL;
    sim->next_spike  = NULL;
    sim->spike_cnt   = 0;

    sim->state_buf = create_buffer(slots, sim->n);
    sim->buffered_state = (state_t *) malloc(sizeof(state_t) * sim->n);
}

void clear_sim(sim_t *sim){
    free_buffer(sim->state_buf);
    if (sim->buffered_state)    free(sim->buffered_state);
    if (sim->top_input)         free_spikes(sim->top_input);
    if (sim->top_spike)         free_spikes(sim->top_spike);
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

    buf_read_all(sim->state_buf, sim->buffered_state);
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

void generate_spike(sim_t *sim, int i, float t, FILE *fd_raster){
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