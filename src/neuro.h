/*
Contains various typedefs used in sim.c .
*/

#ifndef neuro_h

#define neuro_h

typedef struct state_t {
    float t_ela;
    float V_m;
    float g_ex;
    float g_in;
} state_t;

typedef struct neuron_model_t
{
    int state_size;     // number of state variables
    float E_rest;
    float E_L;
    float E_ex;
    float E_in;
    float d_gex;
    float d_gin;
    float d_gstim;
    float E_avg;
    float I_inj;
    float R_L;
    float V_th;
    float tau_ref;
    float tau_ex;
    float tau_in;
    float tau_stim;
    float tau_L;
} neuron_model_t;

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

void state_add(state_t *s1, state_t *s2, state_t *res);
void state_sub(state_t *s1, state_t *s2, state_t *res);

#endif