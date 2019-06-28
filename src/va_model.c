#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "va_model.h"

state_t ZERO_STATE = (state_t){
    .t_ela  = 0.0,
    .V_m    = 0.0,
    .g_ex   = 0.0,
    .g_in   = 0.0
};

const float E_rest     = 0.0;
const float E_L        = -60.0;
const float E_ex       = 0.0;
const float E_in       = -80.0;
const float dg_ex      = 0.27;
const float dg_in      = 4.5;
const float dg_stim    = 0.27;
const float I_inj      = 110.0;
const float R_L        = 0.1;
const float V_th       = 10.0;
const float tau_ref    = 5.0;
const float tau_ex     = 1.0 / 5.0;
const float tau_in     = 1.0 / 10.0;
const float tau_stim   = 1.0 / 10.0;
const float tau_L      = 1.0 / 20.0;

const float STATE_SIZE  = 4;
const int NUM_CONSTANTS = 4;
const int NUM_FACTORS   = 7;

lut_t lut;

void solve_analytic(state_t *state, float *factors){
    /*
    Exponential integration of state variables in 'state'. 'factors' are calculated in
    'calc_factors()' based on the time interval for exponential integration.
    */
    state->t_ela += factors[0];
    state->g_ex   = factors[1] * state->g_ex;
    state->g_in   = factors[2] * state->g_in;
    state->V_m    = factors[3] * state->V_m + factors[4] * state->g_ex + factors[5] * state->g_in + factors[6];
}

void calc_constants(float *constants){
    /*
    Calculates 4 float constants for integration and writes them to constants.
    */
    constants[0] = (E_L - E_ex) * R_L * tau_L / (tau_L - tau_ex);
    constants[1] = (E_L - E_in) * R_L * tau_L / (tau_L - tau_in);
    constants[2] = I_inj * R_L;
    constants[3] = (pow(tau_L, 2) + tau_ex*tau_in + tau_ex*tau_L + tau_in*tau_L)/((tau_ex + tau_L)*(tau_in + tau_L));
}

void calc_factors(float dt, float *factors, float *constants){
    /*
    Calculates 7 float integration factors for time interval dt and writes them to factors.
    */
    factors[0] = dt;
    factors[1] = exp(-tau_ex * dt);
    factors[2] = exp(-tau_in * dt);
    factors[3] = exp(-tau_L  * dt);
    factors[4] = (factors[3] - factors[1]) * constants[0];
    factors[5] = (factors[3] - factors[2]) * constants[1];
    factors[6] = exp(tau_L * dt) * constants[2] * (constants[3] - factors[3]);
}

void generate_lut(float h, int denom){
    /*
    Generates a lookup table for integration factors for the time range [0, h].
    The lookup table has a resolution of 'h / denom' and thus 'denom + 1' entries.
    */
    float t;
    int i;
    float resolution = h / denom;
    lut.t_max = h;
    lut.denom = denom;
    lut.idx_factor = denom/h;
    float *constants = (float *) malloc(sizeof(float) * NUM_CONSTANTS);
    calc_constants(constants);

    lut.values = (float **) malloc(sizeof(float *) * (denom + 1));
    for(t = 0.0, i = 0; i <= denom; t += resolution, i++){
        lut.values[i] = (float *) malloc(sizeof(float) * NUM_FACTORS);
        calc_factors(t, lut.values[i], constants);        
    }
    free(constants);
}

float *lookup(float t){
    /*
    Returns integration factors for time t stored in lut. generate_lut() has to be called before once.
    */
    if(t > lut.t_max || t < 0.0){
        return NULL;
    }
    int idx = lut.idx_factor * t;
    printf("Time: %f Index: %i, Factor: %i\n", t, idx, lut.idx_factor);
    return lut.values[idx];
}

void free_lut(){
    for(int i = 0; i <= lut.denom; i++){
        free(lut.values[i]);
    }
}

float linear_int(float y0, float yh, float yth, float h){
    /*
    Uses linear interpolation of form y(t) = (yh-y0)/h * t + y0 to determine the intersection of y(t) with yth.
    */
   return h*(y0 - yth) / (y0 - yh);
}

float quadratic_int(float y0, float y0_dot, float yh, float yth, float h){
    float denom = (yh - y0 - h * y0_dot);
    float p = pow(h,2) * y0_dot / denom;
    float q = (y0 - yth) / denom;
    float disc = pow((p/2),2) - q;
    if (disc < 0){
        return -1.0;
    }
    float res = -p/2 - sqrt(disc);
    if (res < 0){
        res = -p/2 + sqrt(disc);
    }
    return res;
}
float cubic_int(float y0, float y0_dot, float yh, float yh_dot, float yth, float h){
    return 1.0;
}

float voltage_deriv(float V_m, float g_ex, float g_in){
    return (-V_m + R_L * (g_ex * (E_ex - E_rest) + g_in * (E_in - E_rest) + I_inj))/tau_L;
}

void add_state(state_t *s1, state_t *s2){
    /*
    Adds all state variables: *s1 = *s1 + *s2
    */
    s1->t_ela = s1->t_ela + s2->t_ela;
    s1->V_m   = s1->V_m   + s2->V_m;
    s1->g_ex  = s1->g_ex  + s2->g_ex;
    s1->g_in  = s1->g_in  + s2->g_in;
}

void sub_state(state_t *s1, state_t *s2){
    /*
    Adds all state variables: *s1 = *s1 - *s2
    */
    s1->t_ela = s1->t_ela - s2->t_ela;
    s1->V_m   = s1->V_m   - s2->V_m;
    s1->g_ex  = s1->g_ex  - s2->g_ex;
    s1->g_in  = s1->g_in  - s2->g_in;
}

void print_state(state_t *s){
    printf("t_ela = %.2f, V_m = %.2f, g_ex = %.2f, g_in = %.2f\n",\
            s->t_ela,     s->V_m,     s->g_ex,     s->g_in);
}

void print_factors(float *fac){
    printf("f0 = %.6f, f1 = %.6f, f2 = %.6f, f3 = %.6f, f4 = %.6f, f5 = %.6f, f6 = %.6f\n",\
            fac[0],    fac[1],    fac[2],    fac[3],    fac[4],    fac[5],    fac[6]);
}

