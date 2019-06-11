#include <math.h>
#include "va_model.h"

float state_size = 4;
float E_rest     = 0.0;
float E_L        = -60.0;
float E_ex       = 0.0;
float E_in       = -80.0;
float dg_ex      = 0.27;
float dg_in      = 4.5;
float dg_stim    = 0.27;
float E_avg      = -60.0;
float I_inj      = 120.0;
float R_L        = 0.1;
float V_th       = 10.0;
float tau_ref    = 5.0;
float tau_ex     = 1.0 / 5.0;
float tau_in     = 1.0 / 10.0;
float tau_stim   = 1.0 / 10.0;
float tau_L      = 1.0 / 20.0;

/* Some pre-calculated constants for integration */
int c1 = 2;     // = (E_avg - E_ex) * R_L * tau_L / (tau_L - tau_ex)
int c2 = -2;    // = (E_avg - E_in) * R_L * tau_L / (tau_L - tau_in)
int c3 = 11;    // = I_inj * R_L
int c4 = 1;     // = (tau_L^2 + tau_ex*tau_in + tau_ex*tau_L + tau_in*tau_L)/((tau_ex + tau_L)*(tau_in + tau_L))

void solve_analytic(state_t *state, float *factors){
    /*
    Exponential integration of state variables in 'state'. 'factors' are calculated in
    'calc_factors()' based on the time interval.
    */
    state->t_ela += factors[0];
    state->g_ex   = factors[1] * state->g_ex;
    state->g_in   = factors[2] * state->g_in;
    state->V_m    = factors[3] * state->V_m + factors[4] * state->g_ex + factors[5] * state->g_in + factors[6];
}

void calc_factors(float dt, float *factors){
    /*
    Calculate the integration factors for time interval dt and write them to factors.
    */
    factors[0] = dt;
    factors[1] = exp(-tau_ex * dt);
    factors[2] = exp(-tau_in * dt);
    factors[3] = exp(-tau_L  * dt);
    factors[4] = (factors[3] - factors[1]) * c1;
    factors[5] = (factors[3] - factors[2]) * c2;
    factors[6] = exp(tau_L * dt) * c3 * (c4 - factors[3]);
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

float voltage_deriv(float t, float V_m, float g_ex, float g_in){
    // TODO: check and simplify
    return -V_m/tau_L * exp(-t/tau_L)\
    + g_ex * R_L * tau_ex * (E_ex - E_rest)/(tau_L - tau_ex) * (exp(-t/tau_ex)/tau_ex - exp(-t/tau_L)/tau_L)\
    + g_in * R_L * tau_in * (E_in - E_rest)/(tau_L - tau_in) * (exp(-t/tau_in)/tau_in - exp(-t/tau_L)/tau_L);
}

void state_add(state_t *s1, state_t *s2, state_t *res){
    /*
    Adds all state variables: *res = *s1 + *s2
    */
    res->t_ela = s1->t_ela + s2->t_ela;
    res->V_m = s1->V_m + s2->V_m;
    res->g_ex = s1->g_ex + s2->g_ex;
    res->g_in = s1->g_in + s2->g_in;
}

void state_sub(state_t *s1, state_t *s2, state_t *res){
    /*
    Adds all state variables: *res = *s1 - *s2
    */
    res->t_ela = s1->t_ela - s2->t_ela;
    res->V_m = s1->V_m - s2->V_m;
    res->g_ex = s1->g_ex - s2->g_ex;
    res->g_in = s1->g_in - s2->g_in;
}