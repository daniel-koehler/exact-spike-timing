#ifndef va_model_h
#define va_model_h

/* State variables of Vogels Abbot neuron model */
typedef struct state_t {
    float t_ela;
    float V_m;
    float g_ex;
    float g_in;
} state_t;

extern state_t ZERO_STATE;

/* Parameters of Vogels Abbott neuron model */
extern float state_size;
extern float E_rest;
extern float E_L;
extern float E_ex;
extern float E_in;
extern float dg_ex;
extern float dg_in;
extern float dg_stim;
extern float E_avg;
extern float I_inj;
extern float R_L;
extern float V_th;
extern float tau_ref;
extern float tau_ex;
extern float tau_in;
extern float tau_stim;
extern float tau_L;

/* Precalculated constants to speed up calculation of integration factors */
extern int c1;      // = (E_avg - E_ex) * R_L * tau_L / (tau_L - tau_ex)
extern int c2;      // = (E_avg - E_in) * R_L * tau_L / (tau_L - tau_in)
extern int c3;      // = I_inj * R_L
extern int c4;      // = (tau_L^2 + tau_ex*tau_in + tau_ex*tau_L + tau_in*tau_L)/((tau_ex + tau_L)*(tau_in + tau_L))

/* Exponential integration */
void solve_analytic(state_t *state, float *factors);
void calc_factors(float dt, float *factors);

/* Interpolation for exact spike timing */
float linear_int(float y0, float yh, float yth, float h);
float quadratic_int(float y0, float y0_dot, float yh, float yth, float h);
float cubic_int(float y0, float y0_dot, float yh, float yh_dot, float yth, float h);
float voltage_deriv(float t, float V_m, float g_ex, float g_in);

/* Auxiliary functions for handling states*/
void state_add(state_t *s1, state_t *s2);
void state_sub(state_t *s1, state_t *s2);
#endif