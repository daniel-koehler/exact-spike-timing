#ifndef va_model_h
#define va_model_h

/* State variables of Vogels Abbot neuron model */
typedef struct state {
    float t_ela;
    float V_m;
    float g_ex;
    float g_in;
} state_t;

typedef enum interpolation {
    NONE,
    LINEAR,
    QUADRATIC,
    CUBIC
} interpolation_t;

typedef struct lut{
    float **values;
    float t_max;
    int   denom;
    int   num_entries;
    int   idx_factor;   // used to calculate table index for given time
} lut_t;

extern state_t ZERO_STATE;

/* Parameters of Vogels Abbott neuron model */
extern const float state_size;
extern const float E_rest;
extern const float E_L;
extern const float E_ex;
extern const float E_in;
extern const float dg_ex;
extern const float dg_in;
extern const float dg_stim;
extern const float I_inj;
extern const float R_L;
extern const float V_th;
extern const float tau_ref;
extern const float tau_ex;
extern const float tau_in;
extern const float tau_stim;
extern const float tau_L;

/* Number of constants and factors for exponential integration */
extern const int NUM_CONSTANTS;
extern const int NUM_FACTORS;

/* Exponential integration */
void solve_analytic(state_t *state, float *factors);
void calc_constants(float *constants);
void calc_factors(float dt, float *factors);
void generate_lut(float h, int denom);
void lookup(float t, float * factors);
void free_lut();

/* Interpolation for exact spike timing */
float linear_int(float y0, float yh, float yth, float h);
float quadratic_int(float y0, float y0_dot, float yh, float yth, float h);
float cubic_int(float y0, float y0_dot, float yh, float yh_dot, float yth, float h);
float voltage_deriv(float V_m, float g_ex, float g_in);

/* Auxiliary functions for handling states*/
void add_state(state_t *s1, state_t *s2);
void sub_state(state_t *s1, state_t *s2);
void print_state(state_t *s);
void print_factors(float *fac);
#endif