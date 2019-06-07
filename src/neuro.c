#include "neuro.h"

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