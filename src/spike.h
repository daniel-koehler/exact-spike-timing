#ifndef spike_h
#define spike_h

typedef struct spike_t {
    float   t;
    int     index;
    struct spike_t *next;
} spike_t;

int compare_spikes(const void *p1, const void *p2);
spike_t *sort_spikes(spike_t *head, int spike_cnt);
void print_spike(spike_t *spike);
void print_spikes(spike_t *head);
spike_t *append_spike(spike_t *spike, float t, int idx);
void free_spikes(spike_t *top_input);

#endif