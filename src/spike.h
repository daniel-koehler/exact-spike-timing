#ifndef spike_h
#define spike_h

typedef struct spike {
    int     index;
    float   t;
    float   weight;
    struct spike *next;
} spike_t;

spike_t *new_spike(int index, float t, float weight);
int compare_spikes(const void *p1, const void *p2);
spike_t *sort_spikes(spike_t *head, int spike_cnt);
void print_spike(spike_t *spike);
void print_spikes(spike_t *head);
void append_spike(spike_t **head, spike_t *spike);
void sortin_spike(spike_t **head, spike_t *spike);
void free_spikes(spike_t *top_input);

#endif