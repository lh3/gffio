#ifndef KAGRAPH_H
#define KAGRAPH_H

#include <stdint.h>

typedef struct { // directed graph with forward arcs only
	uint32_t n_v, n_arc;
	uint32_t *nei;
	uint64_t *idx;
} kag_gfor_t;

void kag_gfor_destroy(kag_gfor_t *g);

/**
 * Compute Strongly Connected Components (SCCs) with Tarjan's algorithm
 */
uint64_t *kag_scc_tarjan(const kag_gfor_t *g);

#endif
