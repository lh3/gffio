#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "kagraph.h"

#define kag_calloc(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define kag_realloc(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define kag_push_stack(n_, m_, a_, v_) do { \
		if ((n_) == (m_)) { \
			(m_) += ((m_)>>1) + 16; \
			kag_realloc((a_), (m_)); \
		} \
		(a_)[(n_)++] = (v_); \
	} while (0)

void kag_gfor_destroy(kag_gfor_t *g)
{
	if (g == 0) return;
	free(g->nei); free(g->idx); free(g);
}

/**************************
 * Tarjan's SCC algorithm *
 **************************/

typedef struct {
	uint32_t index, low, on_stack, n_parent;
} tj_node_t;

typedef struct {
	uint32_t index, nb;
	uint32_t nS, mS, nD, mD;
	uint32_t *S; // the Tarjan's stack
	uint64_t *D; // the DFS stack
	const kag_gfor_t *g;
	tj_node_t *a;
	uint64_t *b;
} tj_aux_t;

static void tj_scc1(tj_aux_t *p, uint32_t v0)
{
	const kag_gfor_t *g = p->g;
	kag_push_stack(p->nD, p->mD, p->D, (uint64_t)v0<<32);
	while (p->nD > 0) { // the DFS is not empty
		uint64_t x = p->D[--p->nD]; // pop the DFD stack
		uint32_t v = x>>32, i = (uint32_t)x; // i is the number of children already visited
		uint32_t nv = (uint32_t)g->idx[v], ov = g->idx[v]>>32;
		if (i == 0) { // v hasn't been visited before
			p->a[v].index = p->a[v].low = p->index++;
			p->a[v].on_stack = 1;
			kag_push_stack(p->nS, p->mS, p->S, v);
		}
		if (i == nv) { // done with v
			if (p->a[v].low == p->a[v].index) { // v is the root node
				int32_t i, j = p->nS - 1;
				while (p->S[j] != v) --j;
				for (i = p->nS - 1; i >= j; --i) {
					uint32_t w = p->S[i];
					p->b[p->nb++] = (uint64_t)v<<32 | w;
					p->a[w].on_stack = 0;
				}
				p->nS = j;
			}
			if (p->nD > 0) { // if the DFS stack is not empty, update the top element
				int32_t w = v;
				v = p->D[p->nD-1] >> 32;
				p->a[v].low = p->a[v].low < p->a[w].low? p->a[v].low : p->a[w].low;
			}
		} else { // process v's neighbor
			int32_t w = g->nei[ov + (nv - i - 1)];
			 kag_push_stack(p->nD, p->mD, p->D, (uint64_t)v<<32 | (i+1));
			if (p->a[w].index == (uint32_t)-1)
				kag_push_stack(p->nD, p->mD, p->D, (uint64_t)w<<32);
			else if (p->a[w].on_stack)
				p->a[v].low = p->a[v].low < p->a[w].low? p->a[v].low : p->a[w].low;
		}
	}
}

uint64_t *kag_scc_tarjan(const kag_gfor_t *g)
{
	uint32_t v;
	tj_aux_t aux, *p = &aux;
	memset(p, 0, sizeof(aux));
	p->g = g;
	kag_calloc(p->a, g->n_v);
	kag_calloc(p->b, g->n_v);
	for (v = 0; v < g->n_v; ++v)
		p->a[v].index = p->a[v].low = (uint32_t)-1, p->a[v].on_stack = 0, p->a[v].n_parent = 0;
	for (v = 0; v < g->n_v; ++v) {
		uint32_t j, nv = (uint32_t)g->idx[v], ov = g->idx[v]>>32;
		for (j = 0; j < nv; ++j)
			p->a[g->nei[ov + j]].n_parent++;
	}
	for (v = g->n_v - 1; v != (uint32_t)-1; --v) // start with nodes without parents
		if (p->a[v].n_parent == 0 && p->a[v].index == (uint32_t)-1)
			tj_scc1(p, v);
	for (v = g->n_v - 1; v != (uint32_t)-1; --v) // deal with the rest of nodes
		if (p->a[v].index == (uint32_t)-1)
			tj_scc1(p, v);
	assert(p->nb == g->n_v);
	for (v = 0; v < g->n_v>>1; ++v) { // reverse to get the topological order
		uint32_t t = p->b[g->n_v - v - 1];
		p->b[g->n_v - v - 1] = p->b[v], p->b[v] = t;
	}
	free(p->S);
	free(p->D);
	free(p->a);
	return p->b;
}
