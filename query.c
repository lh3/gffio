#include <stdio.h>
#include "mgf-priv.h"

struct mgf_qbuf_s {
	int32_t ns, ms, nf, mf;
	const mgf_feat_t **stack, **rst;
	int8_t *flag;
	const mgf_gff_t *gff;
};

const mgf_feat_t *mgf_get_by_id(const mgf_gff_t *gff, const char *id)
{
	int32_t k;
	k = mgf_id_get(gff->dict_id, id);
	return k < 0? 0 : &gff->feat[k];
}

mgf_qbuf_t *mgf_qbuf_init(const mgf_gff_t *gff)
{
	mgf_qbuf_t *b;
	MGF_CALLOC(b, 1);
	b->gff = gff;
	MGF_CALLOC(b->flag, b->gff->n_feat);
	return b;
}

void mgf_qbuf_destroy(mgf_qbuf_t *b)
{
	free(b->stack); free(b->rst); free(b->flag); free(b);
}

const mgf_feat_t **mgf_descend(mgf_qbuf_t *b, const mgf_feat_t *f, int32_t *n)
{
	int32_t i;
	b->nf = 0;
	MGF_PUSH_BACK(b->ns, b->ms, b->stack, f);
	while (b->ns > 0) {
		f = b->stack[--b->ns];
		MGF_PUSH_BACK(b->nf, b->mf, b->rst, f);
		for (i = f->n_child - 1; i >= 0; --i) {
			const mgf_feat_t *g = f->child[i];
			int32_t v = g - b->gff->feat;
			if (!b->flag[v])
				MGF_PUSH_BACK(b->ns, b->ms, b->stack, g);
		}
	}
	for (i = 0; i < b->nf; ++i)
		b->flag[b->rst[i] - b->gff->feat] = 0;
	*n = b->nf;
	return b->rst;
}
