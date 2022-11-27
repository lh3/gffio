#include <stdio.h>
#include <string.h>
#include "mgf-priv.h"
#include "ksort.h"

/***************************
 * Query descendants by ID *
 ***************************/

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

/***********************
 * Get BED12-like info *
 ***********************/

#define intv_key(x) ((x).st)
KRADIX_SORT_INIT(mgf_intv, mgf_intv_t, intv_key, 8)

void mgf_mrna_init(mgf_mrna_t *t)
{
	memset(t, 0, sizeof(*t));
}

int32_t mgf_mrna_gen(mgf_qbuf_t *b, const mgf_gff_t *gff, const mgf_feat_t *f, mgf_mrna_t *t)
{
	int32_t i, j, n_fs, n_exon, n_cds;
	const mgf_feat_t **fs;
	const char *s;
	kstring_t str = {0,0,0};

	if (f == 0 || f->feat != MGF_FEAT_MRNA) return -1; // only looking at mRNA or transcript
	fs = mgf_descend(b, f, &n_fs);
	for (i = 1; i < n_fs; ++i)
		if (fs[i]->ctg != fs[0]->ctg)
			break;
	if (i < n_fs) {
		if (mgf_verbose >= 2)
			fprintf(stderr, "[W::%s] descendant features on different contigs\n", __func__);
		return -2;
	}

	n_exon = n_cds = 0;
	for (i = 0; i < n_fs; ++i) {
		const mgf_feat_t *e = fs[i];
		if (e->feat == MGF_FEAT_EXON) ++n_exon;
		if (e->feat == MGF_FEAT_CDS) ++n_cds;
	}
	if (n_exon > 0 && n_exon < n_cds) {
		if (mgf_verbose >= 2)
			fprintf(stderr, "[W::%s] more exons than CDS recards\n", __func__);
		return -3;
	}
	t->n_exon = n_exon > n_cds? n_exon : n_cds;
	t->is_cds = (n_cds > 0);
	if (t->n_exon == 0) { // TODO: this can be relaxed
		if (mgf_verbose >= 2)
			fprintf(stderr, "[W::%s] no exon associated a transcript\n", __func__);
		return -4;
	}
	if (t->n_exon > t->m_exon) {
		t->m_exon = t->n_exon;
		kroundup32(t->m_exon);
		MGF_REALLOC(t->exon, t->m_exon);
	}

	t->st_cds = t->st = INT64_MAX, t->en_cds = t->en = INT64_MIN;
	t->is_cds = 0, t->strand = f->strand;
	for (i = j = 0; i < n_fs; ++i) {
		const mgf_feat_t *e = fs[i];
		if (e->feat == MGF_FEAT_CDS) {
			t->is_cds = 1;
			t->st_cds = t->st_cds < e->st? t->st_cds : e->st;
			t->en_cds = t->en_cds > e->en? t->en_cds : e->en;
		}
		if (e->feat == MGF_FEAT_EXON || e->feat == MGF_FEAT_CDS) {
			t->st = t->st < e->st? t->st : e->st;
			t->en = t->en > e->en? t->en : e->en;
		}
		if (e->feat == MGF_FEAT_EXON || (e->feat == MGF_FEAT_CDS && n_exon == 0))
			t->exon[j].st = e->st, t->exon[j++].en = e->en;
	}
	if (!t->is_cds) t->st_cds = t->st, t->en_cds = t->en;
	assert(t->n_exon == n_exon);
	if (t->n_exon > 1)
		radix_sort_mgf_intv(t->exon, t->exon + t->n_exon);

	// generate the name string
	str.l = 0, str.m = t->m_name, str.s = t->name;
	if (f->n_parent == 1) { // get the gene ID
		const mgf_feat_t *g = f->parent[0];
		if (g->id) mgf_sprintf_lite(&str, "%s:", g->id);
		else mgf_sprintf_lite(&str, "%ld:", g->lineoff);
		s = mgf_attr_find(gff, g, "Name");
		if (s) mgf_sprintf_lite(&str, "%s:", s);
		else mgf_sprintf_lite(&str, ":");
	}
	if (f->id) mgf_sprintf_lite(&str, "%s:", f->id);
	else mgf_sprintf_lite(&str, "%ld:", f->lineoff);
	s = mgf_attr_find(gff, f, "transcript_type");
	if (s == 0)
		s = mgf_attr_find(gff, f, "transcript_biotype");
	if (s) mgf_sprintf_lite(&str, "%s", s);
	t->m_name = str.m, t->name = str.s;
	return 0;
}
