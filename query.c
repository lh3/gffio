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

void mgf_mrna_free(mgf_mrna_t *t)
{
	free(t->name); free(t->exon);
}

int32_t mgf_mrna_gen(mgf_qbuf_t *b, const mgf_gff_t *gff, const mgf_feat_t *f, mgf_mrna_t *t)
{
	int32_t i, j, n_fs, n_exon, n_cds;
	const mgf_feat_t **fs;
	const char *s;
	kstring_t str = {0,0,0};

	if (f == 0 || f->feat != MGF_FEAT_MRNA) return -1; // only looking at mRNA or transcript
	fs = mgf_descend(b, f, &n_fs);
	if (n_fs == 0) return -1; // not retrieved anything
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
	t->has_cds = (n_cds > 0), t->frame = 0;
	if (t->n_exon == 0) { // TODO: this can be relaxed
		if (mgf_verbose >= 2)
			fprintf(stderr, "[W::%s] no exon associated with a transcript\n", __func__);
		return -4;
	}
	if (t->n_exon > t->m_exon) {
		t->m_exon = t->n_exon;
		kroundup32(t->m_exon);
		MGF_REALLOC(t->exon, t->m_exon);
	}

	if (t->has_cds) {
		const mgf_feat_t *c0 = 0, *c1 = 0;
		for (i = 0; i < n_fs; ++i) {
			if (fs[i]->feat != MGF_FEAT_CDS) continue;
			if (c0 == 0 || c1 == 0) c0 = c1 = fs[i];
			if (fs[i]->st < c0->st) c0 = fs[i];
			if (fs[i]->en > c1->en) c1 = fs[i];
		}
		if (f->strand > 0 && c0->frame >= 0)
			t->frame = c0->frame;
		else if (f->strand < 0 && c1->frame >= 0)
			t->frame = c1->frame;
	}

	t->st = f->st, t->en = f->en;
	t->st_cds = INT64_MAX, t->en_cds = INT64_MIN;
	t->strand = f->strand, t->ctg = f->ctg;
	for (i = j = 0; i < n_fs; ++i) {
		const mgf_feat_t *e = fs[i];
		if (e->feat == MGF_FEAT_CDS) {
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
	if (!t->has_cds) t->st_cds = t->st, t->en_cds = t->en;
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
	mgf_sprintf_lite(&str, "%d:", t->frame);
	s = mgf_attr_find(gff, f, "transcript_type");
	if (s == 0)
		s = mgf_attr_find(gff, f, "transcript_biotype");
	if (s) mgf_sprintf_lite(&str, "%s", s);
	t->m_name = str.m, t->name = str.s;

	if ((t->st < f->st || t->en > f->en) && mgf_verbose >= 2)
		fprintf(stderr, "[W::%s] exon coordinates beyond transcript coordinates at %s\n", __func__, t->name);
	return 0;
}

/*********************
 * Extract sequences *
 *********************/

static char *mgf_codon_std = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX";
					// 01234567890123456789012345678901234567890123456789012345678901234
					// KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX <- this is the AGCT order

static unsigned char mgf_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static char mgf_comp_tab[128] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static void mgf_revcomp(int32_t len, char *seq)
{
	int32_t i;
	for (i = 0; i < len>>1; ++i) {
		uint8_t t = seq[len - 1 - i];
		seq[len - 1 - i] = (uint8_t)seq[i] >= 128? seq[i] : mgf_comp_tab[(uint8_t)seq[i]];
		seq[i] = t >= 128? t : mgf_comp_tab[t];
	}
	if (len&1) seq[i] = (uint8_t)seq[i] >= 128? seq[i] : mgf_comp_tab[(uint8_t)seq[i]];
}

int32_t mgf_extract_seq(const mgf_gff_t *gff, const mgf_seqs_t *seq, const mgf_mrna_t *t, int32_t fmt, char **str_, int32_t *cap_)
{
	int32_t j, len = 0, cap = *cap_, seq_id;
	char *str = *str_;
	if ((fmt == MGF_FMT_FA_CDS || fmt == MGF_FMT_FA_PROTEIN) && t->has_cds == 0) return -1; // no CDS
	seq_id = mgf_id_get(seq->h, t->ctg);
	if (seq_id < 0) return -1; // contig name not found; TODO: add a warning
	if (t->en > seq->len[seq_id]) return -1; // beyond the end of ctg
	if (fmt == MGF_FMT_FA_MRNA) {
		for (j = 0, len = 0; j < t->n_exon; ++j)
			len += t->exon[j].en - t->exon[j].st;
	} else if (fmt == MGF_FMT_FA_CDS || fmt == MGF_FMT_FA_PROTEIN) {
		for (j = 0, len = 0; j < t->n_exon; ++j) {
			int64_t st, en;
			st = t->exon[j].st > t->st_cds? t->exon[j].st : t->st_cds;
			en = t->exon[j].en < t->en_cds? t->exon[j].en : t->en_cds;
			if (st >= en) continue;
			len += en - st;
		}
	} else abort();
	if (fmt == MGF_FMT_FA_PROTEIN && len % 3 != 0 && mgf_verbose >= 2)
		fprintf(stderr, "[W::%s] the length of %s CDS is not a multiple of 3\n", __func__, t->name);
	if (len + 1 > cap) {
		cap = len + 1;
		kroundup32(cap);
		MGF_REALLOC(str, cap);
	}
	if (fmt == MGF_FMT_FA_MRNA) {
		for (j = 0, len = 0; j < t->n_exon; ++j) {
			memcpy(&str[len], &seq->seq[seq_id][t->exon[j].st], t->exon[j].en - t->exon[j].st);
			len += t->exon[j].en - t->exon[j].st;
		}
	} else if (fmt == MGF_FMT_FA_CDS || fmt == MGF_FMT_FA_PROTEIN) {
		for (j = 0, len = 0; j < t->n_exon; ++j) {
			int64_t st, en;
			st = t->exon[j].st > t->st_cds? t->exon[j].st : t->st_cds;
			en = t->exon[j].en < t->en_cds? t->exon[j].en : t->en_cds;
			if (st >= en) continue;
			memcpy(&str[len], &seq->seq[seq_id][st], en - st);
			len += en - st;
		}
	}
	if (t->strand < 0) mgf_revcomp(len, str);
	if (fmt == MGF_FMT_FA_PROTEIN) {
		int32_t j, k;
		for (j = t->frame, k = 0; j + 2 < len; j += 3) {
			int32_t c0 = mgf_nt4_table[(uint8_t)str[j]];
			int32_t c1 = mgf_nt4_table[(uint8_t)str[j+1]];
			int32_t c2 = mgf_nt4_table[(uint8_t)str[j+2]];
			str[k++] = c0>3 || c1>3 || c2>3? 'X' : mgf_codon_std[c0<<4|c1<<2|c2];
		}
		len = k;
	}
	str[len] = 0;
	*str_ = str, *cap_ = cap;
	return len;
}
