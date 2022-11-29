#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "gfpriv.h"
#include "kagraph.h"

int gf_verbose = 3;

void gf_attr_append(gf_gff_t *g, gf_feat_t *f, const char *key, const char *val)
{
	if (key != 0 && val != 0) {
		if (f->n_attr == f->m_attr)
			GF_EXPAND(f->attr, f->m_attr);
		f->attr[f->n_attr].key = gf_dict_put(g->dict, key);
		f->attr[f->n_attr].val = gf_dict_put(g->dict, val);
		f->n_attr++;
	}
}

const char *gf_attr_find(const gf_gff_t *g, const gf_feat_t *f, const char *key)
{
	const char *key2;
	int32_t i;
	key2 = gf_dict_get(g->dict, key);
	if (key2 == 0) return 0; // the key doesn't exist
	for (i = 0; i < f->n_attr; ++i)
		if (f->attr[i].key == key2)
			return f->attr[i].val;
	return 0;
}

static void gf_add_id_name(gf_gff_t *gff, gf_feat_t *f, const char *s_id, const char *s_xid, const char *id_name, const char *id_xkey)
{
	int32_t j, j_id = -1, j_xid = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_id) j_id = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_id >= 0 && j_xid < 0)
		gf_attr_append(gff, f, id_xkey, f->attr[j_id].val);
	else if (j_xid >= 0 && j_id < 0)
		gf_attr_append(gff, f, id_name, f->attr[j_xid].val);
	else if (j_xid < 0 && j_id < 0) {
		if (gf_verbose >= 2)
			fprintf(stderr, "WARNING: missing ID for '%s %ld %ld'\n", f->ctg, (long)f->st+1, (long)f->en);
	}
}

static void gf_add_parent(gf_gff_t *gff, gf_feat_t *f, const char *s_xid, const char *s_par)
{
	int32_t j, j_xid = -1, j_par = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_par) j_par = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_par >= 0 || j_xid < 0) return;
	gf_attr_append(gff, f, "Parent", f->attr[j_xid].val);
}

static void gf_check_id(const gf_gff_t *gff)
{
	int32_t i;
	if (gf_verbose < 4) return;
	for (i = 0; i < gff->n_feat; ++i) {
		gf_feat_t *f = &gff->feat[i];
		if (f->id == 0 && f->n_parent == 0)
			fprintf(stderr, "WARNING: missing both ID and Parent at '%s %ld %ld'\n", f->ctg, (long)f->st+1, (long)f->en);
	}
}

static void gf_build_id_dict(gf_gff_t *gff)
{
	int32_t i, t, n_dup = 0;
	gff->dict_id = gf_id_init();
	for (i = 0; i < gff->n_feat; ++i) {
		gf_feat_t *f = &gff->feat[i];
		if (f->id) {
			t = gf_id_put(gff->dict_id, f->id, i);
			if (!t) ++n_dup;
			if (!t && gf_verbose >= 4)
				fprintf(stderr, "WARNING: duplicated ID '%s'\n", f->id);
		}
	}
	if (n_dup > 0 && gf_verbose >= 2)
		fprintf(stderr, "WARNING: there are duplicated IDs. Please use -v4 to see details.\n");
}

void gf_label(gf_gff_t *gff)
{
	int32_t i;
	const char *s_gene, *s_trans, *s_mrna, *s_exon, *s_cds, *s_start, *s_stop;
	const char *s_id, *s_par, *s_name, *s_gid, *s_tid, *s_gname, *s_tname;

	s_gene  = gf_dict_get(gff->dict, "gene");
	s_trans = gf_dict_get(gff->dict, "transcript");
	s_mrna  = gf_dict_get(gff->dict, "mRNA");
	s_exon  = gf_dict_get(gff->dict, "exon");
	s_cds   = gf_dict_get(gff->dict, "CDS");
	s_start = gf_dict_get(gff->dict, "start_codon");
	s_stop  = gf_dict_get(gff->dict, "stop_codon");
	s_id    = gf_dict_get(gff->dict, "ID");
	s_par   = gf_dict_get(gff->dict, "Parent");
	s_name  = gf_dict_get(gff->dict, "Name");
	s_gid   = gf_dict_get(gff->dict, "gene_id");
	s_tid   = gf_dict_get(gff->dict, "transcript_id");
	s_gname = gf_dict_get(gff->dict, "gene_name");
	s_tname = gf_dict_get(gff->dict, "transcript_name");

	for (i = 0; i < gff->n_feat; ++i) {
		gf_feat_t *f = &gff->feat[i];
		if (f->feat_ori == s_gene) {
			f->feat = GF_FEAT_GENE;
			gf_add_id_name(gff, f, s_id, s_gid, "ID", "gene_id");
			gf_add_id_name(gff, f, s_name, s_gname, "Name", "gene_name");
		} else if (f->feat_ori == s_trans || f->feat_ori == s_mrna) {
			f->feat = GF_FEAT_MRNA;
			gf_add_id_name(gff, f, s_id, s_tid, "ID", "transcript_id");
			gf_add_id_name(gff, f, s_name, s_tname, "Name", "transcript_name");
			gf_add_parent(gff, f, s_gid, s_par);
		} else if (f->feat_ori == s_exon) {
			f->feat = GF_FEAT_EXON;
			gf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_cds) {
			f->feat = GF_FEAT_CDS;
			gf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_start) {
			f->feat = GF_FEAT_START;
			gf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_stop) {
			f->feat = GF_FEAT_STOP;
			gf_add_parent(gff, f, s_tid, s_par);
		}
		f->id = gf_attr_find(gff, f, "ID");
	}
	gf_check_id(gff);
	gf_build_id_dict(gff);
}

void gf_connect(gf_gff_t *gff)
{
	int32_t i, k;
	const char *s_par = 0, *par;
	s_par = gf_dict_get(gff->dict, "Parent");
	if (s_par == 0) return;
	for (i = 0; i < gff->n_feat; ++i) {
		gf_feat_t *f = &gff->feat[i];
		f->n_child = f->n_parent = 0;
	}
	for (i = 0; i < gff->n_feat; ++i) { // count children
		gf_feat_t *f = &gff->feat[i];
		par = gf_attr_find(gff, f, "Parent");
		if (par) {
			k = gf_id_get(gff->dict_id, par);
			if (k >= 0) { // TODO: support multiple parents (low priority)
				gff->feat[k].n_child++;
				f->n_parent = 1;
			} else if (gf_verbose >= 2)
				fprintf(stderr, "WARNING: failed to find Parent ID '%s'\n", par);
		}
	}
	for (i = 0; i < gff->n_feat; ++i) { // allocate
		gf_feat_t *f = &gff->feat[i];
		if (f->n_parent + f->n_child > 0) {
			GF_CALLOC(f->parent, f->n_parent + f->n_child);
			f->child = f->parent + f->n_parent;
			f->n_parent = f->n_child = 0;
		}
	}
	for (i = 0; i < gff->n_feat; ++i) { // population ::child and ::parent
		gf_feat_t *f = &gff->feat[i];
		par = gf_attr_find(gff, f, "Parent");
		if (par) {
			k = gf_id_get(gff->dict_id, par);
			if (k >= 0) { // TODO: support multiple parents (low priority)
				gff->feat[k].child[gff->feat[k].n_child++] = f;
				f->parent[f->n_parent++] = &gff->feat[k];
			}
		}
	}
}

// Convert gff grpah to TJ graph. This is not necessary but will make TJ more general
static kag_gfor_t *gf_gff2tj(const gf_gff_t *gff)
{
	int32_t i, j;
	uint64_t n_arc;
	kag_gfor_t *g;
	for (i = 0, n_arc = 0; i < gff->n_feat; ++i)
		n_arc += gff->feat[i].n_child;
	assert(n_arc < UINT32_MAX);
	GF_CALLOC(g, 1);
	g->n_v = gff->n_feat;
	GF_CALLOC(g->idx, g->n_v);
	GF_CALLOC(g->nei, n_arc);
	for (i = 0, g->n_arc = 0; i < gff->n_feat; ++i) {
		const gf_feat_t *f = &gff->feat[i];
		g->idx[i] = (uint64_t)g->n_arc << 32 | f->n_child;
		for (j = 0; j < f->n_child; ++j)
			g->nei[g->n_arc++] = f->child[j] - gff->feat;
	}
	return g;
}

const gf_feat_t **gf_toposort(const gf_gff_t *gff)
{
	uint32_t v;
	uint64_t *b;
	kag_gfor_t *g;
	const gf_feat_t **a;
	g = gf_gff2tj(gff);
	b = kag_scc_tarjan(g);
	GF_CALLOC(a, gff->n_feat);
	for (v = 0; v < gff->n_feat; ++v)
		a[v] = &gff->feat[(uint32_t)b[v]];
	free(b);
	kag_gfor_destroy(g);
	return a;
}

void gf_group(gf_gff_t *gff)
{
	gff->feat_view = gf_toposort(gff);
}
