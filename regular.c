#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "mgf-priv.h"
#include "kagraph.h"

int mgf_verbose = 3;

void mgf_attr_append(mgf_gff_t *g, mgf_feat_t *f, const char *key, const char *val)
{
	if (key != 0 && val != 0) {
		if (f->n_attr == f->m_attr)
			MGF_EXPAND(f->attr, f->m_attr);
		f->attr[f->n_attr].key = mgf_dict_put(g->dict, key);
		f->attr[f->n_attr].val = mgf_dict_put(g->dict, val);
		f->n_attr++;
	}
}

const char *mgf_attr_find(const mgf_gff_t *g, const mgf_feat_t *f, const char *key)
{
	const char *key2;
	int32_t i;
	key2 = mgf_dict_get(g->dict, key);
	if (key2 == 0) return 0; // the key doesn't exist
	for (i = 0; i < f->n_attr; ++i)
		if (f->attr[i].key == key2)
			return f->attr[i].val;
	return 0;
}

static void mgf_add_id_name(mgf_gff_t *gff, mgf_feat_t *f, const char *s_id, const char *s_xid, const char *id_name, const char *id_xkey)
{
	int32_t j, j_id = -1, j_xid = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_id) j_id = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_id >= 0 && j_xid < 0)
		mgf_attr_append(gff, f, id_xkey, f->attr[j_id].val);
	else if (j_xid >= 0 && j_id < 0)
		mgf_attr_append(gff, f, id_name, f->attr[j_xid].val);
	else if (j_xid < 0 && j_id < 0) {
		if (mgf_verbose >= 2)
			fprintf(stderr, "WARNING: missing ID for '%s %ld %ld'\n", f->ctg, (long)f->st+1, (long)f->en);
	}
}

static void mgf_add_parent(mgf_gff_t *gff, mgf_feat_t *f, const char *s_xid, const char *s_par)
{
	int32_t j, j_xid = -1, j_par = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_par) j_par = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_par >= 0 || j_xid < 0) return;
	mgf_attr_append(gff, f, "Parent", f->attr[j_xid].val);
}

static void mgf_check_id(const mgf_gff_t *gff)
{
	int32_t i;
	if (mgf_verbose < 4) return;
	for (i = 0; i < gff->n_feat; ++i) {
		mgf_feat_t *f = &gff->feat[i];
		if (f->id == 0 && f->n_parent == 0)
			fprintf(stderr, "WARNING: missing both ID and Parent at '%s %ld %ld'\n", f->ctg, (long)f->st+1, (long)f->en);
	}
}

static void mgf_build_id_dict(mgf_gff_t *gff)
{
	int32_t i, t, n_dup = 0;
	gff->dict_id = mgf_id_init();
	for (i = 0; i < gff->n_feat; ++i) {
		mgf_feat_t *f = &gff->feat[i];
		if (f->id) {
			t = mgf_id_put(gff->dict_id, f->id, i);
			if (!t) ++n_dup;
			if (!t && mgf_verbose >= 4)
				fprintf(stderr, "WARNING: duplicated ID '%s'\n", f->id);
		}
	}
	if (n_dup > 0 && mgf_verbose >= 2)
		fprintf(stderr, "WARNING: there are duplicated IDs\n");
}

void mgf_label(mgf_gff_t *gff)
{
	int32_t i;
	const char *s_gene, *s_trans, *s_mrna, *s_exon, *s_cds, *s_start, *s_stop;
	const char *s_id, *s_par, *s_name, *s_gid, *s_tid, *s_gname, *s_tname;

	s_gene  = mgf_dict_get(gff->dict, "gene");
	s_trans = mgf_dict_get(gff->dict, "transcript");
	s_mrna  = mgf_dict_get(gff->dict, "mRNA");
	s_exon  = mgf_dict_get(gff->dict, "exon");
	s_cds   = mgf_dict_get(gff->dict, "CDS");
	s_start = mgf_dict_get(gff->dict, "start_codon");
	s_stop  = mgf_dict_get(gff->dict, "stop_codon");
	s_id    = mgf_dict_get(gff->dict, "ID");
	s_par   = mgf_dict_get(gff->dict, "Parent");
	s_name  = mgf_dict_get(gff->dict, "Name");
	s_gid   = mgf_dict_get(gff->dict, "gene_id");
	s_tid   = mgf_dict_get(gff->dict, "transcript_id");
	s_gname = mgf_dict_get(gff->dict, "gene_name");
	s_tname = mgf_dict_get(gff->dict, "transcript_name");

	for (i = 0; i < gff->n_feat; ++i) {
		mgf_feat_t *f = &gff->feat[i];
		if (f->feat_ori == s_gene) {
			f->feat = MGF_FEAT_GENE;
			mgf_add_id_name(gff, f, s_id, s_gid, "ID", "gene_id");
			mgf_add_id_name(gff, f, s_name, s_gname, "Name", "gene_name");
		} else if (f->feat_ori == s_trans || f->feat_ori == s_mrna) {
			f->feat = MGF_FEAT_TRANS;
			mgf_add_id_name(gff, f, s_id, s_tid, "ID", "transcript_id");
			mgf_add_id_name(gff, f, s_name, s_tname, "Name", "transcript_name");
			mgf_add_parent(gff, f, s_gid, s_par);
		} else if (f->feat_ori == s_exon) {
			f->feat = MGF_FEAT_EXON;
			mgf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_cds) {
			f->feat = MGF_FEAT_CDS;
			mgf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_start) {
			f->feat = MGF_FEAT_START;
			mgf_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_stop) {
			f->feat = MGF_FEAT_STOP;
			mgf_add_parent(gff, f, s_tid, s_par);
		}
		f->id = mgf_attr_find(gff, f, "ID");
	}
	mgf_check_id(gff);
	mgf_build_id_dict(gff);
}

void mgf_connect(mgf_gff_t *gff)
{
	int32_t i, k;
	const char *s_par = 0, *par;
	s_par = mgf_dict_get(gff->dict, "Parent");
	if (s_par == 0) return;
	for (i = 0; i < gff->n_feat; ++i) {
		mgf_feat_t *f = &gff->feat[i];
		f->n_child = f->n_parent = 0;
	}
	for (i = 0; i < gff->n_feat; ++i) { // count children
		mgf_feat_t *f = &gff->feat[i];
		par = mgf_attr_find(gff, f, "Parent");
		if (par) {
			k = mgf_id_get(gff->dict_id, par);
			if (k >= 0) { // TODO: support multiple parents (low priority)
				gff->feat[k].n_child++;
				f->n_parent = 1;
			} else if (mgf_verbose >= 2)
				fprintf(stderr, "WARNING: failed to find Parent ID '%s'\n", par);
		}
	}
	for (i = 0; i < gff->n_feat; ++i) { // allocate
		mgf_feat_t *f = &gff->feat[i];
		if (f->n_parent + f->n_child > 0) {
			MGF_CALLOC(f->parent, f->n_parent + f->n_child);
			f->child = f->parent + f->n_parent;
			f->n_parent = f->n_child = 0;
		}
	}
	for (i = 0; i < gff->n_feat; ++i) { // population ::child and ::parent
		mgf_feat_t *f = &gff->feat[i];
		par = mgf_attr_find(gff, f, "Parent");
		if (par) {
			k = mgf_id_get(gff->dict_id, par);
			if (k >= 0) { // TODO: support multiple parents (low priority)
				gff->feat[k].child[gff->feat[k].n_child++] = f;
				f->parent[f->n_parent++] = &gff->feat[k];
			}
		}
	}
}

// Convert gff grpah to TJ graph. This is not necessary but will make TJ more general
static kag_gfor_t *mgf_gff2tj(const mgf_gff_t *gff)
{
	int32_t i, j;
	uint64_t n_arc;
	kag_gfor_t *g;
	for (i = 0, n_arc = 0; i < gff->n_feat; ++i)
		n_arc += gff->feat[i].n_child;
	assert(n_arc < UINT32_MAX);
	MGF_CALLOC(g, 1);
	g->n_v = gff->n_feat;
	MGF_CALLOC(g->idx, g->n_v);
	MGF_CALLOC(g->nei, n_arc);
	for (i = 0, g->n_arc = 0; i < gff->n_feat; ++i) {
		const mgf_feat_t *f = &gff->feat[i];
		g->idx[i] = (uint64_t)g->n_arc << 32 | f->n_child;
		for (j = 0; j < f->n_child; ++j)
			g->nei[g->n_arc++] = f->child[j] - gff->feat;
	}
	return g;
}

const mgf_feat_t **mgf_toposort(const mgf_gff_t *gff)
{
	uint32_t v;
	uint64_t *b;
	kag_gfor_t *g;
	const mgf_feat_t **a;
	g = mgf_gff2tj(gff);
	b = kag_scc_tarjan(g);
	MGF_CALLOC(a, gff->n_feat);
	for (v = 0; v < gff->n_feat; ++v)
		a[v] = &gff->feat[(uint32_t)b[v]];
	free(b);
	kag_gfor_destroy(g);
	return a;
}

void mgf_group(mgf_gff_t *gff)
{
	gff->feat_view = mgf_toposort(gff);
}
