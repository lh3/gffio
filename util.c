#include <stdlib.h>
#include "mgf-priv.h"

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
	}
}
