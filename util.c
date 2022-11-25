#include <stdlib.h>
#include "gio-priv.h"

void gio_attr_append(gio_gff_t *g, gio_feat_t *f, const char *key, const char *val)
{
	if (key != 0 && val != 0) {
		if (f->n_attr == f->m_attr)
			GIO_EXPAND(f->attr, f->m_attr);
		f->attr[f->n_attr].key = gio_dict_put(g->dict, key);
		f->attr[f->n_attr].val = gio_dict_put(g->dict, val);
		f->n_attr++;
	}
}

const char *gio_attr_find(const gio_gff_t *g, const gio_feat_t *f, const char *key)
{
	const char *key2;
	int32_t i;
	key2 = gio_dict_get(g->dict, key);
	if (key2 == 0) return 0; // the key doesn't exist
	for (i = 0; i < f->n_attr; ++i)
		if (f->attr[i].key == key2)
			return f->attr[i].val;
	return 0;
}

static void gio_add_id_name(gio_gff_t *gff, gio_feat_t *f, const char *s_id, const char *s_xid, const char *id_name, const char *id_xkey)
{
	int32_t j, j_id = -1, j_xid = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_id) j_id = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_id >= 0 && j_xid < 0)
		gio_attr_append(gff, f, id_xkey, f->attr[j_id].val);
	else if (j_xid >= 0 && j_id < 0)
		gio_attr_append(gff, f, id_name, f->attr[j_xid].val);
}

static void gio_add_parent(gio_gff_t *gff, gio_feat_t *f, const char *s_xid, const char *s_par)
{
	int32_t j, j_xid = -1, j_par = -1;
	for (j = 0; j < f->n_attr; ++j) {
		if (f->attr[j].key == s_par) j_par = j;
		else if (f->attr[j].key == s_xid) j_xid = j;
	}
	if (j_par >= 0 || j_xid < 0) return;
	gio_attr_append(gff, f, "Parent", f->attr[j_xid].val);
}

void gio_label(gio_gff_t *gff)
{
	int32_t i;
	const char *s_gene, *s_trans, *s_mrna, *s_exon, *s_cds, *s_start, *s_stop;
	const char *s_id, *s_par, *s_name, *s_gid, *s_tid, *s_gname, *s_tname;

	s_gene  = gio_dict_get(gff->dict, "gene");
	s_trans = gio_dict_get(gff->dict, "transcript");
	s_mrna  = gio_dict_get(gff->dict, "mRNA");
	s_exon  = gio_dict_get(gff->dict, "exon");
	s_cds   = gio_dict_get(gff->dict, "CDS");
	s_start = gio_dict_get(gff->dict, "start_codon");
	s_stop  = gio_dict_get(gff->dict, "stop_codon");
	s_id    = gio_dict_get(gff->dict, "ID");
	s_par   = gio_dict_get(gff->dict, "Parent");
	s_name  = gio_dict_get(gff->dict, "Name");
	s_gid   = gio_dict_get(gff->dict, "gene_id");
	s_tid   = gio_dict_get(gff->dict, "transcript_id");
	s_gname = gio_dict_get(gff->dict, "gene_name");
	s_tname = gio_dict_get(gff->dict, "transcript_name");

	for (i = 0; i < gff->n_feat; ++i) {
		gio_feat_t *f = &gff->feat[i];
		if (f->feat_ori == s_gene) {
			f->feat = GIO_FEAT_GENE;
			gio_add_id_name(gff, f, s_id, s_gid, "ID", "gene_id");
			gio_add_id_name(gff, f, s_name, s_gname, "Name", "gene_name");
		} else if (f->feat_ori == s_trans || f->feat_ori == s_mrna) {
			f->feat = GIO_FEAT_TRANS;
			gio_add_id_name(gff, f, s_id, s_tid, "ID", "transcript_id");
			gio_add_id_name(gff, f, s_name, s_tname, "Name", "transcript_name");
			gio_add_parent(gff, f, s_gid, s_par);
		} else if (f->feat_ori == s_exon) {
			f->feat = GIO_FEAT_EXON;
			gio_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_cds) {
			f->feat = GIO_FEAT_CDS;
			gio_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_start) {
			f->feat = GIO_FEAT_START;
			gio_add_parent(gff, f, s_tid, s_par);
		} else if (f->feat_ori == s_stop) {
			f->feat = GIO_FEAT_STOP;
			gio_add_parent(gff, f, s_tid, s_par);
		}
	}
}
