#ifndef MINIGFF_H
#define MINIGFF_H

#include <stdint.h>

#define MGF_FEAT_NONE      0
#define MGF_FEAT_GENE      1
#define MGF_FEAT_MRNA      2
#define MGF_FEAT_EXON      3
#define MGF_FEAT_CDS       4
#define MGF_FEAT_INTRON    5
#define MGF_FEAT_START     6
#define MGF_FEAT_STOP      7

#define MGF_FMT_GFF3       1
#define MGF_FMT_GTF        2
#define MGF_FMT_BED12L     3
#define MGF_FMT_BED12S     4
#define MGF_FMT_BED_EXON   5
#define MGF_FMT_BED_CDS    6
#define MGF_FMT_BED_INTRON 7
#define MGF_FMT_FA_MRNA    8
#define MGF_FMT_FA_CDS     9
#define MGF_FMT_FA_PROTEIN 10

typedef struct {
	const char *key, *val;
} mgf_attr_t;

typedef struct mgf_feat_s {
	int32_t feat;
	int32_t strand:16, frame:16;
	int32_t n_child, n_parent;
	int32_t n_attr, m_attr;
	int64_t st, en;
	double score; 
	int64_t lineoff;
	const char *ctg, *src, *feat_ori, *id;
	mgf_attr_t *attr;
	struct mgf_feat_s **parent, **child;
} mgf_feat_t;

typedef struct {
	int64_t lineoff;
	char *line;
} mgf_comm_t;

struct mgf_dict_s;
typedef struct mgf_dict_s mgf_dict_t;

typedef struct {
	int64_t n_feat, m_feat;
	int64_t n_comm, m_comm;
	mgf_feat_t *feat;
	const mgf_feat_t **feat_view;
	mgf_comm_t *comm;
	mgf_dict_t *dict;
	void *dict_id;
} mgf_gff_t;

typedef struct {
	int32_t n_seq, m_seq;
	int64_t *len;
	char **name, **seq;
	void *h;
} mgf_seqs_t;

typedef struct {
	int64_t st, en;
} mgf_intv_t;

struct mgf_qbuf_s;
typedef struct mgf_qbuf_s mgf_qbuf_t;

typedef struct {
	int32_t n_exon, m_exon;
	int32_t err:16, strand:2, has_cds:2, has_start:2, has_stop:2, cds_fixed:2, frame:6;
	int32_t m_name;
	int64_t st, en, st_cds, en_cds;
	char *name;
	const char *ctg;
	mgf_intv_t *exon;
} mgf_mrna_t;

extern int mgf_verbose;

#ifdef __cplusplus
extern "C" {
#endif

mgf_gff_t *mgf_read(const char *fn);
void mgf_destroy(mgf_gff_t *gff);
void mgf_write(const char *fn, const mgf_gff_t *gff, int32_t fmt);
void mgf_write_list(const char *fn, const mgf_gff_t *gff, int32_t fmt, int32_t n, char **list);
void mgf_write_fasta(const char *fn, const mgf_gff_t *gff, const mgf_seqs_t *seq, int32_t fmt);

void mgf_write_feat(char **str, int32_t *len, int32_t *cap, const mgf_gff_t *gff, const mgf_feat_t *f, int32_t fmt);
void mgf_group(mgf_gff_t *gff);

const mgf_feat_t *mgf_get_by_id(const mgf_gff_t *gff, const char *id);
mgf_qbuf_t *mgf_qbuf_init(const mgf_gff_t *gff);
void mgf_qbuf_destroy(mgf_qbuf_t *b);
const mgf_feat_t **mgf_descend(mgf_qbuf_t *b, const mgf_feat_t *f, int32_t *n);
int32_t mgf_extract_seq(const mgf_gff_t *gff, const mgf_seqs_t *seq, const mgf_mrna_t *t, int32_t fmt, char **str_, int32_t *cap_);

char **mgf_read_list(const char *o, int *n_);

mgf_seqs_t *mgf_seqs_read(const char *fn);
void mgf_seqs_destroy(mgf_seqs_t *s);

#ifdef __cplusplus
}
#endif

#endif
