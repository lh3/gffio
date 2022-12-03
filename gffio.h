#ifndef GFFIO_H
#define GFFIO_H

#include <stdint.h>

#define GF_FEAT_NONE      0
#define GF_FEAT_GENE      1
#define GF_FEAT_MRNA      2
#define GF_FEAT_EXON      3
#define GF_FEAT_CDS       4
#define GF_FEAT_INTRON    5
#define GF_FEAT_START     6
#define GF_FEAT_STOP      7

#define GF_FMT_GFF3       1
#define GF_FMT_GTF        2
#define GF_FMT_BED12L     3
#define GF_FMT_BED12S     4
#define GF_FMT_BED_EXON   5
#define GF_FMT_BED_CDS    6
#define GF_FMT_BED_INTRON 7
#define GF_FMT_FA_MRNA    8
#define GF_FMT_FA_CDS     9
#define GF_FMT_FA_PROTEIN 10

typedef struct {
	const char *key, *val;
} gf_attr_t;

typedef struct gf_feat_s {
	int32_t feat;
	int32_t strand:16, frame:10, has_id:2, has_parent:2, has_name:2;
	int32_t n_child, n_parent;
	int32_t n_attr, m_attr;
	int64_t st, en;
	double score; 
	int64_t lineoff;
	const char *ctg, *src, *feat_ori;
	gf_attr_t *attr;
	const char *id, *parent1, *name, *biotype;
	struct gf_feat_s **parent, **child;
} gf_feat_t;

typedef struct {
	int64_t lineoff;
	char *line;
} gf_comm_t;

struct gf_dict_s;
typedef struct gf_dict_s gf_dict_t;

typedef struct {
	int64_t n_feat, m_feat, n_feat_view;
	int64_t n_comm, m_comm;
	gf_feat_t *feat;
	const gf_feat_t **feat_view;
	gf_comm_t *comm;
	gf_dict_t *dict;
	void *dict_id;
} gf_gff_t;

typedef struct {
	int32_t n_seq, m_seq;
	int64_t *len;
	char **name, **seq;
	void *h;
} gf_seqs_t;

typedef struct {
	int64_t st, en;
} gf_intv_t;

struct gf_qbuf_s;
typedef struct gf_qbuf_s gf_qbuf_t;

typedef struct {
	int32_t n_exon, m_exon;
	int32_t err:16, strand:2, has_cds:2, has_start:2, has_stop:2, cds_fixed:2, frame:6;
	int32_t m_name;
	int64_t st, en, st_cds, en_cds;
	char *name;
	const char *ctg;
	gf_intv_t *exon;
} gf_mrna_t;

extern int gf_verbose;

#ifdef __cplusplus
extern "C" {
#endif

gf_gff_t *gf_read(const char *fn);
void gf_destroy(gf_gff_t *gff);
void gf_write(const char *fn, const gf_gff_t *gff, int32_t fmt);
void gf_write_list(const char *fn, const gf_gff_t *gff, int32_t fmt, int32_t n, char **list);
void gf_write_fasta(const char *fn, const gf_gff_t *gff, const gf_seqs_t *seq, int32_t fmt);

void gf_write_feat(char **str, int32_t *len, int32_t *cap, const gf_gff_t *gff, const gf_feat_t *f, int32_t fmt);
void gf_group(gf_gff_t *gff);

const gf_feat_t *gf_get_by_id(const gf_gff_t *gff, const char *id);
gf_qbuf_t *gf_qbuf_init(const gf_gff_t *gff);
void gf_qbuf_destroy(gf_qbuf_t *b);
const gf_feat_t **gf_descend(gf_qbuf_t *b, const gf_feat_t *f, int32_t *n);
int32_t gf_extract_seq(const gf_gff_t *gff, const gf_seqs_t *seq, const gf_mrna_t *t, int32_t fmt, char **str_, int32_t *cap_);
void gf_mrna_choose_long(gf_gff_t *gff);

char **gf_read_list(const char *o, int *n_);

gf_seqs_t *gf_seqs_read(const char *fn);
void gf_seqs_destroy(gf_seqs_t *s);

#ifdef __cplusplus
}
#endif

#endif
