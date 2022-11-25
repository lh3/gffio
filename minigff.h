#ifndef MINIGFF_H
#define MINIGFF_H

#include <stdint.h>

#define MGF_FEAT_NONE    0
#define MGF_FEAT_GENE    1
#define MGF_FEAT_TRANS   2
#define MGF_FEAT_EXON    3
#define MGF_FEAT_CDS     4
#define MGF_FEAT_INTRON  5
#define MGF_FEAT_START   6
#define MGF_FEAT_STOP    7

#define MGF_FMT_GFF3     1
#define MGF_FMT_GTF      2
#define MGF_FMT_BED12    3

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
	mgf_feat_t *feat, **feat_view;
	mgf_comm_t *comm;
	mgf_dict_t *dict;
	void *dict_id;
} mgf_gff_t;

extern int mgf_verbose;

#ifdef __cplusplus
extern "C" {
#endif

mgf_gff_t *mgf_read(const char *fn);
void mgf_destroy(mgf_gff_t *gff);
void mgf_write(const char *fn, const mgf_gff_t *gff, int32_t fmt);

#ifdef __cplusplus
}
#endif

#endif
