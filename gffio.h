#ifndef GFFIO_H
#define GFFIO_H

#include <stdint.h>

#define GIO_FEAT_NONE    0
#define GIO_FEAT_GENE    1
#define GIO_FEAT_TRANS   2
#define GIO_FEAT_EXON    3
#define GIO_FEAT_CDS     4
#define GIO_FEAT_INTRON  5
#define GIO_FEAT_START   6
#define GIO_FEAT_STOP    7

typedef struct {
	const char *key, *val;
} gio_attr_t;

typedef struct gio_feat_s {
	int32_t feat;
	int32_t is_rev:16, frame:16;
	int32_t n_child, m_child;
	int32_t n_attr, m_attr;
	int64_t start, end;
	double score; 
	int64_t lineoff;
	char *ctg, *src, *feat_ori, *id, *name;
	gio_attr_t *attr;
	struct gio_feat_s *child;
} gio_feat_t;

typedef struct {
	int64_t lineoff;
	char *line;
} gio_comm_t;

struct gio_dict_s;
typedef struct gio_dict_s gio_dict_t;

typedef struct {
	int64_t n_feat, m_feat;
	int64_t n_comm, m_comm;
	gio_feat_t *feat;
	gio_comm_t *comm;
	gio_dict_t *dict;
} gio_gff_t;

#ifdef __cplusplus
extern "C" {
#endif

gio_gff_t *gio_read(const char *fn);
void gio_destroy(gio_gff_t *gff);

#ifdef __cplusplus
}
#endif

#endif
