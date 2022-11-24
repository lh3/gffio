#include <string.h>
#include <math.h>
#include <zlib.h>
#include "gio-priv.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

static void gio_parse_attr(gio_feat_t *f, char *str)
{
}

static void gio_parse_feat(gio_gff_t *gff, char *str)
{
	int32_t i;
	char *p, *q;
	gio_feat_t *f;
	if (gff->n_feat == gff->m_feat)
		GIO_EXPAND(gff->feat, gff->m_feat);
	f = &gff->feat[gff->n_feat++];
	for (p = q = str, i = 0;; ++p) {
		if (*p == '\t' || *p == 0) {
			int32_t c = *p;
			*p = 0;
			if (i == 0) { // contig name
				f->ctg = gio_dict_put(gff->dict, q);
			} else if (i == 1) {
				f->src = gio_dict_put(gff->dict, q);
			} else if (i == 2) { // source
				f->feat_ori = gio_dict_put(gff->dict, q);
			} else if (i == 3) { // start
				f->start = atol(q) - 1;
			} else if (i == 4) { // end
				f->end = atol(q);
			} else if (i == 5) { // score
				char *r;
				f->score = strcmp(q, ".") == 0? nan("1") : strtod(q, &r);
			} else if (i == 6) { // strand
				f->strand = *q == '+'? 1 : *q == '-'? -1 : 0;
			} else if (i == 7) { // frame
				f->frame = *q >= '0' && *q <= '9'? atoi(q) : -1;
			} else if (i == 8) { // attributes
				gio_parse_attr(f, q);
			}
			q = p + 1, ++i;
			if (c == 0 || i == 9) break;
		}
	}
}

gio_gff_t *gio_read_ks(kstream_t *ks)
{
	gio_gff_t *gff = 0;
	kstring_t str = {0,0,0};
	int dret;
	int64_t lineoff = 0;
	GIO_CALLOC(gff, 1);
	gff->dict = gio_dict_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		if (str.s[0] == '#') {
			gio_comm_t *p;
			if (gff->n_comm == gff->m_comm)
				GIO_EXPAND(gff->comm, gff->m_comm);
			p = &gff->comm[gff->n_comm++];
			p->lineoff = lineoff;
			GIO_CALLOC(p->line, str.l + 1);
			memcpy(p->line, str.s, str.l + 1);
		} else {
			gio_parse_feat(gff, str.s);
		}
		++lineoff;
	}
	return gff;
}

gio_gff_t *gio_read(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	gio_gff_t *gff;
	fp = fn || strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	gff = gio_read_ks(ks);
	ks_destroy(ks);
	gzclose(fp);
	return gff;
}

void gio_destroy(gio_gff_t *gff)
{
}
