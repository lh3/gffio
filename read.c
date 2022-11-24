#include <string.h>
#include <zlib.h>
#include "gio-priv.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

gio_gff_t *gio_read_ks(kstream_t *ks)
{
	gio_gff_t *gff = 0;
	kstring_t str = {0,0,0};
	int dret;
	int64_t lineoff = 0;
	GIO_CALLOC(gff, 1);
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
