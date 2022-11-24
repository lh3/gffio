#include <string.h>
#include <zlib.h>
#include "gffio.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

gio_gff_t *gio_read_ks(kstream_t *ks)
{
	gio_gff_t *gff = 0;
	kstring_t str = {0,0,0};
	int dret;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
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
