#include <string.h>
#include <stdio.h>
#include "gffio.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	gio_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c;
	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gffio view [options] <in.gff>\n");
		return 1;
	}
	gff = gio_read(argv[o.ind]);
	gio_write(0, gff, GIO_FMT_GFF3);
	gio_destroy(gff);
	return 0;
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: gffio <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view      read GTF/GFF3\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage(stderr);
	if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
