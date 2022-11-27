#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "minigff.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	mgf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, to_group = 0, fmt = MGF_FMT_GFF3;
	char *id_list = 0;
	while ((c = ketopt(&o, argc, argv, 1, "gv:l:b", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
		else if (c == 'g') to_group = 1;
		else if (c == 'b') fmt = MGF_FMT_BED12;
		else if (c == 'l') id_list = o.arg;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minigff view [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -b        output in the 12-column BED format\n");
		fprintf(stderr, "  -g        group by hierarchy\n");
		fprintf(stderr, "  -l STR    extract records descended from ID list []\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", mgf_verbose);
		return 1;
	}
	gff = mgf_read(argv[o.ind]);

	if (id_list) {
		int32_t n;
		char **list;
		list = mgf_read_list(id_list, &n);
		mgf_write_list(0, gff, fmt, n, list);
	} else {
		if (to_group) mgf_group(gff);
		mgf_write(0, gff, fmt);
	}

	mgf_destroy(gff);
	return 0;
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minigff <command> <arguments>\n");
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
