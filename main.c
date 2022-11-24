#include <string.h>
#include <stdio.h>
#include "gffio.h"

int main_read(int argc, char *argv[])
{
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
	if (strcmp(argv[1], "view") == 0) return main_read(argc-1, argv+1);
	else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
