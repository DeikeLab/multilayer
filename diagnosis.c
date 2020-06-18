#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

#include "utils.h"
/* vector get_vector_field (char * prefix, int layerindex) { */
/*    char s[10]; */
/*    sprintf (s, "%d", layerindex); */
/* 	 strcat (prefix, s); */
/* 	 fprintf (stderr, "%s", prefix); */
/*    return lookup_vector (prefix); */
/* } */
/* scalar get_scalar_field (char * prefix, int layerindex) { */
/*    char s[10]; */
/*    sprintf (s, "%d", layerindex); */
/* 	 strcat (prefix, s); */
/*    return lookup_field (prefix); */
/* } */

vector get_vector_field (int layerindex) {
   char s[10];
   sprintf (s, "u%d", layerindex);
   return lookup_vector (s);
}
scalar get_scalar_field (int layerindex) {
   char s[10];
   sprintf (s, "w%d", layerindex);
   return lookup_field (s);
}

double snapshot_time = 0;
int layer_number = 59;
scalar SLOPE[];

void output ()
{
	char filename[100];
	vector targetv = get_vector_field (layer_number);
	sprintf (filename, "ux_%g_layer%d", snapshot_time, layer_number);
	FILE * fp1 = fopen (filename, "w");
	output_matrix (targetv.x, fp1, linear = true);
  fclose (fp1);
	scalar targets = get_scalar_field (layer_number);
	sprintf (filename, "w_%g_layer%d", snapshot_time, layer_number);
	FILE * fp2 = fopen (filename, "w");
	output_matrix (targets, fp2, linear = true);
	fclose (fp2);
	foreach () {
		SLOPE[] = (eta[1,0]-eta[-1,0])/(2.*Delta);
	}
	sprintf (filename, "eta_%g_layer%d", snapshot_time, layer_number);
	FILE * fp3 = fopen (filename, "w");
	output_matrix (eta, fp3, linear = true);
	fclose (fp3);
	sprintf (filename, "slope_%g_layer%d", snapshot_time, layer_number);
	FILE * fp4 = fopen (filename, "w");
	output_matrix (SLOPE, fp4, linear = true);
	fclose (fp4);

}

int main (int argc, char * argv[])
{
	N = 256;
	nl = 60;
	if (argc > 1)
		snapshot_time = atof(argv[1]);
	if (argc > 2)
		layer_number = atoi(argv[2]);
	if (argc > 3) // Optional argument if total layer number and N are different
		nl = atoi(argv[3]);
	if (argc > 4)
		N = atoi(argv[4]);
	run ();
}

output_movie () {
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0) 
    fprintf (fp, "set term pngcairo font ',9' size 1024,768;"
	     "unset key\n"
	     "set pm3d interpolate 4,4 lighting specular 0.6\n"
	     "set zrange [-5:5]\n"
	     "set cbrange [-2:2]\n"
	     "set xlabel 'x'\n"
	     "set ylabel 'y'\n"
	     );
  fprintf (fp,
	   "set output 'plot%05d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot '-' u 1:2:3:4 w pm3d\n",
	   i, t);
	char s[10];
	sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  output_field ({eta,u_temp.x}, fp, linear = true);
  fprintf (fp, "e\n\n");
  fflush (fp);
	char filename1[50], filename2[50];
	sprintf (filename1, "./surface/eta_matrix_%g", t);
	sprintf (filename2, "./surface/eta_%g", t);
	FILE * feta1 = fopen (filename1, "w");
	// Might need to change to mpi function later
	output_matrix (eta, feta1, N, linear = true);
	fclose (feta1);
	FILE * feta2 = fopen (filename2, "w");
	output_field ({eta}, feta2, linear = true);
	fclose (feta2);
}
event restore_dump (i = 0)
{
	char targetname[100];
	sprintf (targetname, "dump%g", snapshot_time);
	if (!restore (targetname)) {
  	fprintf(ferr, "Not restored!\n");
  	return 1;
  }
  restore (targetname);
	output_movie ();
	return 1;
}
