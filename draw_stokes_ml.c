/**
	 Drawing snapshots of Navier-Stokes solver result. */

#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "remap_test.h"
#include "layered/perfs.h"

double snapshot_time;
int main (int argc, char * argv[])
{
	if (argc > 1)
		snapshot_time = atof(argv[1]);
	if (argc > 2)
		nl = atoi(argv[2]);
	else
		nl = 30;
	origin (-L0/2, -L0/2);
	periodic (right);
	N = 32;
  char targetname[100], imagename[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
  	fprintf(ferr, "Not restored!\n");
  	return 1;
  }
	restore (targetname);
	sprintf (imagename, "ux%g", snapshot_time);
  view (fov = 20, theta = pi/4, phi = -pi/3,  width = 800, height = 600); 
  char s[80];
  sprintf (s, "t = %.2f T0", snapshot_time);
  draw_string (s, size = 100);
  sprintf (s, "u%d.x", nl-1);
  for (double x = -1; x <= 1; x++)
    translate (x)
      squares (s, linear = true, min = -0.15, max = 0.6);
	save ("ux.ppm");
	/* { */
	/* 	static FILE * fp = fopen (imagename, "w"); */
	/* 	save (fp = fp); */
	/* } */
}

