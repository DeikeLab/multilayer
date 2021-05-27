/**
	 Drawing snapshots of Navier-Stokes solver result. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "tag.h"

double snapshot_time;
int main (int argc, char * argv[])
{
	if (argc > 1)
		snapshot_time = atof(argv[1]);
	origin (-L0/2, -L0/2, -L0/2);
	periodic (right);
	periodic (front);
	N = 32;
  char targetname[100], imagename[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
  	fprintf(ferr, "Not restored!\n");
  	return 1;
  }
	restore (targetname);
	sprintf (imagename, "ux%g", snapshot_time);
	view (width = 800, height = 600, theta = pi/4, phi = pi/6, fov = 20, bg = {255,255,255});
  /* for (double x = -2*L0; x <= L0; x += L0) */
  /*   translate (x) { */
	/* 		isosurface("f", 0.5, color = "u.x", min = -0.15, max = 0.6); */
	/* 	} */
	//isosurface("f", 0.5, color = "u.x", min = -0.15, max = 0.6);
	draw_vof ("f", larger = 1.2, edges = false, color = "u.x", min = -0.15, max = 0.6);
  char s[80];
  sprintf (s, "t = %.2f T0", snapshot_time);
  draw_string (s, size = 80);
	save ("ux.ppm");
	/* { */
	/* 	static FILE * fp = fopen (imagename, "w"); */
	/* 	save (fp = fp); */
	/* } */
}

