/**
# 3D breaking Stokes wave (multilayer solver)

The solution is obtained using the layered model and demonstrates its
robustness and a degree of realism even for this complex case. Note
also the interesting longitudinal "scars". */

#include "grid/quadtree.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "input.h"

int read_xy_float(char * fname, scalar s, int dlev){
  unsigned long long int size = (1 << (dimension*dlev));
  /**
  In serial, life is a bit easier.
  */
  float * a = (float*) malloc (sizeof(float)*size);
  FILE * fp = fopen (fname, "rb");
  fread (a, sizeof(float), size, fp);
  int nxoffset = 0;
  int nyoffset = 0;
  unsigned long long int CR = (1 << dlev);
  // o is found to be -3 which is causing the trouble in Antoon's code
  // int o = -BGHOSTS - 1; 
  int o = -BGHOSTS;
  unsigned long long int index;
  /** Loading the data itself is now straightforward*/
  // int iternumber = 0;
  foreach(){
    index = ((nxoffset + point.i + o) + (CR*(nyoffset + point.j + o)));
    // fprintf(ferr, "point.i = %d, point.j = %d\n", point.i, point.j);
    s[] = (double)a[index];
    // iternumber += 1;
  }
  // fprintf(ferr, "iteration number = %d", iternumber);  
  return 1;
}

/**
The initial condition is a externally imported wave field. Some controlling parameters. */

#define g_ 9.8
double ETAE = 0.1; // refinement criteria for eta
int MAXLEVEL = 8; // max level of refinement in adapt_wavelet function
int MINLEVEL = 6; // min level
double TEND = 50.; // t end
int NLAYER = 10; // number of layers

/** 
		The adaptive function. */
/* int my_adapt() { */
/* #if QUADTREE */
/*   /\* scalar eta[]; *\/ */
/*   /\* foreach() *\/ */
/*   /\*   eta[] = h[] > dry ? h[] + zb[] : 0; *\/ */
/*   /\* boundary ({eta}); *\/ */
/* 	fprintf (stderr, "Line1!\n"); */
/* 	fflush (stderr); */
/*   scalar eta[]; */
/*   foreach() { */
/*     eta[] = zb[]; */
/*     for (scalar h in hl) */
/*       eta[] += h[]; */
/*   } */
/* 	astats s = adapt_wavelet ({eta}, (double[]){ETAE}, */
/* 			    MAXLEVEL, MINLEVEL); */
/*   fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/*   return s.nf; */
/* #else // Cartesian */
/*   return 0; */
/* #endif */
/* } */

/**
The domain is periodic in $x$ and resolved using 256$^2$
points and 60 layers. */

int main(int argc, char * argv[])
{
	if (argc > 1)
		NLAYER = atoi(argv[1]);
  if (argc > 2)
    MAXLEVEL = atoi(argv[2]);
  if (argc > 3)
    MINLEVEL = atoi(argv[3]);
  if (argc > 4)
    ETAE = atof(argv[4]);
  if (argc > 5)
    TEND = atof(argv[5]);
	L0 = 50.;
  origin (-L0/2., -L0/2.);
  periodic (right);
	periodic (top);
  N = 128; // start with a grid of 128
  nl = 10;
  G = g_;
  nu = 1/1000000.;
  run();
}

/**
The intial conditions for the free-surface and velocity are given by
the input field. */

event init (i = 0)
{
	/** Import surface position. h of each layer is (H-zb)/nl. */
  /* FILE * feta = fopen ("./pre/eta", "r"); */
	scalar surface[]; // H is not defined previously
	// int LEVEL_data = sqrt(N);
  read_xy_float ("./pre/eta", surface, 7);
  /* fclose (feta); */
	fprintf(stderr, "Read in eta!\n");
	foreach() {
		zb[] = -5; // zb[] is a scalar field already defined in header file
		eta[] = surface[];
		scalar h,w;
		for (h,w in hl,wl) {
			h[] = (eta[]-zb[])/nl;
			w[] = 0;
		}
	}
	/** Import velocity. */
	char filename_u[50], filename_v[50], filename_w[50];
	/* for (i=0; i<nl; i++) { */
	/* 	for (u,w in ul,wl) { */
	/* 		sprintf (filename_u, "./pre/u_layer%d", i); */
	/* 		sprintf (filename_v, "./pre/v_layer%d", i); */
	/* 		sprintf (filename_w, "./pre/w_layer%d", i); */
	/* 		read_xy_float (filename_u, u.x, 7); */
	/* 		fprintf(stderr, "Read in u, i = %d!\n", i); */
	/* 		read_xy_float (filename_v, u.y, 7); */
	/* 		read_xy_float (filename_w, w, 7); */
	/* 	} */
	/* } */
	for (int ii=0; ii<nl; ii++) {
			sprintf (filename_u, "./pre/u_layer%d", ii);
			sprintf (filename_v, "./pre/v_layer%d", ii);
			sprintf (filename_w, "./pre/w_layer%d", ii);
			read_xy_float (filename_u, ul[ii].x, 7);
			read_xy_float (filename_v, ul[ii].y, 7);
			read_xy_float (filename_w, wl[ii], 7);
			fprintf(stderr, "Read in velocity, layer index = %d!\n", ii);
	}

	/* // test limiter */
	/* scalar u_surface = ul[nl-1].x; */
	/* foreach() { */
	/* 	if (u_surface[] > 50) */
	/* 		u_surface[] = 50; */
	/* } */
	/* ul[nl-1].x = u_surface; */
}

/** 
We log the evolution of kinetic and potential energy.
*/

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    scalar h, w;
    vector u;
    double zc = zb[];
    for (h,w,u in hl,wl,ul) {
      double norm2 = sq(w[]);
      foreach_dimension()
			norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
	static FILE * fp = fopen("energy.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe);
  fflush (fp);
}

/**
Note that the movie generation below is very expensive. */

#if 1
event movie (t += 1; t <= 50)
{
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
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f'\n"
	   "splot '-' u 1:2:3:4 w pm3d\n",
	   i/3, t);
  scalar H[];
  foreach() {
    H[] = zb[];
    for (scalar h in hl)
      H[] += h[];
  }
  boundary ({H});
  scalar ux = ul[nl-1].x;
  output_field ({H,ux}, fp, linear = true);
  fprintf (fp, "e\n\n");
  fflush (fp);
	char filename1[50], filename2[50];
	sprintf (filename1, "./surface/eta_matrix_%g", t);
	sprintf (filename2, "./surface/eta_%g", t);
	FILE * feta1 = fopen (filename1, "w");
	// Might need to change to mpi function later
	output_matrix (H, feta1, N, linear = true);
	fclose (feta1);
	FILE * feta2 = fopen (filename2, "w");
	output_field ({H}, feta2, linear = true);
	fclose (feta2);
}

/* event moviemaker (t = end) */
/* { */
/*   system ("for f in plot*.png; do convert $f ppm:-; done | ppm2mp4 movie.mp4"); */
/* } */
#endif

/* event adapt (i++) { */
/* 	fprintf(stderr, "Adapting start!\n"); */
/* 	fflush(stderr); */
/* 	my_adapt(); */
/* } */

event endrun (t = 10) {
	dump();
}

