/**
	 # breaking Stokes wave

	 A steep, third-order Stokes wave is unstable and breaks.

	 ![Animation of the free-surface](stokes/movie.mp4)

	 The solution obtained using the layered model matches the
	 Navier-Stokes/VOF solution remarkably well, even after breaking.

	 ~~~gnuplot Wave evolution: layered (left column) and Navier-Stokes/VOF (right column) { width=100% }
	 unset key
	 unset xtics
	 unset ytics
	 unset border
	 set multiplot layout 1,2
	 set size ratio -1
	 plot for [i = 0:10] 'log' index i u 1:($2-0.15*i) w l lc -1 lt 1
	 plot for [i = 0:10] '../stokes-ns/log' index i u 1:($2-0.15*i) w l lc -1 lt 1
	 unset multiplot
	 ~~~

	 See [Popinet (2020)](/Bibliography#popinet2020) for a more detailed
	 discussion and [stokes-ns.c]() for the Navier-Stokes/VOF code. */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
//#include "layered/remap.h"
#include "remap_test.h"
#include "layered/perfs.h"

double TEND=100;
double ak = 0.2;
double RE = 40000.;
double h_ = 10;
double gpe_base = 0; 
#define g_   9.8

double kp_ = 2.*pi/10.;
double P_ = 0.01;
int N_power_ = 5;
#define N_mode_ 16
double F_kx_[N_mode_], omega[N_mode_], phase[N_mode_];
double kx_[N_mode_];
double dkx_;
// Breaking time and place
double xb_ = 0, tb_ = 40;

void power_input1D () {
	for (int i=0; i<N_mode_; i++)
		kx_[i] = 2.*pi/L0*(i+1);
	dkx_ = kx_[1]-kx_[0];
	// A function that reads in F_kx. Next step is to generate F_kxky_ all inside
	double size = N_mode_;
	for (int i=0;i<size;i++) {
		F_kx_[i] = P_*pow(kx_[i],-2.5)*exp(-0.75*sq(kp_/kx_[i]));
		fprintf(stderr, "%g %g \n", kx_[i], F_kx_[i]);
	}
	// Phase and omega, next focusing phase
	double kmod = 0;
	for (int i=0; i<N_mode_;i++) {
		kmod = kx_[i];
		omega[i] = sqrt(g_*kmod);
		phase[i] = -kx_[i]*xb_ + omega[i]*tb_;
	}
}

int main(int argc, char * argv[])
{
	if (argc > 1)
		RE = atof (argv[1]);
	if (argc > 2)
		h_ = atof (argv[2]);
	if (argc > 3)
		nl = atoi (argv[3]);
	if (argc > 4)
		coeff = atof (argv[4]);
	else
		nl = 60;
	if (argc > 5)
		N = atoi (argv[5]);
	else
		N = 256;
	L0 = 50.;
  origin (-L0/2.);
  periodic (right);
  G = g_;
  nu = 1/RE;
	gpe_base = -0.5*sq(h_)*L0*g_;
  // max_slope = 1.; // a bit less dissipative
  run();
}

//#include "test/stokes.h"

double wave1D (double x, double y)
{
	double eta = 0;
	double ampl = 0, a = 0;
	for (int i=0; i<N_mode_;i++) {
			ampl = sqrt(2.*F_kx_[i]*dkx_);
			a = (kx_[i]*x + phase[i]);
			eta += ampl*cos(a);
	}
	return eta;
}
double u_x1D (double x, double y)
{
	double u_x = 0;
	double ampl = 0, a = 0;
	double z_actual = 0, kmod = 0;
	for (int i=0; i<N_mode_;i++) {
			ampl = sqrt(2.*F_kx_[i]*dkx_);
			z_actual = (y < ampl ? (y) : ampl);
			kmod = kx_[i];
			a = kx_[i]*x + phase[i];
			u_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a);
	}
	return u_x;
}
double u_y1D (double x, double y)
{
	double u_y = 0;
	double ampl = 0, a = 0;
	double z_actual = 0, kmod = 0;
	for (int i=0; i<N_mode_;i++) {
			ampl = sqrt(2.*F_kx_[i]*dkx_);
			z_actual = (y < ampl ? (y) : ampl);
			kmod = kx_[i];
			a = kx_[i]*x + phase[i];
			u_y += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
	}
	return u_y;
}

event init (i = 0)
{
	power_input1D();
	foreach() {
    zb[] = -h_;
		eta[] = wave1D(x, 0);
    double H = wave1D(x, 0) - zb[];
    foreach_layer() {
      h[] = H/nl;
    }
	}
	// remap?
	vertical_remapping (h, tracers);
	foreach() {
		double z = zb[];
		foreach_layer() {
			z += h[]/2.;
      u.x[] = u_x1D(x, z);
      w[] = u_y1D(x, z);
      z += h[]/2.;
		}
  }
}

/* event profiles (t += T0/4.; t <= 2.5*T0) { */
/*   foreach_leaf() { */
/*     double H = zb[]; */
/*     foreach_layer() */
/*       H += h[]; */
/*     fprintf (stderr, "%g %g\n", x, H); */
/*   } */
/*   fprintf (stderr, "\n\n"); */
/* } */

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
				norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
	static FILE * fp = fopen ("budget.dat", "w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe-gpe_base);
	fflush (fp);
}

event movie (t += 1; t < TEND)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
						 "set size ratio -1\n");  
  fprintf (fp,
					 "set output 'plot%04d.png'\n"
					 "set title 't = %.2f'\n"
					 "p [%g:%g][-5:5]'-' u 1:3:2 w filledcu lc 3 t ''",
					 i/3, t, X0, X0 + L0);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}

/* event moviemaker (t = end) */
/* { */
/*   system ("rm -f movie.mp4 && " */
/* 					"ffmpeg -r 25 -f image2 -i plot%04d.png " */
/* 					"-vcodec libx264 -vf format=yuv420p -movflags +faststart " */
/* 					"movie.mp4 2> /dev/null && " */
/* 					"rm -f plot*.png"); */
/* } */

event dumpstart (i = 0) {
	static FILE * fp = fopen ("grid_start.dat", "w");
  fprintf (fp, "x h ux\n");
	foreach_leaf () {
		foreach_layer()
			fprintf(fp, "%g %g %g\n", x, h[], u.x[]);
	}
	fflush (fp);
}

event dumpend (t=end) {
	static FILE * fp = fopen ("grid_end.dat", "w");
  fprintf (fp, "x h ux\n");
	foreach_leaf () {
		foreach_layer()
			fprintf(fp, "%g %g %g\n", x, h[], u.x[]);
	}
	fflush (fp);
}
