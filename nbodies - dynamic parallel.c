#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#include <X11/Xlib.h> // X11 library headers
#include <X11/Xutil.h>
#include <X11/Xos.h>
#define X_RESN 1000 /* x resolution */
#define Y_RESN 1000 /* y resolution */
#define UNIVERSE_SIZE 1000 // Universe width/height in km.
#define G 0.0000000000667 // Gravitational Constant
// #define G 3
#define PI  3.14159265358979323846 // Pi

// Body type.
struct Body {
	long m;
	double fx;
	double fy;
	double vx;
	double vy;
	double posx;
	double posy;
	int depth;
};
typedef struct Body Body;

// function prototype
Display * x11setup(Window *win, GC *gc, int width, int height);
MPI_Datatype createMPIBodyType();
void createBodies(Body *bodies, int n, XColor *colours);
long genRand(long max, long min);
void printBody(Body *a);
void initBody(Body *a, long m, double vx, double vy, double posx, double posy, int depth);
void calcForce(Body *a, Body *b);
void updatePosition(Body *a, double dt);
void setupColours(Display *display, XColor *colours, Colormap screenColourmap);
XColor pickColour(int depth);

int main(int argc, char *argv[])
{
	int rank, nprocs, i, nbodies, x, y;
	Window win; // initialization for a window
	GC gc; // graphics context
	Display *display = NULL;
	unsigned int width = X_RESN, height = Y_RESN; /* window size */
	clock_t start, end, elapsed;
	Colormap screenColourmap;
	XColor colours[14];
	int n = 100;
	Body bodies[n];
	
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Datatype MPI_BODY = createMPIBodyType();
        
	if(rank==0)
	{
		display = x11setup(&win, &gc, width, height);
		screenColourmap = DefaultColormap(display, DefaultScreen(display));
		setupColours(display, colours, screenColourmap);
		time_t t;
		srand((unsigned) time(&t));
		createBodies(bodies, n, colours);
	}
	
	// main loop
	int running = 1; // loop variable
	start = clock();
	int counter = 0;
	while(running) {
		// checks to see if there have been any events,
		// will exit the main loop if any key is pressed
		if(rank==0) {
			if(XPending(display)) {
				printf("Dunno what this does.\n");
				XEvent ev;
				XNextEvent(display, &ev);
				switch(ev.type) {
					case KeyPress:
					printf("Not running\n");
						running = 0;
						break;
				}
			}
		}
		
		

		end = clock();
		elapsed = end - start;
		// Only update the display if > 1 millisecond has passed since the last update.
		if(elapsed / (CLOCKS_PER_SEC/1000) > 1 && rank==0) {
			// printf("Run Count: %d\n", counter++);
 			int probeFlag = 0;
			MPI_Status status;
			int count = 0;
			int tracker[nprocs];
			for (int i = 1; i < nprocs; i++) {
				MPI_Send(bodies, n, MPI_BODY, i, 98, MPI_COMM_WORLD);
			}
			for (int i = 0; i < nprocs-1; i++) {
				MPI_Probe(i+1, 98, MPI_COMM_WORLD, &status);
				int slave = status.MPI_SOURCE;
				MPI_Recv(&probeFlag, 1, MPI_INT, MPI_ANY_SOURCE, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&i, 1, MPI_INT, slave, 99, MPI_COMM_WORLD);
				tracker[slave] = i; 
				// printf("Next: %d, to be processed by: %d\n", tracker[slave], slave);
			}
			for (int i = nprocs-1; i < n; i++) {
				while (!probeFlag) {
					MPI_Iprobe(MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &probeFlag, &status);
					if (probeFlag) {
						int slave = status.MPI_SOURCE;
						MPI_Recv(&bodies[tracker[slave]], 1, MPI_BODY, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Send(&i, 1, MPI_INT, slave, 99, MPI_COMM_WORLD);
						// printf("Body: %d processed by: %d, Next: %d\n", tracker[slave], slave, i);  
						tracker[slave] = i;
						count++;
					}
				}
				probeFlag = 0;
			}
			
			probeFlag = 0;
			for (int i = 1; i < nprocs; i++) {
				MPI_Probe(MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
				int slave = status.MPI_SOURCE;
				// printf("Final result from: %d\n", slave);
				MPI_Recv(&bodies[tracker[slave]], 1, MPI_BODY, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&probeFlag, 1, MPI_INT, i, 98, MPI_COMM_WORLD);
			}
			for (int i = 1; i < nprocs; i++) {
				MPI_Recv(&probeFlag, 1, MPI_INT, i, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("Resignation received from %d\n", i);
			}
			XClearWindow(display, win);
			for (int i = 0; i < n; i++) {
				XColor colour = colours[bodies[i].depth];
				XSetForeground(display, gc, colour.pixel);
				// Draw the bodies. Uses the density of Andesite magma to estimate a size from the mass.
				double diameter = cbrt(((double)bodies[i].m / 2500) / 1000) * 20;
				int posX = (((bodies[i].posx - diameter / 2) / UNIVERSE_SIZE) * X_RESN);
				int posY = (((bodies[i].posy - diameter / 2) / UNIVERSE_SIZE) * X_RESN); // Intentionally X_RESN to force position scaling to be based on screen height.
				XFillArc(display, win, gc, posX, posY, diameter, diameter, 0, 23040);
			}
			start = end;
			XFlush(display);
		} else if (rank != 0) {
			int flag = 1;
			int probeFlag = 0;
			MPI_Status status;
			Body body;
			MPI_Recv(&bodies, n, MPI_BODY, 0, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&flag, 1, MPI_INT, 0, 98, MPI_COMM_WORLD);
			while (flag) {
				MPI_Iprobe(0, 98, MPI_COMM_WORLD, &probeFlag, MPI_STATUS_IGNORE);
				if (probeFlag) {
					MPI_Recv(&flag, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send(&flag, 1, MPI_INT, 0, 98, MPI_COMM_WORLD);
				} else {
					MPI_Iprobe(0, 99, MPI_COMM_WORLD, &probeFlag, MPI_STATUS_IGNORE);
					if (probeFlag) {
						MPI_Recv(&i, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						bodies[i].fx = 0;
						bodies[i].fy = 0;
						for (int j = 0; j < n; j++) {
							if (j != i) {
								calcForce(&bodies[i], &bodies[j]);
							}
						}
						updatePosition(&bodies[i], 200);
						MPI_Send(&bodies[i], 1, MPI_BODY, 0, 99, MPI_COMM_WORLD);
					}
				}
			}
		}
	}

	if(rank==0 && display) {
		XCloseDisplay(display); // close the display window
	}
	
	MPI_Finalize();
	return 0;
}

MPI_Datatype createMPIBodyType() {
	int count = 8;
	int array_of_blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Aint array_of_displacements[8];
	MPI_Datatype array_of_types[8] = {MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Datatype MPI_BODY;
	MPI_Aint longextent, doubleextent, intextent;
	MPI_Type_extent(MPI_LONG, &longextent);
	MPI_Type_extent(MPI_INT, &intextent);
	MPI_Type_extent(MPI_DOUBLE, &doubleextent);
	array_of_displacements[0] = 0;
	array_of_displacements[1] = longextent;
	array_of_displacements[2] = array_of_displacements[1] + doubleextent;
	array_of_displacements[3] = array_of_displacements[2] + doubleextent;
	array_of_displacements[4] = array_of_displacements[3] + doubleextent;
	array_of_displacements[5] = array_of_displacements[4] + doubleextent;
	array_of_displacements[6] = array_of_displacements[5] + doubleextent;
	array_of_displacements[7] = array_of_displacements[6] + doubleextent;
	MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_BODY);
	MPI_Type_commit(&MPI_BODY);
	return MPI_BODY;
}

// Creates n bodies and stores them into the bodies array.
void createBodies(Body *bodies, int n, XColor *colours) {
	for (int i = 0; i < n; i++) {
		double posx = genRand(1000, 0);
		double posy = genRand(1000, 0);
		double vx = genRand(0, 0);
		double vy = genRand(0, 0);
		long m = genRand(5000000, 1000000);
		int depth = (int)genRand(13, 0);
		initBody(&bodies[i], m, vx, vy, posx, posy, depth);
	}
}

// Generates a random long between max and min.
long genRand(long max, long min) {
	double randBase = (double)rand() / (double)RAND_MAX; //Generate random between 0 and 1.
	return (randBase * (max - min)) + min;
}

// Prints out all variables within a Body.
void printBody(Body *a) {
	printf("Position = %f, %f. Force = %f, %f. Velocity = %f, %f. Mass = %d. Colour = %d\n", (*a).posx, (*a).posy, (*a).fx, (*a).fy, (*a).vx, (*a).vy, (*a).m, (*a).depth);
}

// Initialise values within the body.
void initBody(Body *a, long m, double vx, double vy, double posx, double posy, int depth) {
	(*a).m = m;
	(*a).fx = 0;
	(*a).fy = 0;
	(*a).vx = vx;
	(*a).vy = vy;
	(*a).posx = posx;
	(*a).posy = posy;
	(*a).depth = depth;
}

// Calculate the forces on the body and apply to Body a.
void calcForce(Body *a, Body *b) {
	double dx = (*b).posx - (*a).posx;
	double dy = (*b).posy - (*a).posy;
	double r = sqrt(pow(dx, 2) + pow(dy, 2));
	if (r != 0.0) {
		double F = (G * (*a).m * (*b).m) / pow(r, 2);
		double rmod = sqrt((*a).m) / 100 + sqrt((*b).m) / 100; // Modifier to dampen the effect of r values close to 0.
		(*a).fx += F * dx / (r + rmod);
		(*a).fy += F * dy / (r + rmod);
	}
}

// Update the velocity and X/Y positions of Body a with respect to the time interval dt.
void updatePosition(Body *a, double dt) {
	(*a).vx += dt * (*a).fx / (*a).m;
	(*a).vy += dt * (*a).fy / (*a).m;
	(*a).posx += dt * (*a).vx;
	(*a).posy += dt * (*a).vy;
}


// Allocates a set of colours for use by way of a Colourmap.
void setupColours(Display *display, XColor *colours, Colormap screenColourmap) {
	for (int i = 0; i < 14; i++) {
		colours[i] = pickColour(i);
		XAllocColor(display, screenColourmap, &colours[i]);
	}
}

// Picks a colour based on an integer input between 0 and 13.
XColor pickColour(int depth) {
	double tau = PI * 2;
	double tauRainbow = tau / 14;
	double tauPart = tau / 3;
	XColor colour;
	colour.red = (int)(sin(tauRainbow * depth + 1 * tauPart) * 32767 + 32768);
	colour.green = (int)(sin(tauRainbow * depth + 2 * tauPart) * 32767 + 32768);
	colour.blue = (int)(sin(tauRainbow * depth + 3 * tauPart) * 32767 + 32768);
	return colour;
}

Display * x11setup(Window *win, GC *gc, int width, int height)
{
	
	/* --------------------------- X11 graphics setup ------------------------------ */
	Display 		*display;
	unsigned int 	win_x,win_y, /* window position */
					border_width, /* border width in pixels */
					display_width, display_height, /* size of screen */
					screen; /* which screen */
	
	char 			window_name[] = "N-Body Simulation", *display_name = NULL;
	unsigned long 	valuemask = 0;
	XGCValues 		values;
	
	XSizeHints 		size_hints;
	
	//Pixmap 		bitmap;
	//XPoint 		points[800];
	FILE 			*fopen ();//, *fp;
	//char 			str[100];
	
	XSetWindowAttributes attr[1];
	
	if ( (display = XOpenDisplay (display_name)) == NULL ) { /* connect to Xserver */
		fprintf (stderr, "Cannot connect to X server %s\n",XDisplayName (display_name) );
		exit (-1);
	}
	
	screen = DefaultScreen (display); /* get screen size */
	display_width = DisplayWidth (display, screen);
	display_height = DisplayHeight (display, screen);
	
	win_x = 0; win_y = 0; /* set window position */
	
	border_width = 4; /* create opaque window */
	*win = XCreateSimpleWindow (display, RootWindow (display, screen),
			win_x, win_y, width, height, border_width,
			WhitePixel (display, screen), BlackPixel (display, screen));
			
	size_hints.flags = USPosition|USSize;
	size_hints.x = win_x;
	size_hints.y = win_y;
	size_hints.width = width;
	size_hints.height = height;
	size_hints.min_width = 300;
	size_hints.min_height = 300;
	
	XSetNormalHints (display, *win, &size_hints);
	XStoreName(display, *win, window_name);
	
	*gc = XCreateGC (display, *win, valuemask, &values); /* create graphics context */
	
	XSetBackground (display, *gc, BlackPixel (display, screen));
	XSetForeground (display, *gc, WhitePixel (display, screen));
	XSetLineAttributes (display, *gc, 1, LineSolid, CapRound, JoinRound);
	
	attr[0].backing_store = Always;
	attr[0].backing_planes = 1;
	attr[0].backing_pixel = BlackPixel(display, screen);
	
	XChangeWindowAttributes(display, *win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);
	
	XSelectInput(display, *win, KeyPressMask);
	
	XMapWindow (display, *win);
	XSync(display, 0);
	
	/* --------------------------- End of X11 graphics setup ------------------------------ */
	return display;
}