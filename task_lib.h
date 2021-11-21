#ifndef TASK_H
#define TASK_H

//==============================================================================================================
//    GLOBAL CONSTANTS
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    GAME COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define N_BALLS		16				// number of balls
#define N_THREADS	(N_BALLS + 2)			// number of threads
#define N_POCKETS	6				// number of pockets
#define N_COVERS	6				// number of covers
#define N_BUTTN		5				// number of buttons

//--------------------------------------------------------------------------------------------------------------
//    CONVERSION COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define CM_PIX		4.0				// conversion factor	(centimeter --> pixel)
#define PIX_CM		(1 / CM_PIX)			// conversion factor	(pixel --> centimeter)

//--------------------------------------------------------------------------------------------------------------
//    GEOMETRIC COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define WR		10.0				// width of the wooden rails + cushions	[cm]
#define WC		 2.0				// width of the cushions only		[cm]
#define PK		12.5				// dimension of the pocket		[cm]
#define DB		6.25				// diameter of the ball			[cm]

#define HD		200.0				// horizontal dimension of the table	[cm]
#define VD		100.0				// vertical dimension of the table	[cm]
#define HP		((PIX_CM * HR - HD)/2)		// horizontal position of the table	[cm]
#define VP		((PIX_CM * (VR + WP) - VD)/2)	// vertical position of the table	[cm]

//--------------------------------------------------------------------------------------------------------------
//    DYNAMIC COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define DMP		0.5	// default walls damping coeff		(DMP = 0 --> no damp)
#define DMP_MIN		0.0	// damp minimum value
#define DMP_MAX		1.0	// damp maximum value
#define DMP_RES		0.05	// damp resolution			(variable min increment)

#define ACC		-20.0	// default acceleration	[cm/s^2]	(due to rolling friction)
#define ACC_MIN		-100.0	// accel minimum value	[cm/s^2]
#define ACC_MAX		0.0	// accel maximum value	[cm/s^2]
#define ACC_RES		5.0	// accel resolution	[cm/s^2]	(variable min increment)

#define V_SUP		1200.0	// speed supremum	[cm/s]

//--------------------------------------------------------------------------------------------------------------
//    FLAGS COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define N_FLAGS		6	// number of global flag variables

#define I_CLS		0	// index of the array "flag" associated with the button "close"
#define I_RST		1	// index of the array "flag" associated with the button "reset"
#define I_SGT		2	// index of the array "flag" associated with the button "sight"
#define I_AW1		3	// index of the array "flag" associated with the button "arrw1"
#define I_AW2		4	// index of the array "flag" associated with the button "arrw2"
#define I_PLC		5	// index of the array "flag" associated with the function "place_cue_ball" 

//--------------------------------------------------------------------------------------------------------------
//    THREAD COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define B_PER		7	// ball task period		[ms]
#define U_PER		14	// user task period		[ms]
#define G_PER		28	// graphic task period		[ms]

#define B_DL		7	// ball task deadline		[ms]
#define U_DL		14	// user task deadline		[ms]
#define G_DL		28	// graphic task deadline	[ms]

#define B_PRIO		90	// ball task priority
#define U_PRIO		30	// user task priority
#define G_PRIO		10	// graphic task priority

//--------------------------------------------------------------------------------------------------------------
//    AUXILIARY COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define _GNU_SOURCE

#define PATH_LEN	25	// max length of bitmap file names (path)
#define STRG_LEN	35	// max length of a string
#define TRL_LEN		1000	// length of the array "trail"
#define STEP		2	// step of sampling of the trail

//==============================================================================================================
//    STRUCTURE TYPE DECLARATIONS
//==============================================================================================================

struct task_par {

	int	arg; 			// task argument 
	long 	wcet; 			// worst-case execution time	[ms] 
	int 	period; 		// period			[ms] 
	int 	deadline; 		// relative deadline		[ms] 
	int 	priority; 		// priority level		[0, 99] 
	int 	dmiss; 			// number of misses
	struct 	timespec at;		// next activation time		[ms]
	struct 	timespec dl;		// absolute deadline		[ms]
};

struct state {

	float	px;			// x-component of position vector
	float	py;			// y-component of position vector
	float	vx;			// x-component of velocity vector
	float	vy;			// y-component of velocity vector
};

//==============================================================================================================
//    EXTERN GLOBAL VARIABLE DECLARATIONS
//==============================================================================================================

extern struct state		ball[N_BALLS];		// ball state: position and velocity vectors
extern struct state		trail[TRL_LEN];		// array of trail points

extern int			dl_miss;		// current number of deadline missed
extern float			dmp;			// current damping factor value
extern float			acc;			// current acceleration value

extern pthread_mutex_t		ball_mutex[N_BALLS];	// ball mutex
extern pthread_mutex_t		flag_mutex[N_FLAGS];	// flag mutex
extern pthread_mutex_t		trail_mutex;		// trail mutex
extern pthread_mutex_t		dmp_mutex;		// dmp mutex
extern pthread_mutex_t		acc_mutex;		// acc mutex
extern pthread_mutex_t		dl_mutex;		// dl_miss mutex

extern int			flag[N_FLAGS];		// flag variables

			// flag[0]: state of the close button (1 --> pressed, 0 --> not pressed)
			// flag[1]: state of the reset button (1 --> pressed, 0 --> not pressed)
			// flag[2]: state of the sight button (1 --> pressed, 0 --> not pressed)
			// flag[3]: state of the arrw1 button (1 or 2 --> pressed, 0 --> not pressed)
			// flag[4]: state of the arrw2 button (1 or 2 --> pressed, 0 --> not pressed)

//==============================================================================================================
//    FUNCTION PROTOTYPES
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    TIME MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void	time_copy(struct timespec* td, struct timespec ts);
void	time_add_ms(struct timespec* t, int ms);
int 	time_cmp(struct timespec t1, struct timespec t2); 

//--------------------------------------------------------------------------------------------------------------
//    THREADS MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void	init_threads(void);
void	create_thread(void* task_name, int period, int deadline, int priority, int i);
void	wait_for_termination(void);
void	set_activation(int i);
void	wait_for_activation(int i);
int 	deadline_miss(int i);
int	get_task_index(void* arg);

//--------------------------------------------------------------------------------------------------------------
//    MUTEX MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void		init_mutex(void);

void		set_ball_state(int i, struct state value);
void		set_trail(struct state* value);
void		set_flag_variable(int i, int value);
void		set_acc_variable(float value);
void		set_dmp_variable(float value);

struct state	read_ball_state(int i);
struct state*	read_trail(void);
int		read_flag_variable(int i);
float		read_acc_variable(void);
float		read_dmp_variable(void);
float		read_dl_miss_variable(void);

void		switch_flag_variable(int i);
void		increase_dl_miss_variable(void);

//--------------------------------------------------------------------------------------------------------------
//    BALL MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void		init_balls(void);
void*		ball_task(void* arg);
struct state	inter_collisions(struct state b, float acc, float delta_t);
struct state	wall_collision(struct state b, float dmp);
struct state	ball_collision(struct state b, int j);
struct state	in_pocket(struct state b, int i);

//--------------------------------------------------------------------------------------------------------------
//    GRAPHICS MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void*	graphic_task(void* arg);

//--------------------------------------------------------------------------------------------------------------
//    USER MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void*	user_task(void* arg);
void	turn_off(void);
void	reset(void);
void	modify_acc(void);
void	modify_dmp(void);
int	aim(void);
void	hit_cue_ball(struct state* b);
void	place_cue_ball(struct state* b);
void	compute_trail(struct state* b);

//--------------------------------------------------------------------------------------------------------------
//    AUXILIARY FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void	init(void);
void	init_flags(void);
float	min(float* vect, int dim);
float	norm(float x, float y);
float	dot(float x1, float y1, float x2, float y2);
int	triang_root(int n);
void	generate_triang_config(float* x, float* y, int dim);
void	change_balls_sequence(int* seq, int dim);
void	exchange(int* vect, int i, int j);

//--------------------------------------------------------------------------------------------------------------
//    TEST FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

int	steady_state_test(struct state* b);
int	final_state_test(struct state* b);
int	ball_in_play_test(struct state b);
int	within_table_test(float x, float y, float offset);
int	within_region_test(int x, int y, int a, int b, int delta_a, int delta_b);

//--------------------------------------------------------------------------------------------------------------

#endif
