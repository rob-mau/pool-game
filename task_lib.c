#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sched.h>
#include <allegro.h>
#include "task_lib.h"
#include "graphics_lib.h"

//==============================================================================================================
//    GLOBAL VARIABLE DEFINITIONS
//==============================================================================================================

struct sched_param	mypar;				// scheduling parameters
struct task_par		tp[N_THREADS];		// task parameters
struct state		ball[N_BALLS];		// ball state: position and velocity vectors
struct state		trail[TRL_LEN];		// array of trail points

int			dl_miss;					// current number of missed deadlines
float			dmp;					// current damping factor
float			acc;					// current acceleration

pthread_t			tid[N_THREADS];		// thread identifiers
pthread_attr_t		att[N_THREADS];		// thread attributes

pthread_mutex_t		ball_mutex[N_BALLS];	// ball mutex
pthread_mutex_t		flag_mutex[N_FLAGS];	// flag mutex
pthread_mutex_t		trail_mutex;			// trail mutex
pthread_mutex_t		acc_mutex;				// acc mutex
pthread_mutex_t		dmp_mutex;				// dmp mutex
pthread_mutex_t		dl_mutex;				// dl_miss mutex
pthread_mutexattr_t	mutexatt;				// mutex attributes

int			flag[N_FLAGS];		// flag variables

		// flag[0]: state of the close button (1 --> pressed, 0 --> not pressed)
		// flag[1]: state of the reset button (1 --> pressed, 0 --> not pressed)
		// flag[2]: state of the sight button (1 --> pressed, 0 --> not pressed)
		// flag[3]: state of the arrw1 button (1 or 2 --> pressed, 0 --> not pressed)
		// flag[4]: state of the arrw2 button (1 or 2 --> pressed, 0 --> not pressed)

//==============================================================================================================
//    FUNCTION DEFINITIONS
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    TIME MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    time_copy() functon: 
//        copies a source time variable "ts" in a destination variable pointed by "td"
//--------------------------------------------------------------------------------------------------------------

void	time_copy(struct timespec* td, struct timespec ts)
{
	td->tv_sec  = ts.tv_sec;
	td->tv_nsec = ts.tv_nsec;	
}

//-------------------------------------------------------------------------------------------------------------
//    time_add_ms() function: 
//        adds a value "ms", expressed in milliseconds, to the time variable pointed by "t"
//-------------------------------------------------------------------------------------------------------------

void	time_add_ms(struct timespec* t, int ms)
{
	t->tv_sec  +=  ms/1000;
	t->tv_nsec += (ms%1000) * 1000000;  

	if (t->tv_nsec > 1000000000) {		// in case tv_nsec exceed 1 seconds
		t->tv_nsec -= 1000000000;
		t->tv_sec  += 1; 
	}  
}

//--------------------------------------------------------------------------------------------------------------
//    time_cmp() function: 
//        compares two time variables "t1" and "t2" and returns: 0 if t1 = t2, 1 if t1 > t2, -1 if t1 < t2
//--------------------------------------------------------------------------------------------------------------

int	time_cmp(struct timespec t1, struct timespec t2)
{
	if (t1.tv_sec > t2.tv_sec)		return 1;
	if (t1.tv_sec < t2.tv_sec)		return -1;
	if (t1.tv_nsec > t2.tv_nsec)	return 1;
	if (t1.tv_nsec < t2.tv_nsec)	return -1;
	return 0;	
}

//--------------------------------------------------------------------------------------------------------------
//    THREADS MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    init_threads() function: initializes all the threads
//--------------------------------------------------------------------------------------------------------------

void	init_threads(void)
{
int	i;

	for (i = 0; i < N_BALLS; i++) create_thread(ball_task, B_PER, B_DL, B_PRIO, i);
	create_thread(graphic_task, G_PER, G_DL, G_PRIO, N_BALLS);
	create_thread(user_task, U_PER, U_DL, U_PRIO, N_BALLS + 1);
}

//--------------------------------------------------------------------------------------------------------------
//    create_thread() function: creates a thread managed by SCHED_RR policy
//--------------------------------------------------------------------------------------------------------------

void	create_thread(void* task_function, int period, int deadline, int priority, int i)
{
	tp[i].arg = i;
	tp[i].period = period;
	tp[i].deadline = deadline;
	tp[i].priority = priority;
	
	pthread_attr_init(&att[i]);
	pthread_attr_setinheritsched(&att[i], PTHREAD_EXPLICIT_SCHED);
	pthread_attr_setschedpolicy(&att[i], SCHED_RR);	

	mypar.sched_priority = tp[i].priority;
	pthread_attr_setschedparam(&att[i], &mypar);
	
	pthread_create(&tid[i], &att[i], task_function, (void*)&tp[i]);	
}

//--------------------------------------------------------------------------------------------------------------
//    wait_for_termination() function: wait for thread termination
//--------------------------------------------------------------------------------------------------------------

void	wait_for_termination(void)
{
int	k;

	for (k = 0; k < N_THREADS; k++) pthread_join(tid[k], NULL);
}

//--------------------------------------------------------------------------------------------------------------
//    set_activation() function: 
//        computes the next activation time "at" and the absolute deadline of the task "dl"
//--------------------------------------------------------------------------------------------------------------

void	set_activation(int i)
{
struct	timespec t;

	clock_gettime(CLOCK_MONOTONIC, &t);

	time_copy(&(tp[i].at), t);
	time_copy(&(tp[i].dl), t);

	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].deadline);
}

//--------------------------------------------------------------------------------------------------------------
//    wait_for_activation() function:
//        suspends the calling thread until the next activation and updates the variables "at" and "dl"
//--------------------------------------------------------------------------------------------------------------

void	wait_for_activation(int i)
{
	clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &(tp[i].at), NULL);

	time_add_ms(&(tp[i].at), tp[i].period);		
	time_add_ms(&(tp[i].dl), tp[i].period);
}

//--------------------------------------------------------------------------------------------------------------
//    deadline_miss() function: 
//        if the thread is executing when is reactivated, it increases the "dmiss" and return 1, else return 0
//--------------------------------------------------------------------------------------------------------------

int	deadline_miss(int i)
{
struct	timespec now;

	clock_gettime(CLOCK_MONOTONIC, &now);

	if (time_cmp(now, tp[i].dl) > 0) {
		tp[i].dmiss++;
		return 1;
	}
	return 0;
}

//--------------------------------------------------------------------------------------------------------------
//    get_task_index() function: returns the index of the task
//--------------------------------------------------------------------------------------------------------------

int	get_task_index(void* arg)
{
struct	task_par* tp;
	
	tp = (struct task_par*) arg;
	return tp->arg;
}

//--------------------------------------------------------------------------------------------------------------
//    MUTEX MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    init_mutex() function: initializes mutual exclusions
//--------------------------------------------------------------------------------------------------------------

void	init_mutex(void)
{
int	k;

	pthread_mutexattr_init(&mutexatt);
	pthread_mutexattr_setprotocol(&mutexatt, PTHREAD_PRIO_INHERIT);
	
	for (k = 0; k < N_BALLS; k++) pthread_mutex_init(&ball_mutex[k], &mutexatt);
	for (k = 0; k < N_FLAGS; k++) pthread_mutex_init(&flag_mutex[k], &mutexatt);
	pthread_mutex_init(&acc_mutex, &mutexatt);
	pthread_mutex_init(&dmp_mutex, &mutexatt);
	pthread_mutex_init(&dl_mutex, &mutexatt);
	pthread_mutex_init(&trail_mutex, &mutexatt);
}

//--------------------------------------------------------------------------------------------------------------
//    set_ball_state() function: writes on the global variable "ball[i]"
//--------------------------------------------------------------------------------------------------------------

void	set_ball_state(int i, struct state value)
{
	pthread_mutex_lock(&ball_mutex[i]);
	ball[i] = value;
	pthread_mutex_unlock(&ball_mutex[i]);
}

//--------------------------------------------------------------------------------------------------------------
//    set_trail() function: writes on the global variable "trail"
//--------------------------------------------------------------------------------------------------------------

void	set_trail(struct state* value)
{
int	k;

	pthread_mutex_lock(&trail_mutex);
	for (k = 0; k < TRL_LEN; k++) trail[k] = value[k];
	pthread_mutex_unlock(&trail_mutex);
}

//--------------------------------------------------------------------------------------------------------------
//    set_flag_variable() function: writes on the global variable "flag[i]"
//--------------------------------------------------------------------------------------------------------------

void	set_flag_variable(int i, int value)
{
	pthread_mutex_lock(&flag_mutex[i]);
	flag[i] = value;
	pthread_mutex_unlock(&flag_mutex[i]);
}

//--------------------------------------------------------------------------------------------------------------
//    set_acc_variable() function: writes on the global variable "acc"
//--------------------------------------------------------------------------------------------------------------

void	set_acc_variable(float value)
{
	pthread_mutex_lock(&acc_mutex);
	acc = value;
	pthread_mutex_unlock(&acc_mutex);
}

//--------------------------------------------------------------------------------------------------------------
//    set_dmp_variable() function: writes on the global variable "dmp"
//--------------------------------------------------------------------------------------------------------------

void	set_dmp_variable(float value)
{
	pthread_mutex_lock(&dmp_mutex);
	dmp = value;
	pthread_mutex_unlock(&dmp_mutex);
}

//--------------------------------------------------------------------------------------------------------------
//    read_ball_state() function: reads from the global variable "ball[i]"
//--------------------------------------------------------------------------------------------------------------

struct state	read_ball_state(int i)
{
struct state	value;

	pthread_mutex_lock(&ball_mutex[i]);
	value = ball[i];
	pthread_mutex_unlock(&ball_mutex[i]);

	return value;
}

//--------------------------------------------------------------------------------------------------------------
//    read_trail() function: reads from the global variable "trail"
//--------------------------------------------------------------------------------------------------------------

struct state*	read_trail(void)
{
static struct state	value[TRL_LEN];
int k;

	pthread_mutex_lock(&trail_mutex);
	for (k = 0; k < TRL_LEN; k++) value[k] = trail[k];
	pthread_mutex_unlock(&trail_mutex);

	return value;
}

//--------------------------------------------------------------------------------------------------------------
//    read_flag_variable() function: reads from the global variable "flag[i]"
//--------------------------------------------------------------------------------------------------------------

int	read_flag_variable(int i)
{
int	value;

	pthread_mutex_lock(&flag_mutex[i]);
	value = flag[i];
	pthread_mutex_unlock(&flag_mutex[i]);

	return value;
}

//--------------------------------------------------------------------------------------------------------------
//    read_acc_variable() function: reads from the global variable "acc"
//--------------------------------------------------------------------------------------------------------------

float	read_acc_variable(void)
{
float	value;

	pthread_mutex_lock(&acc_mutex);
	value = acc;
	pthread_mutex_unlock(&acc_mutex);

	return value;
}

//--------------------------------------------------------------------------------------------------------------
//    read_dmp_variable() function: reads from the global variable "dmp"
//--------------------------------------------------------------------------------------------------------------

float	read_dmp_variable(void)
{
float	value;

	pthread_mutex_lock(&dmp_mutex);
	value = dmp;
	pthread_mutex_unlock(&dmp_mutex);

	return value;
}

//--------------------------------------------------------------------------------------------------------------
//    read_dl_miss_variable() function: reads from the global variable "dl_miss"
//--------------------------------------------------------------------------------------------------------------

float	read_dl_miss_variable(void)
{
float	value;

	pthread_mutex_lock(&dl_mutex);
	value = dl_miss;
	pthread_mutex_unlock(&dl_mutex);

	return value;
}
//--------------------------------------------------------------------------------------------------------------
//    switch_flag_variable() function: switches the value of the global variable "flag[i]"
//--------------------------------------------------------------------------------------------------------------

void	switch_flag_variable(int i)
{
	pthread_mutex_lock(&flag_mutex[i]);
	flag[i] = (flag[i] + 1) % 2; 
	pthread_mutex_unlock(&flag_mutex[i]);
}

//--------------------------------------------------------------------------------------------------------------
//    increase_dl_miss() function: increases the value of the global variable "dl_miss"
//--------------------------------------------------------------------------------------------------------------

void	increase_dl_miss_variable(void)
{
	pthread_mutex_lock(&dl_mutex);
	dl_miss++;
	pthread_mutex_unlock(&dl_mutex);
}

//--------------------------------------------------------------------------------------------------------------
//    BALL MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    init_balls() function: assigns initial values to the array of structures "ball"
//--------------------------------------------------------------------------------------------------------------

void	init_balls(void)
{
struct	state start_state[N_BALLS];	// initial states
int	triang_dim;			// dimension of the triangular configuration
float	triang_x[N_BALLS - 1];		// x-component of the position on the triang. config.
float	triang_y[N_BALLS - 1];		// y-component of the position on the triang. config.
int	sequence[N_BALLS - 1];		// sequence according to which balls are placed
int	k;				// loop variable	
	
	triang_dim = triang_root(N_BALLS - 1);	// note: (N_BALLS - 1) must be a triangular number
	generate_triang_config(triang_x, triang_y, triang_dim); // init "triang_x" and "triang_y"

	for (k = 1; k < N_BALLS; k++) sequence[k - 1] = k;	// init "sequence"
	change_balls_sequence(sequence, triang_dim); // rules for 8-ball pool. Don't use if N_BALLS != 16 

	// assignment of initial states to the balls
	start_state[0].px = HP + HD / 4;		
	start_state[0].py = VP + VD / 2;
	start_state[0].vx = 0;
	start_state[0].vy = 0;

	set_ball_state(0, start_state[0]);

	for (k = 1; k < N_BALLS; k++) {		
		start_state[k].px = HP + 3 * HD / 4 + triang_x[k - 1];		
		start_state[k].py = VP + 1 * VD / 2 + triang_y[k - 1];
		start_state[k].vx = 0;	
		start_state[k].vy = 0;

		set_ball_state(sequence[k - 1], start_state[k]);
	}	
}

//--------------------------------------------------------------------------------------------------------------
//    ball_task() function: manages single ball kinematics
//--------------------------------------------------------------------------------------------------------------

void*	ball_task(void* arg)
{
int	i;			// index of the ball using this function (task index)
int	j;			// index of the other balls
float	delta_t;		// sample time in seconds
float	acc_copy;		// copy of the global variable "acc"
float	dmp_copy;		// copy of the global variable "dmp"
struct	state ball_copy;	// copy of the global variable "ball[i]"

	delta_t = (float) B_PER / 1000;

	i = get_task_index(arg);
	set_activation(i);
	while(1) {

		ball_copy = read_ball_state(i);
		if (norm(ball_copy.vx, ball_copy.vy) > 0) {

			acc_copy = read_acc_variable();
			dmp_copy = read_dmp_variable();

			// kinematics calculation
			ball_copy = inter_collisions(ball_copy, acc_copy, delta_t);
			ball_copy = wall_collision(ball_copy, dmp_copy);
			for (j = 0; j < N_BALLS; j++) {
				if (i != j) ball_copy = ball_collision(ball_copy, j);
			}
			ball_copy = in_pocket(ball_copy, i);

			// updating of the state "ball[i]"
			set_ball_state(i, ball_copy);
		}

		if (deadline_miss(i)) {
			increase_dl_miss_variable();
			printf("deadline missed by the ball task (n°%2d)\n", i);
		}
		wait_for_activation(i);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    inter_collisions() function: implements the law of motion of the ball between two consecutive collisions
//--------------------------------------------------------------------------------------------------------------

struct	state inter_collisions(struct state b, float acc, float delta_t)
{
float	v_norm;			// norm of the initial velocity vector
float	ax, ay;			// acceleration	vector
float	delta_vx, delta_vy;	// variation of the velocity after "delta_t" seconds
float	delta_px, delta_py;	// variation of the position after "delta_t" seconds
float	new_vx, new_vy;		// velocity vector after "delta_t" seconds

	v_norm = norm(b.vx, b.vy);

	ax = b.vx / v_norm * acc;
	ay = b.vy / v_norm * acc;

	delta_vx = ax * delta_t;
	delta_vy = ay * delta_t;

	new_vx = b.vx + delta_vx;
	new_vy = b.vy + delta_vy;
	
	if (b.vx * new_vx >= 0) {
		delta_px = (b.vx + new_vx) / 2 * delta_t;
		b.px += delta_px;
		b.vx = new_vx;
	} 
	else 	b.vx = 0;

	if (b.vy * new_vy >= 0) {
		delta_py = (b.vy + new_vy) / 2 * delta_t;
		b.py += delta_py;
		b.vy = new_vy;
	}
	else 	b.vy = 0;

	return b;
}

//--------------------------------------------------------------------------------------------------------------
//    wall_collision() function: implements the collision between the ball and the table cushions
//--------------------------------------------------------------------------------------------------------------

struct	state wall_collision(struct state b, float dmp)
{
int	k;			// loop variable
float	r = DB / 2;		// ball radius
float	dx, dy;			// difference between two position vectors (ball and pocket edge)
float	dist[2 * N_POCKETS];	// distance: norm of the vector [dx, dy]
float	wx, wy;			// b-velocity component along the straight line joining ball and pocket edge
float	w_coeff;		// auxiliary coefficient used to compute wx and wy
float	corr_px, corr_py;	// correction of the position after the collision
float	corr_coeff;		// correction coefficient

// positions of the pocket edges
float	x[2 * N_POCKETS] = {
		
		HP + PK / sqrt(2),
		HP - PK / sqrt(2) + HD,
		HP + (HD - PK) / 2,
		HP + (HD + PK) / 2,

		HP + PK / sqrt(2),
		HP - PK / sqrt(2) + HD,
		HP + (HD - PK) / 2,
		HP + (HD + PK) / 2,

		HP,
	 	HP, 
		HP + HD,
		HP + HD
	};
			
float	y[2 * N_POCKETS] = {
		
		VP,
		VP, 
		VP,
		VP, 

		VP + VD,
		VP + VD,
		VP + VD,
		VP + VD,

		VP + PK / sqrt(2),
		VP - PK / sqrt(2) + VD,
		VP + PK / sqrt(2),
		VP - PK / sqrt(2) + VD
	};

	// walls management
	if (b.px < HP + r) {
		if (b.py > y[10] && b.py < y[11]) {
			b.vx *= (dmp - 1);
			b.px = HP + r;		// correction of the interpenetration effect
						// between ball and table
		}
	}

	if (b.px > HP + HD - r) {
		if (b.py > y[10] && b.py < y[11]) {
			b.vx *= (dmp - 1);
			b.px = HP + HD - r;	// correction of the interpenetration effect
		}
	}

	if (b.py < VP + r) {
		if ((b.px > x[0] && b.px < x[2]) || (b.px > x[3] && b.px < x[1])) {
			b.vy *= (dmp - 1);
			b.py = VP + r;		// correction of the interpenetration effect
		}
	}

	if (b.py > VP + VD - r) {
		if ((b.px > x[0] && b.px < x[2]) || (b.px > x[3] && b.px < x[1])) {
			b.vy *= (dmp - 1);
			b.py = VP + VD - r;	// correction of the interpenetration effect
		}
	}

	// edges management
	for (k = 0; k < (2 * N_POCKETS); k++) {
		dx = b.px - x[k];
		dy = b.py - y[k];
		dist[k] = norm(dx, dy);	

		if (dist[k] > 0 && dist[k] < r) {

			w_coeff = dot(b.vx, b.vy, dx, dy) / (dist[k] * dist[k]);			
			wx = w_coeff * dx;
			wy = w_coeff * dy;

			b.vx += (dmp - 2) * wx;
			b.vy += (dmp - 2) * wy;

			// correction of the interpenetration effect between ball and edge
			corr_coeff = DB / (2 * dist[k]) - 1;	// (if dist = DB/2 --> corr_coeff = 0)
			corr_px = corr_coeff * dx; 
			corr_py = corr_coeff * dy;

			b.px += corr_px;
			b.py += corr_py;
		}
	}
	return b;
}

//--------------------------------------------------------------------------------------------------------------
//    ball_collision() function: implements the collision between two balls ("b1" and "b2")
//--------------------------------------------------------------------------------------------------------------

struct	state ball_collision(struct state b1, int j)
{
struct	state b2;		// copy of the state of the ball that has just been hit 
float 	dx, dy;			// difference between two position vectors (ball b1 and ball b2)
float	dist;			// distance: norm of the vector [dx, dy]
float	wx_1, wy_1;		// b1-velocity component along the straight line joining two ball centers
float	wx_2, wy_2;		// b2-velocity component along the straight line joining two ball centers
float	w_coeff_1;		// auxiliary coefficient used to compute wx_1 and wy_1
float	w_coeff_2;		// auxiliary coefficient used to compute wx_2 and wy_2
float	delta_vx, delta_vy;	// variation of the velocity after the collision
float	corr_coeff;		// correction coefficient 
float	corr_px, corr_py;	// correction of the position after the collision

	b2 = read_ball_state(j);
	dx = b2.px - b1.px;
	dy = b2.py - b1.py;
	dist = norm(dx, dy);

	if (dist > 0 && dist < DB) {
 
		w_coeff_1 = dot(b1.vx, b1.vy, dx, dy) / (dist * dist);
		wx_1 = w_coeff_1 * dx;
		wy_1 = w_coeff_1 * dy; 
	
		w_coeff_2 = dot(b2.vx, b2.vy, dx, dy) / (dist * dist);
		wx_2 = w_coeff_2 * dx;
		wy_2 = w_coeff_2 * dy;
	
		delta_vx = wx_2 - wx_1;
		delta_vy = wy_2 - wy_1;

		b1.vx += delta_vx;
		b1.vy += delta_vy;

		b2.vx -= delta_vx;
		b2.vy -= delta_vy;

		// correction of the interpenetration effect between balls
		corr_coeff = (DB / dist - 1) / 2;	// (if dist = DB --> corr_coeff = 0)
		corr_px = corr_coeff * dx; 
		corr_py = corr_coeff * dy;

		b1.px -= corr_px; 
		b1.py -= corr_py;

		b2.px += corr_px; 
		b2.py += corr_py;

		// updating of the state "ball[j]"
		set_ball_state(j, b2);
	}
	return b1; 
}

//--------------------------------------------------------------------------------------------------------------
//    in_pocket() function: modifies the state of a ball that was just pocketed
//--------------------------------------------------------------------------------------------------------------

struct	state in_pocket(struct state b, int i)
{
float	r = DB / 2;	// ball radius

	if (!within_table_test(b.px, b.py, 5/4 * r)) {

		b.px = (HR - D_SLOT * (N_BALLS - i)) * PIX_CM; 
		b.py = (WP / 2) * PIX_CM;
		b.vx = 0;
		b.vy = 0;
	}
	return b;
}

//--------------------------------------------------------------------------------------------------------------
//    GRAPHICS MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    graphic_task() function: manages the sequence of frames
//--------------------------------------------------------------------------------------------------------------

void*	graphic_task(void* arg)
{
int	i;	// task index

	i = get_task_index(arg);
	set_activation(i);
	while(1) {
		create_frame();

		if (deadline_miss(i)) {
			increase_dl_miss_variable();
			printf("deadline missed by the graphic task\n");
		}
		wait_for_activation(i);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    USER MANAGEMENT FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    user_task() function: manages all user activities
//--------------------------------------------------------------------------------------------------------------

void*	user_task(void* arg)
{
int	i;				// task index
int	k;				// loop variable
int	aim_flag;			// flag variable
struct	state ball_copy[N_BALLS];	// copy of the global variable "ball"

	i = get_task_index(arg);
	set_activation(i);
	while(1) {

		for (k = 0; k < N_BALLS; k++) ball_copy[k] = read_ball_state(k);

		turn_off();
		reset();
		modify_acc();
		modify_dmp();
		aim_flag = aim();

		if (steady_state_test(ball_copy)) {
			if (final_state_test(ball_copy)) init_balls();
			else{
				if (ball_in_play_test(ball_copy[0])) {
					if (aim_flag) compute_trail(ball_copy);
					if (mouse_y > WP) hit_cue_ball(ball_copy);
				}
				else place_cue_ball(ball_copy);
			}	
		}

		if (deadline_miss(i)) {
			increase_dl_miss_variable();
			printf("deadline missed by the user task\n");
		}
		wait_for_activation(i);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    turn_off() function: quits the game
//--------------------------------------------------------------------------------------------------------------

void	turn_off(void)
{
int	k;			// loop variable
int	x = mouse_x;		// x-component of the mouse position
int	y = mouse_y;		// y-component of the mouse position
static	int mouse_pressed = 0;	// this variable is required to activate functions associated 
 				// with this button only at the first click of the mouse.

	if (!mouse_pressed) {
		if (mouse_b & 1){	// mouse button is being pressed
			if (within_region_test(x, y, H_CLS, V_CLS, BUTTN, BUTTN)) {
				mouse_pressed = 1; 
				set_flag_variable(I_CLS, 1);
			}
		}
	}
	else if (!(mouse_b & 1)) {	// mouse button is being released
		mouse_pressed = 0;
		for (k = 0; k < N_THREADS; k++) pthread_cancel(tid[k]);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    reset() function: restarts the game
//--------------------------------------------------------------------------------------------------------------

void	reset(void)
{
int	k;			// loop variable
int	x = mouse_x;		// x-component of the mouse position
int	y = mouse_y;		// y-component of the mouse position
static	int mouse_pressed = 0;	// this variable is required to activate functions associated 
 				// with this button only at the first click of the mouse.

	if (!mouse_pressed) {
		if (mouse_b & 1) {	// mouse button is being pressed
			if (within_region_test(x, y, H_RST, V_RST, BUTTN, BUTTN)) {
				mouse_pressed = 1;
				set_flag_variable(I_RST, 1);
					
				for (k = 0; k < N_BALLS; k++) pthread_cancel(tid[k]);
				init_balls();
				for (k = 0; k < N_BALLS; k++) {
					create_thread(ball_task, B_PER, B_DL, B_PRIO, k);
				}
			}
		}
	}
	else if (!(mouse_b & 1)) {	// mouse button is being released
		mouse_pressed = 0;
		set_flag_variable(I_RST, 0);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    aim() function: activates the aiming mode
//--------------------------------------------------------------------------------------------------------------

int	aim(void)
{
int	x = mouse_x;		// x-component of the mouse position
int	y = mouse_y;		// y-component of the mouse position
static	int aim_flag = 0;	// flag variable
static	int mouse_pressed = 0;	// this variable is required to activate functions associated 
 				// with this button only at the first click of the mouse.

 	if (!mouse_pressed) {	
		if (mouse_b & 1) {	// mouse button is being pressed
			if (within_region_test(x, y, H_SGT, V_SGT, BUTTN, BUTTN)) {
				mouse_pressed = 1;
				aim_flag = (aim_flag + 1) % 2;
				switch_flag_variable(I_SGT);
			}
		}	
	}
	else if (!(mouse_b & 1)) mouse_pressed = 0;	// mouse button is being released

	return aim_flag;
}

//--------------------------------------------------------------------------------------------------------------
//    hit_cue_ball() function: allows the user to hit the cue ball
//--------------------------------------------------------------------------------------------------------------

void	hit_cue_ball(struct state* b)
{
float	mouse_x_cm;		// x-component of the mouse position converted in cm
float	mouse_y_cm;		// y-component of the mouse position converted in cm
float	scale;			// scale factor
static	int mouse_pressed = 0;	// this variable is required to activate functions associated 
 				// with this button only at the first click of the mouse.

 	if (!mouse_pressed) {
		if (mouse_b & 1) {	// mouse button is being pressed
			mouse_pressed = 1;
			mouse_x_cm = mouse_x * PIX_CM;
			mouse_y_cm = mouse_y * PIX_CM;
			scale = V_SUP / norm(HR * PIX_CM, VR * PIX_CM);
				
			b[0].vx = (b[0].px - mouse_x_cm) * scale;
			b[0].vy = (b[0].py - mouse_y_cm) * scale;

			set_ball_state(0, b[0]);
		}
	}
	else if (!(mouse_b & 1)) mouse_pressed = 0;	// mouse button is being released
}

//--------------------------------------------------------------------------------------------------------------
//    place_cue_ball() function: allows user to place the cue ball on the table
//--------------------------------------------------------------------------------------------------------------

void	place_cue_ball(struct state* b)
{
struct	state zero[TRL_LEN] = {0, 0, 0, 0};	// init value for the variable "trail"
float	dx, dy;			// difference between two position vectors (mouse pointer and cue ball)
float	dist[N_BALLS - 1];	// distance: norm of the vector [dx, dy]
float	dmin;			// minimum distance
float	mouse_x_cm;		// x-component of the mouse position converted in cm
float	mouse_y_cm;		// y-component of the mouse position converted in cm
int	k;			// loop variable

	set_trail(zero);
	set_flag_variable(I_PLC, 1);
	while (!(mouse_b & 1) || mouse_y < WP) {
		mouse_x_cm = mouse_x * PIX_CM;
		mouse_y_cm = mouse_y * PIX_CM;

		if (within_table_test(mouse_x_cm, mouse_y_cm, - DB / 2)) {
	
			for (k = 1; k < N_BALLS; k++) {
				dx = mouse_x_cm - b[k].px;
				dy = mouse_y_cm - b[k].py;
				dist[k - 1] = norm(dx, dy);
			}
			dmin = min(dist, N_BALLS - 1);

			if (dmin >= DB) {
				b[0].px = mouse_x_cm;
				b[0].py = mouse_y_cm;	
				b[0].vx = 0;
				b[0].vy = 0;

				set_ball_state(0, b[0]);
			}
		}
		wait_for_activation(N_BALLS + 1);
	}
	set_flag_variable(I_PLC, 0);
}

//--------------------------------------------------------------------------------------------------------------
//    compute_trail() function: computes positions of trail points. In input receives ball structures.
//--------------------------------------------------------------------------------------------------------------

void	compute_trail(struct state* b)
{
int	k;			// loop variable
int	index = 1;		// index of the array "trail_copy"
int	count = 1;		// counter of the while loop

float	delta_t;		// sample time in seconds
float	acc_copy, dmp_copy;	// copies of global variables "acc" and "dmp"
float	dx, dy;			// difference between two position vectors (cue ball and an other ball)
float	dist[N_BALLS - 1];	// distance: norm of the vector [dx, dy]
float	dmin = DB;		// minimum distance
float	scale;			// scale factor

struct	state sample;		// trail point
struct	state sample_buff;	// buffer variable
struct	state trail_copy[TRL_LEN] = {0, 0, 0, 0}; // copy of the global variable "trail"

	acc_copy = read_acc_variable();
	dmp_copy = read_dmp_variable();

	delta_t = (float) B_PER / 1000;
	scale = V_SUP / norm(HR * PIX_CM, VR * PIX_CM);
		
	sample.px =  b[0].px;
	sample.py =  b[0].py;
	sample.vx = (b[0].px - mouse_x * PIX_CM) * scale;
	sample.vy = (b[0].py - mouse_y * PIX_CM) * scale;

	trail_copy[0] = sample;
	while (dmin >= DB && index < TRL_LEN) {

		if (ball_in_play_test(sample)) {
			sample_buff = inter_collisions(sample, acc_copy, delta_t);
 			sample_buff = wall_collision(sample_buff, dmp_copy);
			sample = in_pocket(sample_buff, 0);
		}
		for (k = 1; k < N_BALLS; k++) {
			dx = sample.px - b[k].px;
			dy = sample.py - b[k].py;
			dist[k - 1] = norm(dx, dy);
		}
		dmin = min(dist, N_BALLS - 1);
	
		if (count % STEP == 0) {	
			trail_copy[index] = sample;
			index++;	
		}
		count++;	
	}
	set_trail(trail_copy);
}

//--------------------------------------------------------------------------------------------------------------
//    modify_acc() function: modifies the value of the global variable "acc"
//--------------------------------------------------------------------------------------------------------------

void	modify_acc(void)
{
int	x = mouse_x;	// x-component of the mouse position
int	y = mouse_y;	// y-component of the mouse position

static	int	count_max = (ACC_MAX - ACC_MIN) / ACC_RES;
static	int	count = (ACC - ACC_MIN) / ACC_RES;
static	int 	button_pressed = 0;

	if (!button_pressed) {
		if (mouse_b & 1) {	// mouse button is being pressed
			if (within_region_test(x, y, H_AW1, V_AW1, BUTTN/2, BUTTN/2)) {
				button_pressed = 1;
				set_flag_variable(I_AW1, 1);
				if (count < count_max) set_acc_variable(++count * ACC_RES + ACC_MIN);
			}
			if (within_region_test(x, y, H_AW1, V_AW1 + BUTTN/2, BUTTN/2, BUTTN/2)) {
				button_pressed = 1;
				set_flag_variable(I_AW1, 2);	
				if (count > 0) set_acc_variable(--count * ACC_RES + ACC_MIN);
			}
		}
	}
	else if (!(mouse_b & 1)) {	// mouse button is being released
			button_pressed = 0;
			set_flag_variable(I_AW1, 0);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    modify_dmp() function: modifies the value of the global variable "dmp"
//--------------------------------------------------------------------------------------------------------------

void	modify_dmp(void)
{
int	x = mouse_x;	// x-component of the mouse position
int	y = mouse_y;	// y-component of the mouse position

static	int	count_max = (DMP_MAX - DMP_MIN) / DMP_RES;
static	int	count = (DMP - DMP_MIN) / DMP_RES;
static	int 	button_pressed = 0;

	if (!button_pressed) {
		if (mouse_b & 1) {	// mouse button is being pressed
			if (within_region_test(x, y, H_AW2, V_AW2, BUTTN/2, BUTTN/2)) {
				button_pressed = 1;
				set_flag_variable(I_AW2, 1);
				if (count < count_max) set_dmp_variable(++count * DMP_RES + DMP_MIN);
			}
			if (within_region_test(x, y, H_AW2, V_AW2 + BUTTN/2, BUTTN/2, BUTTN/2)) {
				button_pressed = 1;
				set_flag_variable(I_AW2, 2);
				if (count > 0) set_dmp_variable(--count * DMP_RES + DMP_MIN);
			}
		}
	}
	else if (!(mouse_b & 1)) {	// mouse button is being released
			button_pressed = 0;
			set_flag_variable(I_AW2, 0);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    AUXILIARY FUNCTIONS
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//    init() function: initializes the game
//--------------------------------------------------------------------------------------------------------------

void	init(void)
{	
	srand(time(0));		// init seed for generation of pseudo-random numbers

	acc = ACC;		// init variable "acc"
	dmp = DMP;		// init variable "dmp"
	dl_miss = 0;		// init variable "dl_miss"

	init_graphics();
	init_mutex();
	init_flags();
	init_balls();
	init_threads();	
}

//--------------------------------------------------------------------------------------------------------------
//    init_flags() function: initializes the variables "flag[i]"
//--------------------------------------------------------------------------------------------------------------

void	init_flags(void)
{
int k;

	for(k = 0; k < N_FLAGS; k++) set_flag_variable(k, 0);
}

//--------------------------------------------------------------------------------------------------------------
//    min() function: returns the minimum element in the array "vect" of "dim" elements
//--------------------------------------------------------------------------------------------------------------

float	min(float* vect, int dim)
{
int	k;
float	min;

	min = vect[0];
	for (k = 1; k < dim ; k++) {
		if (min > vect[k]) min = vect[k];
	}
	return min;
}

//--------------------------------------------------------------------------------------------------------------
//    norm() function: returns the norm of a 2-D vector, receiving its components "x" and "y"
//--------------------------------------------------------------------------------------------------------------

float	norm(float x, float y)
{
float	n = sqrt(x * x + y * y);

	return n;
}

//--------------------------------------------------------------------------------------------------------------
//    dot() function: returns the dot product of two 2-D vectors, receiving their components
//--------------------------------------------------------------------------------------------------------------

float	dot(float x1, float y1, float x2, float y2)
{
float	d = x1 * x2 + y1 * y2;

	return d;
}
//--------------------------------------------------------------------------------------------------------------
//    triang_root() function: returns the triangular root of a triangular number "n"
//--------------------------------------------------------------------------------------------------------------

int	triang_root(int n)
{
int	rt = round((sqrt(1 + 8 * n) - 1) / 2);

	return rt;
}

//--------------------------------------------------------------------------------------------------------------
//    generate_triang_config() function: sets up the input arrays ("x" and "y") with position components
//        of points on the triangular configuration. The triang. config. has dimension "dim".
//--------------------------------------------------------------------------------------------------------------

/* positions on the triang. config. and their corresponding "index + 1" values (case: N_BALLS = 16, dim = 5):

					                 (11)
					             (7)
					        (4)      (12)
					    (2)      (8)
					(1)     (5)      (13)
					    (3)      (9)
					        (6)      (14)        	
					            (10)
					                 (15)
*/


void	generate_triang_config(float* x, float* y, int dim)
{
int	index = 0;	// index of arrays
int	i, j;		// loop variables

	for (j = 0; j < dim; j++) {
		for(i = 0; i <= j; i++) {
			x[index] = DB / 2 * sqrt(3) * j;
			y[index] = DB / 2 * (2 * i - j);
			index++;
		}
	}
}

//--------------------------------------------------------------------------------------------------------------
//    change_balls_sequence() function: changes the layout of the balls on the triangular configuration
//--------------------------------------------------------------------------------------------------------------

void	change_balls_sequence(int* sequence, int dim)
{
int	vert_1, vert_2, centre;		// indices of vertex and centre of the triangular config.
int	striped, solid;			// indices of striped and solid balls
int	n_group = (N_BALLS - 2) / 2;	// number of balls of the same group (solids or stripes)
int	swap_flag;			// flag variable
int	k;				// loop variable
int	n;				// random number in [1, 15]

	centre = 5;			// centre on position (5)
	vert_1 = N_BALLS - dim;		// vertex on position (11)
	vert_2 = N_BALLS - 1;		// vertex on position (15)

	// displacement of a striped ball from position "striped" to position "vertex2", through an exchange
	striped = rand() % n_group + (N_BALLS - n_group);
	exchange(sequence, striped, vert_2);

	// displacement of a solid ball from position "solid" to position "vertex1", through an exchange
	solid = rand() % n_group + 1;
	exchange(sequence, solid, vert_1);

	// exchange between balls on position "vertex1" and on position "vertex2" (random event)
	swap_flag = rand() % 2;
	if (swap_flag) exchange(sequence, vert_1, vert_2);

	// displacement of the ball #8 from position (8) to the central position, through an exchange
	exchange(sequence, n_group + 1, centre);

	// casual placement of remaining balls, through exchanges
	for (k = 1; k < N_BALLS; k++) {
		if (k != centre && k != vert_1 && k != vert_2) {
			n = rand() % (N_BALLS - 1) + 1;
			if (n != centre && n != vert_1 && n != vert_2) {
				exchange(sequence, k, n);
			}
		}
  	}
}

//--------------------------------------------------------------------------------------------------------------
//    exchange() function: exchanges the "i"-th element of the vector "vect" with its "j"-th element
//--------------------------------------------------------------------------------------------------------------

void  exchange(int* vect, int i, int j)
{
int elem_buffer;

	elem_buffer = vect[i - 1];
	vect[i - 1] = vect[j - 1];
	vect[j - 1] = elem_buffer;
}

//--------------------------------------------------------------------------------------------------------------
//    steady_state_test() function: returns 1 if all balls have zero speed, 0 otherwise
//--------------------------------------------------------------------------------------------------------------

int	steady_state_test(struct state* b)
{
int	k;
int	count = 0;

	for (k = 0; k < N_BALLS; k++) {
		if (b[k].vx == 0 && b[k].vy == 0) count++;
	}	
	if (count == N_BALLS) return 1;
	return 0;
}

//--------------------------------------------------------------------------------------------------------------
//    final_state_test() function: returns 1 if all balls have remained pocketed, 0 otherwise
//--------------------------------------------------------------------------------------------------------------

int	final_state_test(struct state* b)
{
int	k;
int	count = 0;

	for (k = 1; k < N_BALLS; k++) {
		if (b[k].py < (WP * PIX_CM)) count++;
	}	
	if (count == N_BALLS - 1) return 1;
	return 0;
}

//--------------------------------------------------------------------------------------------------------------
//   ball_in_play_test() function: returns 1 if ball "b" is not pocketed, 0 otherwise
//--------------------------------------------------------------------------------------------------------------

int	ball_in_play_test(struct state b)
{
	if (b.py > (WP * PIX_CM)) return 1;
	return 0;
}

//--------------------------------------------------------------------------------------------------------------
//    within_table_test() function: returns 1 if the object placed in (x, y) is within a rectangular region
//         defined by the position of table borders + an offset value, 0 otherwise
//--------------------------------------------------------------------------------------------------------------

int	within_table_test(float x, float y, float offset)
{
	if (x >= HP - offset && x <= HP + HD + offset) {
		if (y >= VP - offset && y <= VP + VD + offset) return 1;
	}
	return 0;
}

//--------------------------------------------------------------------------------------------------------------
//    within_region_test() function: returns 1 if the object placed in (x, y) is within a rectangular region
//         defined by its bottom left corner (a, b) and its dimensions (delta_a, delta_b), 0 otherwise
//--------------------------------------------------------------------------------------------------------------

int	within_region_test(int x, int y, int a, int b, int delta_a, int delta_b)
{
	if (x >= a && x <= a + delta_a) {
		if (y >= b && y <= b + delta_b) return 1;
	}
	return 0;
}
