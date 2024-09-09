#include <stdio.h>
#include <math.h>
#include <allegro.h>
#include "task_lib.h"
#include "graphics_lib.h"

//==============================================================================================================
//    GLOBAL VARIABLE DEFINITIONS
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    BITMAPS: variables representing a pointer to a bitmap file
//--------------------------------------------------------------------------------------------------------------

static BITMAP*	sprite_ball[N_BALLS];		// sprites of balls
static BITMAP*	sprite_pocket[N_POCKETS];	// sprites of pockets
static BITMAP*	sprite_cover[N_COVERS];		// sprites of covers

static BITMAP*	sprite_close_0;			// "close" button released
static BITMAP*	sprite_close_1;			// "close" button pressed

static BITMAP*	sprite_reset_0;			// "reset" button released
static BITMAP*	sprite_reset_1;			// "reset" button pressed

static BITMAP*	sprite_sight_0;			// "sight" button released
static BITMAP*	sprite_sight_1;			// "sight" button pressed

static BITMAP*	sprite_arrow_0;			// "arrow" button released
static BITMAP*	sprite_arrow_1;			// "arrow" button pressed (up arrow)
static BITMAP*	sprite_arrow_2;			// "arrow" button pressed (down arrow)

//--------------------------------------------------------------------------------------------------------------
//    COLORS: variables representing a color in the RGB space
//--------------------------------------------------------------------------------------------------------------

static int	white;
static int	red;
static int	grey_1;
static int	grey_2;
static int	grey_3;
static int	grey_4;
static int	green_1;
static int	green_2;
static int	brown_1;
static int	brown_2;

//--------------------------------------------------------------------------------------------------------------
//    CONVERTED GEOMETRIC COSTANTS
//--------------------------------------------------------------------------------------------------------------

static int	wr, wc, pk, db;	// converted values of costants WR, WC, PK, DB (in pixel)
static int	hd, vd, hp, vp;	// converted values of costants HD, VD, HP, VP (in pixel)

//==============================================================================================================
//    FUNCTION DEFINITIONS
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    init_graphics() function:	- set graphics options
//				- installs the input managers
//				- converts costant values
//				- loads all bitmap files
//--------------------------------------------------------------------------------------------------------------

void	init_graphics(void)
{
	set_color_depth(COL_DEPTH);
	set_gfx_mode(GFX_AUTODETECT_WINDOWED, HR, VR, 0, 0);

	constant_conversion();
	create_colors();
	create_sprites();

	install_keyboard();
	install_mouse();
	enable_hardware_cursor();
	show_mouse(screen);
}

//--------------------------------------------------------------------------------------------------------------
//    create_colors() function: creates a set of colors, assigning their respective RGB-components
//--------------------------------------------------------------------------------------------------------------

void	create_colors(void)
{
	white 	= makecol(255, 255, 255);
	red 	= makecol(255,   0,   0);
	grey_1	= makecol( 62,  59,  53);
	grey_2	= makecol( 42,  38,  34);
	grey_3	= makecol( 40,  40,  40);
	grey_4	= makecol(180, 180, 180);
	green_1	= makecol(  0, 120,   0);
	green_2	= makecol(  0,  80,   0);
	brown_1	= makecol( 50,  10,  10);
	brown_2	= makecol(100,  50,   0);
}

//--------------------------------------------------------------------------------------------------------------
//    create_sprites() function: loads all bitmap files contained in folders "balls", "pockets" and "buttons"
//--------------------------------------------------------------------------------------------------------------

void	create_sprites(void)
{
// variables containing the path of a bitmap file
char	  ball_file_name[PATH_LEN] = {'\0'}; 
char	pocket_file_name[PATH_LEN] = {'\0'};
char	 cover_file_name[PATH_LEN] = {'\0'};

int	k;	// loop variable

	for (k = 0; k < N_BALLS; k++) {
		sprintf(ball_file_name, "sprites/balls/ball_%d.bmp", k);
		sprite_ball[k] = load_bitmap(ball_file_name, NULL);
	}

	for (k = 0; k < N_POCKETS; k++) {
		sprintf(pocket_file_name, "sprites/pockets/pocket_%d.bmp", k + 1);
		sprite_pocket[k] = load_bitmap(pocket_file_name, NULL);
	}

	for (k = 0; k < N_COVERS; k++) {
		sprintf(cover_file_name, "sprites/pockets/cover_%d.bmp", k + 1);
		sprite_cover[k] = load_bitmap(cover_file_name, NULL);
	}

	sprite_close_0 = load_bitmap("sprites/buttons/close_0.bmp", NULL);
	sprite_close_1 = load_bitmap("sprites/buttons/close_1.bmp", NULL);

	sprite_reset_0 = load_bitmap("sprites/buttons/reset_0.bmp", NULL);
	sprite_reset_1 = load_bitmap("sprites/buttons/reset_1.bmp", NULL);

	sprite_sight_0 = load_bitmap("sprites/buttons/sight_0.bmp", NULL);
	sprite_sight_1 = load_bitmap("sprites/buttons/sight_1.bmp", NULL);

	sprite_arrow_0 = load_bitmap("sprites/buttons/arrow_0.bmp", NULL);
	sprite_arrow_1 = load_bitmap("sprites/buttons/arrow_1.bmp", NULL);
	sprite_arrow_2 = load_bitmap("sprites/buttons/arrow_2.bmp", NULL);
}	


//--------------------------------------------------------------------------------------------------------------
//    create_frame() function: draws all the graphic entities and copies them on the bitmap variable "screen"
//--------------------------------------------------------------------------------------------------------------

void	create_frame(void)
{
// buffer on which all graphic entities are drawn
BITMAP*	bmp_buff = create_bitmap(HR, VR);

	scare_mouse();	
	clear_to_color(bmp_buff, brown_1);

	draw_panel(bmp_buff);	
	draw_table(bmp_buff);
	draw_trail(bmp_buff);
	draw_balls(bmp_buff);
	draw_covers(bmp_buff);

	// "bmp_buff" value is copied on "screen" variable (to avoid flickering)
	blit(bmp_buff, screen, 0, 0, 0, 0, HR, VR); 
	destroy_bitmap(bmp_buff);

	unscare_mouse();
}

//--------------------------------------------------------------------------------------------------------------
//    draw_panel() function: draws the graphic element "panel" on the input bitmap variable "buff"
//--------------------------------------------------------------------------------------------------------------

void	draw_panel(BITMAP* buff)
{
// local variables containing values read from global variables
int	close, reset, sight, arrw1, arrw2, place, dlm_copy; 
float	acc_copy, dmp_copy;

char	string_buff[STRG_LEN] = {'\0'};	// buffer on which are printed strings
int	k;				// loop variable

	// reading from global variables
	close = read_flag_variable(I_CLS);
	reset = read_flag_variable(I_RST);
	sight = read_flag_variable(I_SGT);
	arrw1 = read_flag_variable(I_AW1);
	arrw2 = read_flag_variable(I_AW2);
	place = read_flag_variable(I_PLC);

	acc_copy = read_acc_variable();
	dmp_copy = read_dmp_variable();
	dlm_copy = read_dl_miss_variable();

	// drawing of the panel background
	rectfill(buff, 0, 0, HR, WP, grey_1);

	// drawings of buttons (in different status: pressed or not pressed)
	if (close == 0)	stretch_sprite(buff, sprite_close_0, H_CLS, V_CLS, BUTTN, BUTTN);
	else		stretch_sprite(buff, sprite_close_1, H_CLS, V_CLS, BUTTN, BUTTN);

	if (reset == 0)	stretch_sprite(buff, sprite_reset_0, H_RST, V_RST, BUTTN, BUTTN);
	else		stretch_sprite(buff, sprite_reset_1, H_RST, V_RST, BUTTN, BUTTN);

	if (sight == 0)	stretch_sprite(buff, sprite_sight_0, H_SGT, V_SGT, BUTTN, BUTTN);
	else		stretch_sprite(buff, sprite_sight_1, H_SGT, V_SGT, BUTTN, BUTTN);

	if (arrw1 == 0)		stretch_sprite(buff, sprite_arrow_0, H_AW1, V_AW1, BUTTN / 2 + 1, BUTTN);
	else if (arrw1 == 1)	stretch_sprite(buff, sprite_arrow_1, H_AW1, V_AW1, BUTTN / 2 + 1, BUTTN);
	else			stretch_sprite(buff, sprite_arrow_2, H_AW1, V_AW1, BUTTN / 2 + 1, BUTTN);

	if (arrw2 == 0)		stretch_sprite(buff, sprite_arrow_0, H_AW2, V_AW2, BUTTN / 2 + 1, BUTTN);
	else if (arrw2 == 1)	stretch_sprite(buff, sprite_arrow_1, H_AW2, V_AW2, BUTTN / 2 + 1, BUTTN);
	else			stretch_sprite(buff, sprite_arrow_2, H_AW2, V_AW2, BUTTN / 2 + 1, BUTTN);

	// drawings of ball slots	
	for (k = 0; k < N_BALLS; k++)
		circlefill(buff, HR - D_SLOT * (N_BALLS - k), WP / 2, db / 2, grey_2);

	// acceleration panel subsection
	sprintf(string_buff, " % 1.2f", acc_copy / 100);
	textout_ex(buff, font, "Acceleration:", H_TXT1, V_TXT1, red,   -1);
	textout_ex(buff, font, string_buff,  H_ACC,  V_ACC,  white, -1);
	textout_ex(buff, font, " [m/s^2]", H_TXT2, V_TXT2, white, -1);

	// damping panel subsection
	sprintf(string_buff, "%3.0f %%", dmp_copy * 100);
	textout_ex(buff, font, "Damping:", H_TXT3, V_TXT3, red,   -1);
	textout_ex(buff, font, string_buff, H_DMP,  V_DMP,  white, -1);

	// deadline-miss panel subsection
	sprintf(string_buff, "DL-MISS: %3d", dlm_copy);
	textout_ex(buff, font, string_buff, H_DLM, V_DLM, white, -1);

	// "PLACE THE CUE BALL ON THE TABLE" warning
	if (place == 1)	{
		sprintf(string_buff, "PLACE THE CUE BALL ON THE TABLE");
		textout_ex(buff, font, string_buff, H_TXT4, V_TXT4, white, -1);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    draw_table() function: draws the graphic element "table" on the input bitmap variable "buff"
//--------------------------------------------------------------------------------------------------------------

void	draw_table(BITMAP* buff)
{
// arrays representing positions of pocket sprites
int	x[N_POCKETS] = {

		hp - wr,
		hp - wr,
		hp - wr + hd / 2,
		hp - wr + hd / 2,
		hp - wr + hd,
		hp - wr + hd
	};

int	y[N_POCKETS] = {

		vp - wr,
		vp - wr + vd,
		vp - wr,
		vp - wr + vd,
		vp - wr,
		vp - wr + vd
	};

int	a, b, w;	// position vector on a rectangular grid (a, b) and offset value (w)
int 	i, j, k;	// loop variables
	
	// drawings of wooden rails, cushions and bed of the table
	rectfill(buff, hp - wr, vp - wr, hp + hd - 1 + wr, vp + vd - 1 + wr, brown_2);	// wooden rails
	rectfill(buff, hp - wc, vp - wc, hp + hd - 1 + wc, vp + vd - 1 + wc, green_1);	// cushions
	rectfill(buff, hp,      vp,      hp + hd - 1,      vp + vd - 1,      green_2);	// bed

	// drawings of diamonds on the wooden rails (rectangular pattern)
	w = (wr + wc) / 2;	// offset: if w == 0, diamonds are drawn on the bed border

	for (i = 0; i < H_DMD; i++) {
		for (j = 0; j < V_DMD; j++) {
			a = hp + hd * i / (H_DMD - 1);
			b = vp + vd * j / (V_DMD - 1);
			if(i == 0)
				circlefill(buff, a - w, b, R_DMD, grey_4); // vertical left side
			if(i == (H_DMD - 1))
				circlefill(buff, a + w, b, R_DMD, grey_4); // vertical right side
			if(j == 0)
				circlefill(buff, a, b - w, R_DMD, grey_4); // horizontal up side
			if(j == (V_DMD - 1)) 
				circlefill(buff, a, b + w, R_DMD, grey_4); // horizontal bottom side
		}
	}

	// drawings of the head spot and the foot spot
	circlefill(buff, hp + 1 * hd / 4, vp + vd / 2, R_SP, grey_3);	// head spot
	circlefill(buff, hp + 3 * hd / 4, vp + vd / 2, R_SP, grey_3);	// foot spot

	// drawings of the six pockets
	for (k = 0; k < N_POCKETS; k++) 
		stretch_sprite(buff, sprite_pocket[k], x[k], y[k], 2 * wr, 2 * wr);
}

//--------------------------------------------------------------------------------------------------------------
//    draw_balls() function: draws the graphic elements "ball" on the input bitmap variable "buff"
//--------------------------------------------------------------------------------------------------------------

void	draw_balls(BITMAP* buff)
{
struct	state ball_copy;	// copy of the global variable "ball[k]"
int	x, y;			// coordinates of the top left corner of the sprite
int	k;			// loop variable
	
	for (k = 0; k < N_BALLS; k++) {

		ball_copy = read_ball_state(k);		
		// assignment of the position
		x = cm_to_pixel(ball_copy.px) - db / 2;
		y = cm_to_pixel(ball_copy.py) - db / 2;

		// drawings of the balls
		stretch_sprite(buff, sprite_ball[k], x, y, db, db);
	}
}

//--------------------------------------------------------------------------------------------------------------
//    draw_covers() function: draws the graphic elements "covers" on the input bitmap variable "buff"
//--------------------------------------------------------------------------------------------------------------

void	draw_covers(BITMAP* buff)
{
// arrays representing positions of cover sprites
int	x[N_COVERS] = {

		hp - wr,
		hp - wr,
		hp - wr + hd / 2,
		hp - wr + hd / 2,
		hp - wr + hd,
		hp - wr + hd
	};

int	y[N_COVERS] = {

		vp - wr,
		vp - wr + vd,
		vp - wr,
		vp - wr + vd,
		vp - wr,
		vp - wr + vd
	};

int	k;

	// drawings of the six covers
	for (k = 0; k < N_COVERS; k++)
		stretch_sprite(buff, sprite_cover[k], x[k], y[k], 2 * wr, 2 * wr);	
}

//--------------------------------------------------------------------------------------------------------------
//    draw_trail() function: draws the graphic element "trail" on the input bitmap variable "buff"
//--------------------------------------------------------------------------------------------------------------

void	draw_trail(BITMAP* buff)
{
struct	state ball_copy[N_BALLS];	// copy of the global variable "ball"
struct	state* trail_copy;		// copy of the global variable "trail"
int	sight;				// copy of the global variable "flag[I_SGT]"
int	blue_component, color;		// color variables
int	x, y;				// position of a trail point, converted in pixel
int	k;				// loop variable

	sight = read_flag_variable(I_SGT);
	if (sight == 1) {
		for (k = 0; k < N_BALLS; k++) ball_copy[k] = read_ball_state(k);

		if (steady_state_test(ball_copy) && ball_in_play_test(ball_copy[0])) {
			trail_copy = read_trail();

			// drawing of the trail
			for(k = 0; k < TRL_LEN; k++) {
				if (within_table_test(trail_copy[k].px, trail_copy[k].py, 0)) {
					x = cm_to_pixel(trail_copy[k].px);
					y = cm_to_pixel(trail_copy[k].py);
					
					if (norm(trail_copy[k].vx, trail_copy[k].vy) > 0) {
						blue_component = 255 * k / (TRL_LEN - 1);
						color = makecol(255, 255, blue_component);
						putpixel(buff, x, y, color);
					}
					// drawing of the final point of the trail
					else circlefill(buff, x, y, R_TRL, white); 
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------------------
//    constant_conversion() function: applies cm_to_pixel() to useful constant
//--------------------------------------------------------------------------------------------------------------

void	constant_conversion(void)
{
	wr = cm_to_pixel(WR);
	wc = cm_to_pixel(WC);
	pk = cm_to_pixel(PK);
	db = cm_to_pixel(DB);
	hd = cm_to_pixel(HD);
	vd = cm_to_pixel(VD);
	hp = cm_to_pixel(HP);
	vp = cm_to_pixel(VP);
}

//--------------------------------------------------------------------------------------------------------------
//    cm_to_pixel() function: converts to pixel a quantity expressed in cm
//--------------------------------------------------------------------------------------------------------------

int	cm_to_pixel(float x)
{
int	n = (int) (x * CM_PIX);

	return n;
} 

