#ifndef GRAPHICS_H
#define GRAPHICS_H

//==============================================================================================================
//    GLOBAL CONSTANTS
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    GRAPHIC COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define COL_DEPTH	32	// color depth

#define HR		1024	// horizontal resolution of the window			[pixel]
#define VR		768     // vertical resolution of the window			[pixel]
#define WP		50      // width of the panel		        		[pixel]

#define H_DMD		9	// number of diamonds on the horizontal rail (included those hidden by covers)
#define V_DMD		5	// number of diamonds on the vertical rail (included those hidden by covers)
#define R_DMD		2	// radius of diamonds
#define R_SP		4	// radius of head and foot spots
#define R_TRL		2	// radius of the last point of the trail

//--------------------------------------------------------------------------------------------------------------
//    PANEL COSTANTS
//--------------------------------------------------------------------------------------------------------------

#define BUTTN	35		// dimension of the button				[pixel]
#define D_SLOT	30		// distance between two ball slots			[pixel]

#define H_CLS	8		// horizontal position of the close button		[pixel]
#define V_CLS	8		// vertical position of the close button  		[pixel]

#define H_RST	48		// horizontal position of the reset button		[pixel]
#define V_RST	8		// vertical position of the reset button		[pixel]

#define H_SGT	88		// horizontal position of the sight button		[pixel]
#define V_SGT	8		// vertical position of the sight button		[pixel]

#define H_AW1	168		// horizontal position of the 1st arrow button		[pixel]
#define V_AW1	8		// vertical position of the 1st arrow button		[pixel]

#define H_AW2	328		// horizontal position of the 2nd arrow button		[pixel]
#define V_AW2	8		// vertical position of the 2nd arrow button  		[pixel]

#define H_DLM	14		// horizontal position of the DL-MISS text		[pixel]
#define V_DLM	(WP + 8)	// vertical position of the DL-MISS text		[pixel]

#define H_ACC	(H_AW1 + 30)	// horizontal position of the print of variable "acc"	[pixel]
#define V_ACC	(V_AW1 + 15)	// vertical position of the print of variable "acc"	[pixel]

#define H_DMP	(H_AW2 + 30)	// horizontal position of the print of variable "dmp"	[pixel]
#define V_DMP	(V_AW2 + 15)	// vertical position of the print of variable "dmp"	[pixel]

#define H_TXT1	(H_AW1 + 30)	// horizontal position of the text "Acceleration"	[pixel]
#define V_TXT1	(V_AW1 + 1)	// vertical position of the text "Acceleration"		[pixel]

#define H_TXT2	(H_AW1 + 30)	// horizontal position of the text "[m/s^2]"		[pixel]
#define V_TXT2	(V_AW1 + 25)	// vertical position of the text "[m/s^2]"		[pixel]

#define H_TXT3	(H_AW2 + 30)	// horizontal position of the text "Damping"		[pixel]
#define V_TXT3	(V_AW2 + 1)	// vertical position of the text "Damping"		[pixel]

#define H_TXT4	531		// horizontal position of the text "PLACE THE CUE..."	[pixel]
#define V_TXT4	(WP + 8)	// vertical position of the text "PLACE THE CUE..."	[pixel]

//==============================================================================================================
//    FUNCTION PROTOTYPES
//==============================================================================================================
//--------------------------------------------------------------------------------------------------------------
//    GRAPHIC FUNCTIONS
//--------------------------------------------------------------------------------------------------------------

void	init_graphics(void);

void	create_colors(void);
void	create_sprites(void);
void	create_frame(void);

void	draw_panel(BITMAP* buff);
void	draw_table(BITMAP* buff);
void	draw_balls(BITMAP* buff);
void	draw_trail(BITMAP* buff);
void	draw_covers(BITMAP* buff);

void	constant_conversion(void);
int	cm_to_pixel(float x);

//--------------------------------------------------------------------------------------------------------------

#endif
