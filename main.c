//=================================================================================================
//    POOL GAME
//
//    Real Time Systems
//    Roberto Mauceri
//    2018/19
//=================================================================================================

#include <pthread.h>
#include <allegro.h>
#include "task_lib.h"

int	main()
{
	allegro_init();
	init();
	wait_for_termination();
	allegro_exit();
	
	return 0;	
}
