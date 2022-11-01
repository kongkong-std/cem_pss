//==========info==========//
/*
 * main function
 */

static const char help[] = "hello, hhhh\
			    test test test"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/main_implement.h"

int main()
{
    puts( "here we go..." );

    puts( "\n==========loading linear system==========" );
    linear_system();

    return 0;
}
