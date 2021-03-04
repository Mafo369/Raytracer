#ifndef COLOR_MSG_H
#define COLOR_MSG_H

#define BLACK           "\033[0;30m"
#define DARK_GREY       "\033[1;30m"
#define RED             "\033[0;31m"
#define LIGHT_RED       "\033[1;31m"
#define GREEN           "\033[0;32m"
#define LIGHT_GREEN     "\033[1;32m"
#define BROWN_ORANGE    "\033[0;33m"
#define YELLOW          "\033[1;33m"
#define BLUE            "\033[0;34m"
#define LIGHT_BLUE      "\033[1;34m"
#define PURPLE          "\033[0;35m"
#define LIGHT_PURPLE    "\033[1;35m"
#define CYAN            "\033[0;36m"
#define LIGHT_CYAN      "\033[1;36m"
#define LIGHT_GRAY      "\033[0;37m"
#define WHITE           "\033[1;37m"
#define NC              "\033[0m"


#define logc(color, file, message, ...) fprintf(file, color message "\033[0m", ##__VA_ARGS__)


#endif
