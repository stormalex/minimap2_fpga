#ifndef __USER_COMMON_H__
#define __USER_COMMON_H__

#define debug 1
#define warning 1


#if debug
#define DEBUG(format, ...)          fprintf(stderr, "\033[40;34mDEBUG\033[0m: [%s %d]"format"\n", __FUNCTION__, __LINE__, ##__VA_ARGS__);
#else
#define DEBUG(format, ...)
#endif

#if warning
#define WARNING(format, ...)        fprintf(stderr, "\033[40;33mWARNING\033[0m: [%s %d]"format"\n", __FUNCTION__, __LINE__, ##__VA_ARGS__);
#else
#define WARNING(format, ...)
#endif

#define ERROR(format, ...)          fprintf(stderr, "\033[40;31mERROR\033[0m: [%s %d]"format"\n", __FUNCTION__, __LINE__, ##__VA_ARGS__);


#endif //__USER_COMMON_H__