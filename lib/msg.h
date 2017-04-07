#ifndef MSG_H
#define MSG_H 1

enum LogLevel {msg_debug=0, msg_verbose, msg_info, msg_warn, msg_error, msg_fatal, msg_silent};


/*
#ifdef NDEBUG
#define assert_double(x, y, eps)
#else
#define assert_double(x, y, eps) \
  ((void) (msg_assert_double(__FILE__, __LINE__, x, y, eps)))
#endif

void assert_double(const char file[], const unsigned int line, const double x, const double x_expected, const double eps);
*/

//void msg_init(void);
void msg_set_loglevel(const LogLevel level);
void msg_set_prefix(const char prefix[]);

void msg_printf(const LogLevel level, char const * const fmt, ...);
void msg_abort(char const * const fmt, ...);


#endif
