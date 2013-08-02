/**
 * @file constants.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define DOLLAR 'U'
#define BEGIN 'O'  /* it is oh, not zero */
#define SIGMA 26U /**< size of alphabet */
#define NROW 2U /**< number of rows in dynamic programming table, always 2 */

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

enum {ERROR = -1, DOL_END = -2};
enum {FALSE = 0, TRUE = 1, MAYBE = 2};
enum {NO = 0, YES = 1};

#endif /* _CONSTANTS_H_ */
