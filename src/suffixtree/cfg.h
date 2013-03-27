#ifndef CFG_H_
#define CFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "elib.h"

#define CFG_MAX_LINE_LEN 400
#define CFG_MAX_VAL_LEN 200
#define COMMENT '#'

int getCfgVal(char *cFile, char *key);

#endif /* CFG_H_ */
