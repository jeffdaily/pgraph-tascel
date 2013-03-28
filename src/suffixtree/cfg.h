/**
 * @file cfg.h
 *
 * @author andy.cj.wu@gmail.com
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2010 Washington State University. All rights reserved.
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 */
#ifndef CFG_H_
#define CFG_H_

#define CFG_MAX_LINE_LEN 400
#define CFG_MAX_VAL_LEN 200
#define COMMENT '#'

int get_config_val(const char *config_file, const char *key);

#endif /* CFG_H_ */
