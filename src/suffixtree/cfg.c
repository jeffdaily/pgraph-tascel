#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "cfg.h"
#include "elib.h"

int get_config_val(const char *config_file, const char *key)
{
    FILE *fp = NULL;
    char line[CFG_MAX_LINE_LEN];
    int val;
    char *p = NULL;

    fp = efopen(config_file, "r");

    while (fgets(line, CFG_MAX_LINE_LEN, fp)) {
        /* comment line starts with '#' */
        if (strchr(line, COMMENT)) continue;

        /* empty line */
        if (line[0] == '\n') continue;

        /* config line */
        if (strstr(line, key)) {
            p = line;
            while(isspace(*p++));

            /* %*s used to skip the first item */
            sscanf(p, "%*s %d\n", &val);
            fclose(fp);
            return val;
        }  
    }

    fclose(fp);
    printf("cannot find config value for key: [%s]\n", key);

    return -1;
}
