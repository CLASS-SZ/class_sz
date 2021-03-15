/**
 * This file should not be modified
 * */

#ifdef __ASSIGN_DEFAULT_TSZ__
#define class_sz_ptsz_parameter(NAME,TYPE,DEF_VALUE)          \
ptsz->NAME = DEF_VALUE;
#endif
#ifdef __ALLOCATE_TSZ_PARAMETER__
#define class_sz_ptsz_parameter(NAME,TYPE,DEF_VALUE)          \
TYPE NAME;
#endif
#ifdef __PARSE_TSZ_PARAMETER__
#define class_sz_ptsz_parameter(NAME,TYPE,DEF_VALUE)          \
class_read_ ## TYPE(#NAME,ptsz->NAME);
#endif


#ifdef __ASSIGN_DEFAULT_TSZ__
#define class_sz_string_parameter(NAME,DIR,STRING)   \
sprintf(ptsz->NAME,__CLASSDIR__);                  \
strcat(ptsz->NAME,DIR);
#endif
#ifdef __ALLOCATE_TSZ_PARAMETER__
#define class_sz_string_parameter(NAME,DIR,STRING)    \
FileName NAME;
#endif
#ifdef __PARSE_TSZ_PARAMETER__
#define class_sz_string_parameter(NAME,DIR,STRING)     \
class_read_string(STRING,ptsz->NAME);
#endif


#ifdef __ASSIGN_DEFAULT_TSZ__
#define class_sz_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
ptsz->NAME = DEF_VAL;
#endif
#ifdef __ALLOCATE_TSZ_PARAMETER__
#define class_sz_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
REAL_TP NAME;
#endif
#ifdef __PARSE_TSZ_PARAMETER__
#define class_sz_type_parameter(NAME,READ_TP,REAL_TP,DEF_VAL) \
class_read_ ## READ_TP(#NAME,ptsz->NAME);
#endif
