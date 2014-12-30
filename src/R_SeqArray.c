// ===========================================================
//     _/_/_/   _/_/_/  _/_/_/_/    _/_/_/_/  _/_/_/   _/_/_/
//      _/    _/       _/             _/    _/    _/   _/   _/
//     _/    _/       _/_/_/_/       _/    _/    _/   _/_/_/
//    _/    _/       _/             _/    _/    _/   _/
// _/_/_/   _/_/_/  _/_/_/_/_/     _/     _/_/_/   _/_/
// ===========================================================
//
// R_SeqArray.c: link to R_GDS2.h
//
// Copyright (C) 2014-2015    Xiuwen Zheng
//
// This file is part of CoreArray.
//

#include <R_GDS.h>

// do not modify this file, R_GDS2.h is from the gdsfmt package

#if (defined(GDSFMT_R_VERSION) && (GDSFMT_R_VERSION>=0x010103))
#   include <R_GDS2.h>
#else
#   include <R_GDS.c>
#endif
