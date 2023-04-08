
#ifndef LSSP_PC_VBILUT_H
#define LSSP_PC_VBILUT_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"
#include "pc-bilut.h"

#if USE_ITSOL

void lssp_pc_vbilut_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
#endif
