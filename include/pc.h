
#ifndef LSSP_PC_H
#define LSSP_PC_H

#include "pc-iluk.h"
#include "pc-ilut.h"
#include "pc-biluk.h"
#include "pc-bilut.h"
#include "pc-vbiluk.h"
#include "pc-vbilut.h"
#include "pc-arms.h"
#include "pc-amg.h"
#include "pc-mi20.h"
#include "pc-sxamg.h"

void lssp_pc_create(LSSP_PC &pc, LSSP_PC_TYPE type);
void lssp_pc_destroy(LSSP_PC &pc);
void lssp_pc_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
