// Link-specific controls keep the Release public driver object identical across
// original Leopard2, runtime OFF/ON, and traced runtime qualification links.
#ifndef LEO_ADJACENT_PUBLIC_LINK_H
#define LEO_ADJACENT_PUBLIC_LINK_H
#include "avx2_adjacent_control.h"
const char* LeoAdjacentPublicIdentity();
bool LeoAdjacentPublicScheduleValid(const char* schedule);
bool LeoAdjacentPublicSelect(char state);
unsigned LeoAdjacentPublicState(); // 0/1 runtime, 2 original/native
bool LeoAdjacentPublicDefault();
bool LeoAdjacentPublicMeasurable();
bool LeoAdjacentPublicTraced();
void LeoAdjacentPublicReset();
LeoAdjacentCounts LeoAdjacentPublicCounts();
#endif
