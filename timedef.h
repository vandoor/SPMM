#pragma once
#ifndef __linux__
#include <windows.h>
#define time GetTickCount
#endif