/*******************************************************************************
 * include/tools/mem_monitor.h
 *
 * Copyright (C) 2016-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef _MEM_MONITOR_H_
#define _MEM_MONITOR_H_

#include "sys/types.h"
#include "sys/sysinfo.h"

unsigned long printRamUsage(){

    struct sysinfo memInfo;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double kb = 1024.0;
    const double mb = kb*1024;
    [[maybe_unused]] const double gb = mb*1024;

    sysinfo (&memInfo);
    long long totalVirtualMem = memInfo.totalram;
    //Add other values in next statement to avoid int overflow on right hand side...
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;

    long long totalPhysMem = memInfo.totalram;
    //Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;

    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;

    unsigned long long freeRam = memInfo.freeram;
    freeRam *= memInfo.mem_unit;

    unsigned long long sharedRam = memInfo.sharedram;
    sharedRam *= memInfo.mem_unit;    
    
    unsigned long long buffRam = memInfo.bufferram;
    buffRam *= memInfo.mem_unit; 

    auto parseLine = [](char* line){
        // This assumes that a digit will be found and the line ends in " Kb".
        int i = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') p++;
        line[i-3] = '\0';
        i = atoi(p);
        return i;
    };

    auto getValue = [&](){ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmRSS:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    };

    std::cout<< rank <<  ": totalPhysMem: " << (totalPhysMem/mb) << 
                " MB, physMemUsed: " << physMemUsed/mb << 
                " MB, free ram: " << freeRam/mb <<
                " MB, shared ram: " << sharedRam/mb <<
                " MB, buffered   ram: " << buffRam/mb <<
                " MB, I am using: " << getValue()/kb << " MB" << std::endl;

    return freeRam;
}

#endif