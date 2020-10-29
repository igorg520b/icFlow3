#ifndef FRAMEINFO_H
#define FRAMEINFO_H

#include <iostream>
#include <map>
#include <vector>
#include <string>

namespace icy {

struct FrameInfo
{
public:
    int nNodes;
    int nElems;

    int nActiveNodes;
    int nNonZeroEntries;
    unsigned nodeOffset, elemOffset; // in data file

    // time
    int StepNumber;
    double SimulationTime;
    int TimeScaleFactor;

    unsigned long SimulationIntegerTime;  // // measures time in 1/Parts intervals of InitialTimeStep
    double TimeStep;// time difference between current frame and last frame
    int StepsWithCurrentFactor; // time steps left with current factor (if TSF > 1)

    // solution analysis
    int count_iterations;
    int count_attempts;
    int count_solves; // number of times the solver was invoked (used as a kill switch)};
    double RelativeError;
    double Error0;

    // current progress
    bool solution_reached;
    int solverProgress;

    // benchmarking (in milliseconds)
    long b_total;       // excludes discarded
    long b_prepare;
    long b_clear_ls;
    long b_force_elem;
    long b_force_buoyancy;
    long b_force_collision;
    long b_create_structure;
    long b_assemble;
    long b_solve;
    long b_pull_from_ls;

    // fracture-related
    long b_split;
    long b_local_substep;
    long b_compute_fracture_directions;
    long b_infer_support;
    long b_identify_regions;

    void BenchmarkingClear()
    {
        b_total = 0;
        b_prepare = 0;
        b_clear_ls = 0;
        b_force_elem = 0;
        b_force_buoyancy = 0;
        b_force_collision = 0;
        b_create_structure = 0;
        b_assemble = 0;
        b_solve = 0;
        b_pull_from_ls = 0;

        b_split = 0;
        b_local_substep = 0;
        b_compute_fracture_directions = 0;
        b_infer_support = 0;
        b_identify_regions = 0;
    }

    void Reset() {
        nNodes = nElems = 0;
        nActiveNodes = nNonZeroEntries = 0;
        nodeOffset = elemOffset = 0;
        StepNumber = TimeScaleFactor= 0;
        SimulationTime = 0;
        SimulationIntegerTime = 0;
        TimeStep = 0;
        StepsWithCurrentFactor = 0;
        count_iterations = count_attempts = count_solves = 0;
        RelativeError = Error0 = 0;
        solution_reached = false;
        solverProgress = 0;
        BenchmarkingClear();
    }

    static void BenchmarkingSummarize(std::vector<icy::FrameInfo> &stepStats,
                                      std::vector<std::pair<std::string, long>> &timing_motion,
                                      long &total,
                                      long &n_solves)
    {
        long t_other = 0;

        long t_prepare = 0;
        long t_clear_ls = 0;
        long t_force_elem = 0;
        long t_force_buoyancy = 0;
        long t_force_collision = 0;
        long t_create_structure = 0;
        long t_assemble = 0;
        long t_solve = 0;
        long t_pull_from_ls = 0;

        long t_split = 0;
        long t_local_substep = 0;
        long t_compute_fracture_directions = 0;
        long t_infer_support = 0;
        long t_identify_regions = 0;

        total = 0;
        n_solves = 0;

        for(const icy::FrameInfo &f : stepStats)
        {
            total += f.b_total;
            t_prepare += f.b_prepare;
            t_clear_ls += f.b_clear_ls;
            t_force_elem += f.b_force_elem;
            t_force_buoyancy += f.b_force_buoyancy;
            t_force_collision += f.b_force_collision;
            t_create_structure += f.b_create_structure;
            t_assemble += f.b_assemble;
            t_solve += f.b_solve;
            t_pull_from_ls += f.b_pull_from_ls;

            t_split += f.b_split;
            t_local_substep += f.b_local_substep;
            t_compute_fracture_directions += f.b_compute_fracture_directions;
            t_infer_support += f.b_infer_support;
            t_identify_regions += f.b_identify_regions;

            n_solves += f.count_solves;
        }

        t_other = total - (t_prepare+t_clear_ls+t_force_elem+
                           t_force_buoyancy+t_force_collision+
                           t_create_structure+t_assemble+t_solve+
                           t_split+t_local_substep+t_compute_fracture_directions+
                           t_infer_support+t_identify_regions);

        // t_total is computed for all frames
        // everything else is presented per frame

        t_other /= 1000;
        t_prepare /= 1000;
        t_clear_ls /= 1000;
        t_force_elem /= 1000;
        t_force_buoyancy /= 1000;
        t_force_collision /= 1000;
        t_create_structure /= 1000;
        t_assemble /= 1000;
        t_solve /= 1000;
        t_pull_from_ls /= 1000;

        t_split/=1000;
        t_local_substep/=1000;
        t_compute_fracture_directions/=1000;
        t_infer_support/=1000;
        t_identify_regions/=1000;

        timing_motion.clear();

        timing_motion.push_back({"prep",t_prepare});
        timing_motion.push_back({"clear_ls",t_clear_ls});
        timing_motion.push_back({"elem_frc",t_force_elem});
        timing_motion.push_back({"buoy_frc",t_force_buoyancy});
        timing_motion.push_back({"coll_frc",t_force_collision});
        timing_motion.push_back({"ls_strct",t_create_structure});
        timing_motion.push_back({"assemble",t_assemble});
        timing_motion.push_back({"solve",t_solve});
//        timing_motion.push_back({"pull",t_pull_from_ls));

        timing_motion.push_back({"split",t_split});
        timing_motion.push_back({"substep",t_local_substep});
        timing_motion.push_back({"frac_dir",t_compute_fracture_directions});
        timing_motion.push_back({"support",t_infer_support});
        timing_motion.push_back({"regions",t_identify_regions});

        timing_motion.push_back({"other",t_other});
    }

};
}

#endif // FRAMEINFO_H
