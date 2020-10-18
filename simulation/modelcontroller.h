#ifndef IMPLICITMODEL5_H
#define IMPLICITMODEL5_H

#include <vector>
#include <chrono>
#include <QObject>
#include <QDebug>

#include "model.h"
#include "frameinfo.h"
#include "serializer.h"

namespace icy { class ModelController; }

class icy::ModelController : public QObject
{
    Q_OBJECT
public:
    icy::Model model;
    SimParams prms;
    Serializer serializer;
    std::vector<icy::FrameInfo> stepStats;  // information about each time step
    icy::FrameInfo ts;  // tentative step info (the step itself may fail)
    bool requestToStop = false; // ModelController request backgroundworker to stop calling Step()

    void SaveAs(QString fileName);
    void Load(QString fileName);        // load initial setup from file
    void ImportFloePatch(QString fileName);
    void Remesh();
    void Reset();                       // reset simulation to pristine state; erase model
    void GoToStep(int step);            // load data from file for presentation
    void Trim();                        // remove subsequent steps
    void Prepare();                     // compute constant matrices - call once before first Step()
    void Step();                        // perform one computation step
    void Fracture();
    void RequestAbort();                // cancel current step; invoked by GUI thread

    // progress summary
    int getTotalSteps() {
        std::size_t n = stepStats.size();
        if(n>0) n--;
        return n;
    }
    int getCurrentStep() { return current_step; }

private:
    unsigned current_step = 0;
    bool abort_requested = false;
    void Aborting();       // called before exiting Step() if aborted

    // parts of Step();
    void _BeginStep();  // make initial guess, initialize FrameInfo
    void _Write();

signals:
    void stepCompleted();
    void fractureUpdated();
    void stepAborted();
    void progressUpdated();
};

#endif // IMPLICITMODEL5_H
