#include "backgroundworker.h"

BackgroundWorker::BackgroundWorker(icy::ModelController *_model) : model(_model)
{
    this->start();
}

// resume the worker thread
void BackgroundWorker::Resume()
{
    condition.wakeOne();
}

// cancel current step and pause the worker thread
void BackgroundWorker::Pause()
{
    if(!running) return;
    timeToPause = true;
    model->RequestAbort();
}

// exit the worker thread
void BackgroundWorker::Finalize()
{
    qDebug() << "BackgroundWorker::Finalize()";
    model->RequestAbort();
    kill=true;
    condition.wakeOne();
    wait();
    qDebug() << "BackgroundWorker::Finalize() terminated";
    running = false;
}

void BackgroundWorker::run()
{
    while(!kill) {
        if (timeToPause) {
            mutex.lock();
            running = false;
            emit workerPaused();
            condition.wait(&mutex);
            timeToPause = false;
            running = true;
            mutex.unlock();
        }
        if(kill) break;
        model->Step();
        if(model->getCurrentStep()==model->prms.MaxSteps || model->requestToStop) timeToPause = true;
    }
    qDebug() << "BackgroundWorker::run() terminated";
}
