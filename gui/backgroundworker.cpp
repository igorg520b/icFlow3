#include "backgroundworker.h"

BackgroundWorker::BackgroundWorker(ModelControllerInterface *controller_) : controller(controller_)
{
    this->start();
}

// resume the worker thread
void BackgroundWorker::Resume()
{
    controller->Prepare();
    condition.wakeOne();
}

// cancel current step and pause the worker thread
void BackgroundWorker::Pause()
{
    if(!running) return;
    timeToPause = true;
    controller->RequestAbort();
}

// exit the worker thread
void BackgroundWorker::Finalize()
{
    controller->RequestAbort();
    kill=true;
    condition.wakeOne();
    bool result = wait();
    qDebug() << "BackgroundWorker::Finalize() terminated; " << result;
}

void BackgroundWorker::run()
{
    while(!kill)
    {
        if (timeToPause)
        {
            timeToPause = false;
            running = false;
            Q_EMIT workerPaused();
            mutex.lock();
            condition.wait(&mutex);
            mutex.unlock();
            running = true;
        }
        if(kill) break;

        bool result = controller->Step();
        if(!result) timeToPause = true;
        Q_EMIT stepCompleted();
    }
    qDebug() << "BackgroundWorker::run() terminated";
}
