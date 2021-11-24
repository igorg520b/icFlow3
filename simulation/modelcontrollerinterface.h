#ifndef MODELCONTROLLERINTERFACE_H
#define MODELCONTROLLERINTERFACE_H

#include <string>

class ModelControllerInterface
{
public:
    virtual void Prepare() = 0;
    virtual bool Step() = 0;
    virtual void RequestAbort() = 0;

    virtual void GoToStep(int step) = 0;
    virtual void SaveAs(std::string fileName) = 0;
    virtual void Load(std::string fileName) = 0;
};


#endif // MODELCONTROLLERINTERFACE_H
