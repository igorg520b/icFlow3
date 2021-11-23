#ifndef TESTSET_SIMPLEOBJECTPOOL_H
#define TESTSET_SIMPLEOBJECTPOOL_H

#include <vector>

namespace icy {
    template<class T>
    class SimpleObjectPool;
}


template<class T>
class icy::SimpleObjectPool {

public:
    SimpleObjectPool(int initialSize);
    ~SimpleObjectPool();
    SimpleObjectPool& operator=(SimpleObjectPool&) = delete;

    T* take();
    void release(T* obj);
    void release(std::vector<T*> &vec);
    void releaseAll();

private:
    std::vector<T*> available;      // items that are free to use
    std::vector<T*> registry;       // all items of the pool
};

template<class T>
icy::SimpleObjectPool<T>::SimpleObjectPool(int initialSize)
{
    available.reserve(initialSize*2);
    registry.reserve(initialSize*2);

    for(int i=0;i<initialSize;i++)
    {
        T* obj = new T;
        available.push_back(obj);
        registry.push_back(obj);
    }
}

template<class T>
icy::SimpleObjectPool<T>::~SimpleObjectPool()
{
    for(auto &x : registry) delete x;
}

template <class T>
T* icy::SimpleObjectPool<T>::take()
{
    T *obj;
    if(available.size() == 0) {
        obj = new T;
        registry.push_back(obj);
    }
    else {
        obj = available.back();
        available.pop_back();
    }
    return obj;
}

template<class T>
void icy::SimpleObjectPool<T>::release(T* obj)
{
    available.push_back(obj);
}

template<class T>
void icy::SimpleObjectPool<T>::release(std::vector<T*> &vec)
{
    for(T* p : vec) available.push_back(p);
    vec.clear();
}


template<class T>
void icy::SimpleObjectPool<T>::releaseAll()
{
    available.clear();
    available.insert(available.end(), registry.begin(), registry.end());
}



#endif //TESTSET_SIMPLEOBJECTPOOL_H
