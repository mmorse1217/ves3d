/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _MEMORYMANAGER_H_
#define _MEMORYMANAGER_H_

#include <memory>
#include "Error.h"
#include "tr1.h"

namespace memory{

    template<typename T, typename MM>
    class MemoryTester;

    template<typename T,typename RM>
    class ResourcePtr : public tr1::shared_ptr<T>
    {
      public:
	friend MemoryTester<T,RM>;

        explicit ResourcePtr() : tr1::shared_ptr<T>() {}
        explicit ResourcePtr(T *p) : tr1::shared_ptr<T>(p) {}
        ResourcePtr(ResourcePtr const &r) : tr1::shared_ptr<T>(r) {}

        ResourcePtr& operator=(ResourcePtr const &r){
	    tr1::shared_ptr<T>::operator=(r);
	    return *this;
        }

        const ResourcePtr& operator=(ResourcePtr const &r) const{
	    tr1::shared_ptr<T>::operator=(r);
	    return *this;
        }

        ResourcePtr& operator=(T* p){
	    tr1::shared_ptr<T>::operator=(p);
	    return *this;
        }

        const ResourcePtr& operator=(T* p) const{
	    tr1::shared_ptr<T>::operator=(p);
	    return *this;
        }

	~ResourcePtr(){
            if (tr1::shared_ptr<T>::is_last())
		free();
	}

      protected:
	void free() const{
	    if (tr1::shared_ptr<T>::ptr_ != NULL){
		COUTDEBUG("Freeing "<<tr1::shared_ptr<T>::ptr_);
		RM::instance().deallocate(this);
		tr1::shared_ptr<T>::ptr_ = NULL;
	    }
      	}
    };

    template<typename DT, const DT &DEVICE>
    class MemoryManager{
      public:
      	typedef void value_type;
	typedef void* raw_pointer;
      	typedef size_t size_type;

	// frequently used types
	typedef ResourcePtr<void  , MemoryManager> void_ptr;
	typedef ResourcePtr<char  , MemoryManager> char_ptr;
	typedef ResourcePtr<int   , MemoryManager> int_ptr;
	typedef ResourcePtr<float , MemoryManager> float_ptr;
	typedef ResourcePtr<double, MemoryManager> double_ptr;

	static MemoryManager& instance();

	template<typename T>
	Error_t allocate(size_type n, ResourcePtr<T, MemoryManager> *p) const;

	template<typename T>
	Error_t deallocate(const ResourcePtr<T, MemoryManager> *p) const;

	size_type max_size() const;

      protected:
	static MemoryManager instance_;

      private:
      	MemoryManager(){}
	MemoryManager(MemoryManager&){}
    };

    template<typename DT, const DT &DEVICE>
    MemoryManager<DT,DEVICE> MemoryManager<DT,DEVICE>::instance_;

    template<typename DT, const DT &DEVICE>
    MemoryManager<DT,DEVICE>& MemoryManager<DT,DEVICE>::instance()
    {
    	return MemoryManager::instance_;
    }

    template<typename DT, const DT &DEVICE>
    template<typename T>
    Error_t MemoryManager<DT,DEVICE>::allocate(size_type n, ResourcePtr<T, MemoryManager> *p) const
    {
	void *mem(DEVICE.Malloc(n * sizeof(T)));
	*p = (T*) mem;
	return ErrorEvent::Success;
    }

    template<typename DT, const DT &DEVICE>
    template<typename T>
    Error_t MemoryManager<DT,DEVICE>::deallocate(const ResourcePtr<T, MemoryManager> *p) const
    {
	DEVICE.Free(const_cast<raw_pointer>((raw_pointer) p->get()));
	return ErrorEvent::Success;
    }

    template<typename DT, const DT &DEVICE>
    MemoryManager<DT,DEVICE>::size_type MemoryManager<DT,DEVICE>::max_size() const
    { return 0;}
}

#endif //_MEMORYMANAGER_H_
