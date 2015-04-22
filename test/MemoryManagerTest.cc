/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test
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

#include "MemoryManager.h"
#include "Device.h"
#include "Logger.h"

typedef Device<CPU> DevCPU;
extern const DevCPU cpu_dev(0);
typedef memory::MemoryManager<DevCPU, cpu_dev> CPUMem;

namespace memory{
    template<typename T, typename MM>
    class MemoryTester
    {
      public:
	typedef T value_type;
	typedef memory::ResourcePtr<T, MM> pointer;

	bool RunAll(MM& mgr){
	    bool retval(false);
	    if (test_allocate(mgr)    &&
		test_deallocate(mgr)  &&
		test_ResourcePtr(mgr)){
		COUT(emph<<"* all tests for MemoryManager passed"<<emph);
	    }

	    return retval;
	}

	bool test_allocate(MM& mgr){

	    T** ref;
	    {
		pointer a;
		mgr.allocate(10, &a);
		ASSERT(typeid(a.get())==typeid(T*), "Bad type" );
		ASSERT(a.is_last(), "more than one pointer");
		ref = &a.ptr_;
	    }
	    ASSERT(*ref==NULL, "pointer isn't set to null");

	    COUT("* test_allocate passed");
	    return true;
	}

	bool test_deallocate(MM& mgr){
	    T* ref;
	    pointer a;
	    mgr.allocate(10, &a);
	    ref = a.get();
	    // dealloc shouldn't be called on ResourcePtr
	    // directy. Done only for testing.
	    mgr.deallocate(&a);
	    ASSERT(a.get()==ref, "dealloc sets the pointer to null");
	    a.ptr_ = NULL; /* to avoid segfault */

	    COUT("* test_deallocate passed");
	    return true;
	}

	bool test_ResourcePtr(MM& mgr){
	    T* p(NULL);
	    T **pa, **pb;
	    {
		pointer b;
		{
		    pointer a;
		    mgr.allocate(10, &a);
		    p = a.get();
		    pa = &a.ptr_;
		    ASSERT(b.get()==NULL, "Bad pointer");
		    b = a;
		    ASSERT(a.is_last()==false, "check number of pointers");
		    ASSERT(b.is_last()==false, "check number of pointers");
		}
		ASSERT(b.is_last(), "check number of pointers");
		ASSERT(b.get()==p, "bad pointer");
		ASSERT(*pa==NULL, "bad pointer");
		pb = &b.ptr_;
	    }
	    ASSERT(*pb==NULL, "bad pointer");
	    COUT("* test_ResourcePtr passed");
	    return true;
	}

    };
}

int main(int, char**){
    CPUMem& MM(CPUMem::instance());
    memory::MemoryTester<double,CPUMem> T;
    T.RunAll(MM);
}
