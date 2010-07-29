/**
 * @file   WorkSpace.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 15:39:53 2010
 * 
 * @brief  The declaration for the work space manager.
 */

#ifndef _WORKSPACE_H_
#define _WORKSPACE_H_

#include<stack>

template<typename stackType> 
class WorkSpace
{
    WorkSpace();
    WorkSpace(int stack_size_in);
    ~WorkSpace();

    stackType* Pop();
    stackType* Pop(int n_in);
    void Push(stackType* work_space_in);
    void Push(stackType* work_space_in, int n_in);
    void Expand(int factor_in = 2);
    void Shrink(int factor_in = 2);
    
    void Extend(int n_in);
    void Curtail(int n_in);
};

#include "WorkSpace.cc"
#endif
