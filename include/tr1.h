/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief poor man's way of getting tr1::shared pointer. This doesn't
 * fully comply with tr1, but it surves the need until C++11 supported
 * by compilers.
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

#ifndef _TR1_H_
#define _TR1_H_

#include <cstddef>
#include "Logger.h"

namespace tr1{

    /**
     * place-holder class for the true shared_ptr from tr1. If boost
     * is available or the compiler supports C++11, those
     * implementation should be used.
     *
     * @warning This not an efficient implemenation of shared_ptr. It
     * is intended to be used only for heavy objects with long life
     * span.
     */
    template<typename T>
    class shared_ptr_base{
      public:
        typedef T element_type;

	explicit shared_ptr_base():
            ptr_(NULL),fore(NULL),aft(NULL) {}

        explicit shared_ptr_base(T *p):
            ptr_(p),fore(NULL),aft(NULL) {}

        shared_ptr_base(shared_ptr_base const &r):
            ptr_(r.ptr_)
        { join(&r); }

        shared_ptr_base& operator=(shared_ptr_base const &r){
            leave();
            join(&r);
            ptr_ = r.ptr_;
	    return *this;
        }

        const shared_ptr_base& operator=(shared_ptr_base const &r) const{
            leave();
            join(&r);
            ptr_ = r.ptr_;
	    return *this;
        }

        shared_ptr_base& operator=(T* p){
            leave();
            ptr_ = p;
	    return *this;
        }

        const shared_ptr_base& operator=(T* p) const{
            leave();
            ptr_ = p;
	    return *this;
        }

        virtual ~shared_ptr_base()
        {
            COUTDEBUG("Destructing shared_ptr to "<<this->ptr_);
            leave();
        }

        T* operator->() {return ptr_;}
        const T* operator->() const {return ptr_;}
	T* get() {return ptr_;}
	const T* get() const {return ptr_;}

	operator bool() const
        { return(bool(ptr_)); }

      protected:
        mutable T* ptr_;
        mutable shared_ptr_base const* fore;
        mutable shared_ptr_base const* aft;

	void join(shared_ptr_base const* before) const{
            fore      = before;
            aft       = before->aft;
            fore->aft = this;

            if (aft)
                aft->fore = this;
        }

	void leave() const{
            if (fore)
                fore->aft = aft;

            if (aft)
                aft->fore = fore;

            if ( is_last() && ptr_){
                COUTDEBUG("Deleting pointer to "<<ptr_);
		delete ptr_;
            } else {
		ptr_ = NULL;
	    }
        }

	bool is_last() const{
	    return (fore == NULL && aft == NULL);
	}
    };

    template<typename T>
    class shared_ptr : public shared_ptr_base<T>{
      public:
	explicit shared_ptr(){};
        explicit shared_ptr(T *p) : shared_ptr_base<T>(p) {}
	shared_ptr(shared_ptr const &r) : shared_ptr_base<T>(r) {}

        shared_ptr& operator=(shared_ptr const &r){
	    shared_ptr_base<T>::operator=(r);
	    return *this;
        }

	const shared_ptr& operator=(shared_ptr const &r) const{
	    shared_ptr_base<T>::operator=(r);
	    return *this;
        }

        shared_ptr& operator=(T* p){
	    shared_ptr_base<T>::operator=(p);
	    return *this;
        }

        const shared_ptr& operator=(T* p) const{
	    shared_ptr_base<T>::operator=(p);
	    return *this;
        }

        T& operator*() {return *shared_ptr_base<T>::ptr_;}
	const T& operator*() const {return *shared_ptr_base<T>::ptr_;}
    };

    template<>
    class shared_ptr<void> : public shared_ptr_base<void>{
      public:
	explicit shared_ptr(){};
        explicit shared_ptr(void *p) : shared_ptr_base<void>(p) {}
	shared_ptr(shared_ptr const &r) : shared_ptr_base<void>(r) {}

        shared_ptr& operator=(shared_ptr const &r){
	    shared_ptr_base<void>::operator=(r);
	    return *this;
        }

        const shared_ptr& operator=(shared_ptr const &r) const{
	    shared_ptr_base<void>::operator=(r);
	    return *this;
        }

        shared_ptr& operator=(void* p){
	    shared_ptr_base<void>::operator=(p);
	    return *this;
        }

        const shared_ptr& operator=(void* p) const{
	    shared_ptr_base<void>::operator=(p);
	    return *this;
        }
    };
};

#endif //_TR1_H_
