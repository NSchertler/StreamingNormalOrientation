/*
* Copyright (c) 2010-2011, NVIDIA Corporation
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*   * Neither the name of NVIDIA Corporation nor the
*     names of its contributors may be used to endorse or promote products
*     derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! \file priority_queue.h
*   \brief A CUDA-compatible, fixed-size priority queue
*/

#pragma once

#include "types.h"

namespace nih {

	///
	/// A fixed-size priority queue
	///
	template <typename Key, uint32 SIZE>
	struct priority_queue
	{
		/// constructor
		///
		FORCE_INLINE NIH_HOST_DEVICE priority_queue(Key* store);

		/// is queue empty?
		///
		FORCE_INLINE NIH_HOST_DEVICE bool  empty() const;

		/// return queue size
		///
		FORCE_INLINE NIH_HOST_DEVICE uint32 size() const;

		/// push an element
		///
		FORCE_INLINE NIH_HOST_DEVICE void push(const Key key);

		/// pop an element
		///
		FORCE_INLINE NIH_HOST_DEVICE void pop();

		/// top of the queue
		///
		FORCE_INLINE NIH_HOST_DEVICE Key top() const;

		/// return the i-th element in the heap
		///
		FORCE_INLINE NIH_HOST_DEVICE const Key& operator[] (const unsigned int i) const;

		uint32  m_size;
		Key*    m_queue;
	};

} // namespace nih

#include "priority_queue_inline.h"
