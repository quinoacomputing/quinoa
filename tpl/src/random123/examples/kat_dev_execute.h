/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "kat.h"
KAT_KERNEL void dev_execute_tests(KAT_GLOBAL kat_instance *tests, KAT_UINT ntests){
    size_t i;
    for(i=0; i<ntests; ++i){
        KAT_GLOBAL kat_instance *ti = &tests[i];
        switch(tests[i].method){
            //case philox2x32_e: ti->philox2x32_data.computed = philox2x32_R(ti->rounds, ti->philox2x32_data.ctr, ti->philox2x32_data.key);
#define RNGNxW_TPL(base, N, W) case base##N##x##W##_e: ti->u.base##N##x##W##_data.computed = base##N##x##W##_R(ti->nrounds, ti->u.base##N##x##W##_data.ctr, base##N##x##W##keyinit(ti->u.base##N##x##W##_data.ukey)); break;
#include "rngNxW.h"
#undef RNGNxW_TPL
        case unused: ; // silence a warning
        }
    }
}
