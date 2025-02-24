#ifndef TENSOR_H
#define TENSOR_H 1

/*****************************************************************************
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale

*****************************************************************************
$Header: /repo/nml/include/tensor.h,v 1.3 2003/10/08 16:19:25 hook Exp $
*****************************************************************************/


#include "vector.h"

#define d_1l(p) *int_ptr((p-1))
#define d_1h(p) *int_ptr((p-2))
#define d_2l(p) *int_ptr((p-3))
#define d_2h(p) *int_ptr((p-4))
#define d_3l(p) *int_ptr((p-5))
#define d_3h(p) *int_ptr((p-6))

#endif /* TENSOR_H */
