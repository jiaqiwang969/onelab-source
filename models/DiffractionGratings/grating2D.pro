// Copyright (C) 2020 Guillaume Demésy
//
// This file is part of the model grating2D.pro.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with This program. If not, see <https://www.gnu.org/licenses/>.


Include "grating2D_data.geo";
Include "grating2D_materials.pro";
If (flag_polar<2)
    Include "grating2D_scalar.pro";
Else
    Include "grating2D_conical.pro";
EndIf
