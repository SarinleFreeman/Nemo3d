#!/usr/bin/env perl

# *****************************************************************************
# The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
# Copyright (C) 2002 California Institute of Technology (Caltech)
#
# This application is free software, which you can redistribute and/or modify
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; see the file COPYING. If not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA  02111-1307  USA
#
# For additional information, please contact
#   Gerhard Klimeck (gekco@jpl.nasa.gov)
#   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
#
# Written by:  Chris Bowen
#              Gerhard Klimeck
#              Fabiano Oyafuso
#              Seungwon Lee
#              Olga Lazarenkova
#              Hook Hua
#	       Marta Prada
#
# This product includes software developed by the Apache Software Foundation
# (http://www.apache.org/).
#
# *****************************************************************************
# $Header: /repo/nemo3d/misc/update-nemo3d-input.pl,v 1.7 2007/04/23 23:31:38 musman Exp $
# *****************************************************************************
 
################################################################################
#                           update-nemo3d-input.pl
################################################################################
#
#   DESCRIPTION:
#   This script updates the grammar of a NEMO-3D input file to that of the
#   most current version by successively updating a legal xml input for
#   version i to a legal xml input for version i+1 until the current version
#   is reached.
#
#   FORMAT:  update-nemo3d-input.pl <input file>
#
#   OPTIONS: none
#
#   HOW to MODIFY when a NEW NEMO3D INPUT is RELEASED:  Two changes are 
#   necessary and sufficient.  First one must append the most current version
#   ID to the global array @versionList.  Second, a new subroutine must be
#   written to propagate from the previous latest version of the input file
#   to the current latest version.  IMPORTANT:  The SAME NAMING CONVENTION should
#   be used (i.e. update_version_XXX, where XXX is the version ID with dots
#   replaced by underscores.  ANOTHER IMPORTANT NOTE:  ALL changes to the
#   grammar of the input must be accounted for!  These changes encompass
#   everything that is not merely a change in <value>, including changes
#   to <desc>, <cTag>, etc.
#
#   ATTENTION!!: Update to the 0.16-th version changes the VALUE of 
#                strain_VFF array and increases its length from 3 to 6 elements.
#                This is important change of the value, which affects the work 
#                of the code.
#
#   ALGORITHM NOTE:  The insertion and deletion algorithms used here were
#   designed to facilitate programming, but are horribly inefficient 
#   (~O(n^2) for a each operation).  If you have large inputs, you may
#   want to rewrite these subroutines.
#
#   FINAL IMPORTANT NOTE:  This script WILL FAIL if data structures are
#   modified that are also arrays of arbitray length (e.g. Shape, Material, etc.).
#   The reason for this failure is that the code only does string subsitution 
#   (insertion or deletion), but has no concept of a class.  Therefore, it 
#   cannot identify all objects (class instances) of a class.  If you want
#   this feature, you must modify the code yourself.
#
################################################################################



#-------------------------------------------------------------------------------
# list of functions that may need to be appended when changes in NEMO3D are made
#-------------------------------------------------------------------------------

@versionList = ("0.12", "0.13", "0.14", "0.15", "0.16", "0.17", "0.18", "0.19");


sub update_version_0_12 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);

  my($ins1) = 
"                <param type=\"boolean\">
                    <name>Calculate strain distribution?</name>
                    <cTag>CalcStr</cTag>
                    <desc> Calculate strain distribution?</desc>
                    <value>true</value>
                </param>
";

  my($ins2) =
"                        <param type=\"boolean\">
                            <name>Update lattice constants before strain calculation?</name>
                            <cTag>UpdateLatt</cTag>
                            <desc>update lattice constants before strain calculation?</desc>
                            <value>false</value>
                        </param>
                        <param type=\"array_real\">
                            <name>If UpdateLatt is true, what are the new lattice constants?</name>
                            <cTag>UpdateLattConst</cTag>
                            <desc>If UpdateLatt is true, what are the new lattice constants?</desc>
                            <value>{ 0.56532 0.56532 0.56532 }</value>
                        </param>
                        <param type=\"real\">
                            <name>What is the scaling factor for dE/da?</name>
                            <cTag>SCALE_FCTR</cTag>
                            <desc>What is the scaling factor (S) for dE/da?
                                  Use number between 1e-2 to 1e-5 depending 
                                  on the number of atoms (N) in strain calculation.
                                  An ideal S is about 1.0/sqrt(N). 
                            </desc>
                            <value>1e-3</value>
                        </param>
";

  my($ins3) =
"                        <param type=\"integer\">
                            <name>number of eigenvalues requested</name>
                            <cTag>NumEigVal</cTag>
                            <desc>number of eigenvalues requested</desc>
                            <value>6</value>
                        </param>
";

  my($ins4) = 
"                        <param type=\"boolean\">
                            <name>calculate eigenvectors</name>
                            <cTag>CalEigVec</cTag>
                            <desc>calculate eigenvectors?</desc>
                            <value>true</value>
                        </param>
";

  my($ins5) =
"                            <desc>reference eigenvalue</desc>
";

  update_version_ID($fil,$fil_out);

  @start = ("ExecParam","<param");
  @end = ("param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("CalcPhonDisp","<param");
  @end = ("param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("Strain calculation","<option");
  @end = ("option>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("ExecParam","<desc");
  insert_txt($fil_out, \@strList, $ins1);

  @strList = ("MinWrtLatt","</param");
  insert_txt($fil_out, \@strList, $ins2);

  @strList = ("PhonParpack","</param");
  insert_txt($fil_out, \@strList, $ins3);

  @strList = ("PhonParpack","Tolerance","</param");
  insert_txt($fil_out, \@strList, $ins4);

  @start = ("PhonParpack","reference eigenvalue","<desc");
  @end = ("desc>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("PhonParpack","reference eigenvalue","</cTag");
  insert_txt($fil_out, \@strList,$ins5);


}


sub update_version_0_13 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);

  my($ins1)=
"                    <param type=\"option\">
                        <name>Strain Model</name>
                        <cTag>StrainModel</cTag>
                        <desc>strain model</desc>
                        <value>VFF_keating_strained</value>
                        <option>None</option>
                        <option>VFF_keating</option>
                        <option>VFF_keating_strained</option>
                        <option>Read</option>
                    </param>
";

  my($ins2) =
"                            <desc>reference eigenvalue, eV</desc>
";

  update_version_ID($fil,$fil_out);

  @start = ("Strain calculation","<param");
  @end = ("param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("PhonParpack","MaxIter","<param");
  @end = ("param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("PhonParpack","Tolerance","<param");
  @end = ("param>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("Strain calculation","<desc");
  insert_txt($fil_out, \@strList, $ins1);

  @start = ("PhonParpack","reference eigenvalue","<desc");
  @end = ("desc>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("PhonParpack","reference eigenvalue","</cTag");
  insert_txt($fil_out, \@strList,$ins2);
}


sub update_version_0_14 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);

  my($ins1)=
"                        <desc>Print Atom-Type Information</desc>
";

  my($ins2)=
"                    <param type=\"boolean\">
                        <name>NbrIndx</name>
                        <cTag>NbrIndx</cTag>
                        <desc>Print Neighbor Index</desc>
                        <value>true</value>
                    </param>
";

  update_version_ID($fil,$fil_out);

  @start = ("AtomInfo","<desc");
  @end = ("desc>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("AtomInfo","</cTag");
  insert_txt($fil_out, \@strList,$ins1);

  @strList = ("AtomInfo","</param");
  insert_txt($fil_out, \@strList,$ins2);
}


sub update_version_0_15 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);

  my($ins_dummy)=
"                                <value>{ 7.61  4.78  6.45  7.61  4.78  6.45 }</value>
";
  my($ins_GaAs)=
"                        <value>{7.2 7.62 6.4 5.38 7.62 6.83}</value>
";
  my($ins_InAs)=
"                        <value>{7.61 4.78 6.45 4.27 4.78 6.95}</value>
";
  my($ins_InGaAs)=
"                        <value>{ 7.446  5.916  6.43  4.714  5.916  6.91 }</value>
";
  my($ins_Si)=
"                        <value>{7.67 4.67 4.53  4.77 4.67 5.34 }</value>
";
  my($ins_Ge)=
"                        <value>{8.33 4.80 5.08 4.82 4.80 6.18}</value>
";
  my($ins_SiGe)=
"                        <value>{8.085  4.735  4.74 4.795  4.735  5.76 }</value>
";
  my($ins_PS)=
"                        <desc>What part of the phonon spectrum to compute?</desc>
";
  my($ins_anh)=
"                    <param type=\"option\">
                        <name>Optimization of the anharmonicity corrections</name>
                        <cTag>anh</cTag>
                        <inputTag>anh</inputTag>
                        <desc>By default the Lower part of the phonon
                        spectrum computes with the same set of
                        anharmonicity corrections as that used for
                        strain calculations (first 3 values of strain_VFF), while
                        the Higher part of the phonon spectrum
                        computes with the set of anharmonicity
                        corrections optimised for the optical phonons
                        in Gamma point (second 3 values of strain_VFF)</desc>
                        <value>Default</value>
                        <option>Default</option>
                        <option>AcousticalSet</option>
                        <option>OpticalSet</option>
                        <option>LinearInterpolation</option>
                    </param>
";

  update_version_ID($fil,$fil_out);

  @start = ("Add a material<","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material<","VFF constants");
  insert_txt($fil_out, \@strList,$ins_dummy);

  @start = ("material_1","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1","VFF constants");
  insert_txt($fil_out, \@strList,$ins_GaAs);

  @start = ("material_2","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2","VFF constants");
  insert_txt($fil_out, \@strList,$ins_InAs);

  @start = ("material_3","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_3","VFF constants");
  insert_txt($fil_out, \@strList,$ins_InGaAs);

  @start = ("material_4","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_4","VFF constants");
  insert_txt($fil_out, \@strList,$ins_Si);

  @start = ("material_5","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_5","VFF constants");
  insert_txt($fil_out, \@strList,$ins_Ge);

  @start = ("material_6","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_6","VFF constants");
  insert_txt($fil_out, \@strList,$ins_SiGe);

  @start = ("material_7","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_7","VFF constants");
  insert_txt($fil_out, \@strList,$ins_SiGe);

  @start = ("material_8","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_8","VFF constants");
  insert_txt($fil_out, \@strList,$ins_GaAs);

  @start = ("material_9","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_9","VFF constants");
  insert_txt($fil_out, \@strList,$ins_InAs);

  @start = ("material_10","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_10","VFF constants");
  insert_txt($fil_out, \@strList,$ins_InGaAs);

  @start = ("material_11","VFF constants","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_11","VFF constants");
  insert_txt($fil_out, \@strList,$ins_GaAs);

  @start = ("PS","<desc");
  @end = ("desc>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("PS","</inputTag>");
  insert_txt($fil_out, \@strList,$ins_PS);

  @strList = ("PS","</param>");
  insert_txt($fil_out, \@strList,$ins_anh);

}

sub update_version_0_16 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);

  my($Crystal)=
"                    <desc>crystal structure</desc>
";

  my($Band_Model)=
"                    <value>Bands_20_sp3d5ss_spin</value>
                    <option>Bands_1_s_nospin</option>
                    <option>Bands_10_sp3ss_spin</option>
                    <option>Bands_10_sp3d5ss_nospin</option>
                    <option>Bands_20_sp3d5ss_spin</option>
                </param>
";

  my($Strain_Const)=
"                            <option>AllFree</option>
                            <option>FixOne</option>
                            <option>FixSurfaceAll</option>
                            <option>FixInterior</option>
                            <option>FixZmin</option>
                        </param>
";

  my($Band_Dummy)=
"                                <value>Bands_10_sp3ss_spin</value>
";

  my($Band_20)=
"                        <value>Bands_20_sp3d5ss_spin</value>
";

  my($Band_10)=
"                        <value>Bands_10_sp3ss_spin</value>
";


  my($Band_1)=
"                        <value>Bands_1_s_nospin</value>
";

  my($ShapeInfo)=
"                    <param type=\"boolean\">
                        <name>ShapeInfo</name>
                        <cTag>ShapeInfo</cTag>
                        <desc>Print Shape Information</desc>
                        <value>false</value>
                    </param>
";

  my($New_Materials)=
"                <group type=\"obj\">
                    <name>material_12</name>
                    <cTag>Material</cTag>
                    <desc>Si bulk properties (spds* without spin)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>Si</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.42116</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>-0.00305492</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0665337</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.402511</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>12.8464</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.543095</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>48.5</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>13.8</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.67 4.67 4.53  4.77 4.67 5.34 }</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{4.6637}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_10_sp3d5ss_nospin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{ -2.14463 4.2363 -2.14463 4.2363
                            19.12352 19.12352 13.79658 13.79658 0.00000
                            0.00000 -1.83493 -4.36635 -1.16274 -1.16274
                            2.90562 2.90562 2.79465 2.79465 -2.42485
                            -2.42485 -1.02993 -1.02993 4.10364 -1.51801
                            -1.35554 -1.35554 2.38479 2.38479 -1.68136
                            2.58880 -1.81400 0.543095  1 1 1 }</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.13203 0.19237 0.56247 2.36548 0.49435
                            1.84385 2.5672 0.34492 2.23636 4.5125
                            2.00203 4.66836 2.30238 0.92391 1.08601 27
                            27 0.83942 0.38885 0.89216 0.89216 0.24004
                            0.24004 1.7944 1.7944 0 0 0.17104 0.17104
                            2.0332 0 0 2.64866 }</value>
                    </param>
                </group>
                <group type=\"obj\">
                    <name>material_13</name>
                    <cTag>Material</cTag>
                    <desc>CdSe bulk properties (sps*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>CdSe</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>CdSe</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.7</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>0.0</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.12</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.45</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>10.16</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>unstrnd_cubic_cell_length</desc>
                        <value>0.43</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>NONE</desc>
                        <value>1.0</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>NONE</desc>
                        <value>1.0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{ 1.0 1.0 1.0 }</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{18.6655 13.1111}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_10_sp3ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{ -9.63 1.47 0.03 4.73 7.53
                                  5.72 -4.64 2.64 5.36 4.57
                                  5.54 3.05 2.49 0.16 0.0755 }</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0.0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{  2 2 2 2 2 2 2 2 2    2 2 2 2 2 0 27.0
                            27.0 1                             1 1 1 1
                            1 1 1 1 1 1 1 1 1 1 0 }</value>
                    </param>
                </group>
";

  update_version_ID($fil,$fil_out);

  @start = ("Crystal structure","<desc");
  @end = ("/desc>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("CrystalStruct");
  insert_txt($fil_out, \@strList, $Crystal);

  @start = ("Band Model","<value");
  @end = ("/param>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("Band Model","<desc");
  insert_txt($fil_out, \@strList,$Band_Model);

  @start = ("Strain constraint","<option>");
  @end = ("/param>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("Strain constraint","</value>");
  insert_txt($fil_out, \@strList,$Strain_Const);

  @strList = ("AtomPos", "</param>");
  insert_txt($fil_out, \@strList, $ShapeInfo);

  @start = ("Add a material<","band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("Add a material<","band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_Dummy);

  @start = ("material_1<","GaAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1<","GaAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_2<","InAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2<","InAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_3<", "InGaAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_3<", "InGaAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_4<", "Si", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_4<", "Si", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_5<", "Ge", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_5<", "Ge", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_6<", "SiGe", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_6<", "SiGe", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_7<", "GeSi", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_7<", "GeSi", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_20);

  @start = ("material_8<", "GaAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_8<", "GaAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_10);

  @start = ("material_9<", "InAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_9<", "InAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_10);

  @start = ("material_10<", "InGaAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_10<", "InGaAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_10);

  @start = ("material_11<", "GaAs", "band model","<value");
  @end = ("value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_11<", "GaAs", "band model", "</desc>");
  insert_txt($fil_out, \@strList,$Band_1);
  
#  print $New_Materials;

  @strList = ("material_11<", "mb_strain", "</group>");
  insert_txt($fil_out, \@strList, $New_Materials);

}
#  Update May '06

sub update_version_0_17 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);


 my($tilted_hkl)=
"                <param type=\"boolean\">
                    <name>include tilt angle in a quantum well?</name>
                    <cTag>tilt_hkl</cTag>
                    <desc>hkl defined below, describe the direction of the tilt: hkl=1 0 0 means no tilt and 28 0 1 means a ~2-DEG tilt in plane XZ with respect to a_z</desc>
                    <value>false</value>
                </param>
                <param type=\"integer\">
                    <name>tilt_h</name>
                    <cTag>tilt_h</cTag>
                    <desc>tilt_h</desc>
                    <value>1</value>
                </param> 
                <param type=\"integer\">
                    <name>tilt_k</name>
                    <cTag>tilt_k</cTag>
                    <desc>tilt_k</desc>
                    <value>0</value>
                </param>
                <param type=\"integer\">
                    <name>tilt_l</name>
                    <cTag>tilt_l</cTag>
                    <desc>tilt_l</desc>
                    <value>0</value>
                </param>
                <param type=\"integer\">
                    <name>NaddALY</name>
                    <cTag>NaddALY</cTag>
                    <desc># of additional layers to include in y-direction: use this option if you want to describe an structure with a fraction of the unit cell, which would be measured in units of eg. a_y/4 for Zincblende</desc>
                    <value>0</value>
                </param>
                <param type=\"integer\">
                    <name>NaddALZ</name>
                    <cTag>NaddALZ</cTag>
                    <desc># of additional layers to include in z-direction</desc>
                    <value>0</value>
                </param>
";
 my($new_mats)=
"                <group type=\"obj\">
                    <name>material_14</name>
                    <cTag>Material</cTag>
                    <desc>AlAs bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>AlAs</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.42116</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>-0.00305492</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0665337</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.402511</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>12.8464</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.56532</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>NONE</desc>
                        <value>41.49</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>NONE</desc>
                        <value>8.94</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.2 7.62 6.4 5.38 7.62 6.83}</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{11.5773 12.441}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{ -5.17012 4.39708 0.79695 6.63291 19.80474 24.16587
                                13.13880 12.92122 0.17386 0.01589 -1.64584 -2.84245
                                -1.88341 -2.78690 3.02223 2.95309 1.92174 1.30469
                                -3.03196 -2.64111 -1.84300 -1.73510 4.53156 -1.86816
                                -2.47345 -1.02836 2.52741 2.86419 -1.97058 1.67733
                                -1.58868 0.56532 1 1 1 }</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.00000 0.21266 2.06001 1.38498 2.68497
                            1.31405 1.89889 1.39930 1.81235 2.37964
                            2.93686 1.72443 1.97253 1.89672 1.78540 27.0
                            27.0                             0.58696
                            0.48609 0.88921 0.77095 0.75979
                            1.45891 0.81079 1.21202 1.07015 0.38053
                            1.03256 1.31726 0.00000 1.61350 0.00000
                            1.26262 }</value>
                    </param>
                </group>
                <group type=\"obj\">
                    <name>material_15</name>
                    <cTag>Material</cTag>
                    <desc>InSi bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>InSi</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                   </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>0.594213</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>0.224299</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0270249</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.409741</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>14.6</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.56214</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>NONE</desc>
                        <value>35.18</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>NONE</desc>
                        <value>5.49</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.61 4.78 6.45 4.27 4.78 6.95}</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{19.0663 12.441}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{-3.722525 4.293685 -1.241435 6.09932   19.517055
                                19.043995 13.514135 13.7219  0.196115   0.17978
                                -1.58257 -4.1786 -1.7399 -1.26218 3.04547
                                2.640445  2.5684  1.6978 -2.34781 -2.14968
                                -0.85373 -0.6238  4.3008 -1.46767 -1.58044
                                -1.74144  2.4078  2.4533 -1.51364  2.63045
                                -1.13406  0.56214  1  1 1}</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.00000 0.34902 1.87635 1.25042 2.73407
                            1.32992 1.65530 1.11454 2.32014 2.57895
                            1.94235 2.00323 1.84572 2.07788 3.00000 27.0
                            27.0                             0.25331
                            0.00000 2.97410 2.99811 0.93102
                            3.00000 0.00000 2.24414 0.00001 3.00000
                            0.00000 0.00000 2.36738 0.15676 0.49725
                            0.16244 }</value>
                    </param>
                </group>
               <group type=\"obj\">
                    <name>material_16</name>
                    <cTag>Material</cTag>
                    <desc>SiAs bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>SiAs</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>0.594213</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>0.224299</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0270249</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.409741</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>14.6</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.56214</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>NONE</desc>
                        <value>35.18</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>NONE</desc>
                        <value>5.49</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.61 4.78 6.45 4.27 4.78 6.95}</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{19.0663 12.441}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{-4.281823 4.27948  -1.174687  6.33547 19.6149
                                20.283836 13.38665 13.43075 0.221523 0.16045
                                -1.57007 -4.01147  -1.895893 -1.279756  2.918623
                                2.74707   2.371103  1.474436 -2.426396 -2.20665
                                -0.77855 -0.371473  4.250823 -1.454263 -1.6783867
                                -1.79084  2.448313  2.48525  -1.432413  2.5887533
                                -1.03995  0.56214  1  1 1}</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.00000 0.34902 1.87635 1.25042 2.73407
                            1.32992 1.65530 1.11454 2.32014 2.57895
                            1.94235 2.00323 1.84572 2.07788 3.00000 27.0
                            27.0                             0.25331
                            0.00000 2.97410 2.99811 0.93102
                            3.00000 0.00000 2.24414 0.00001 3.00000
                            0.00000 0.00000 2.36738 0.15676 0.49725
                            0.16244 }</value>
                    </param>
                </group>
                <group type=\"obj\">
                    <name>material_17</name>
                    <cTag>Material</cTag>
                    <desc>SiP bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>SiP</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.42116</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>-0.00305492</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0665337</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.402511</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>12.8464</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.543095</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>48.5</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>13.8</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.67 4.67 4.53  4.77 4.67 5.34 }</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{4.6637 4.6637}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{ -5.159 3.3462 -2.14463 4.2363
                            20.49575 19.12352 12.2064 13.79658 0.02485
                            0.01989 -1.7226 -3.6074 -1.6241 -1.53325
                            2.711 2.8802 2.2259 2.24215 -2.65465
                            -2.2415 -0.7118 -0.69 4.11715 -1.37085
                            -1.73305 -1.81725 1.8667 2.1578 -1.28735
                            2.3251 -1.9294 0.543095  1 1 1 }</value>
                    </param>
                    <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.13203 0.19237 0.56247 2.36548 0.49435
                            1.84385 2.5672 0.34492 2.23636 4.5125
                            2.00203 4.66836 2.30238 0.92391 1.08601 27
                            27 0.83942 0.38885 0.89216 0.89216 0.24004
                            0.24004 1.7944 1.7944 0 0 0.17104 0.17104
                            2.0332 0 0 2.64866 }</value>
                    </param>
        </group>
        <group type=\"obj\">
                    <name>material_18</name>
                    <cTag>Material</cTag>
                    <desc>PSi bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>PSi</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                    <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.42116</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>-0.00305492</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0665337</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.402511</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>12.8464</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.543095</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                         <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>48.5</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>13.8</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.67 4.67 4.53  4.77 4.67 5.34 }</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{4.6637 4.6637}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{ -2.14463 4.2363 -5.159 3.3462
                            19.12352 20.49575 13.79658 12.2064 0.01989
                            0.02485 -1.7226 -3.6074 -1.53325 -1.6241
                            2.8802 2.711 2.24215 2.2259 -2.2415
                            -2.65465 -0.69 -0.7118 4.11715 -1.37085
                            -1.81725 -1.73305 2.1578 1.8667 -1.28735
                            2.3251 -1.9294 0.543095  1 1 1 }</value>
                    </param>
                     <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.13203 0.19237 0.56247 2.36548 0.49435
                            1.84385 2.5672 0.34492 2.23636 4.5125
                            2.00203 4.66836 2.30238 0.92391 1.08601 27
                            27 0.83942 0.38885 0.89216 0.89216 0.24004
                            0.24004 1.7944 1.7944 0 0 0.17104 0.17104
                            2.0332 0 0 2.64866 }</value>
                    </param>
                </group>
                <group type=\"obj\">
                    <name>material_19</name>
                    <cTag>Material</cTag>
                    <desc>P bulk properties (spds*)</desc>
                    <param type=\"deletebranch\">
                        <name>Delete This Material</name>
                        <cTag>DUMMY</cTag>
                        <value>dumval</value>
                        <desc>None</desc>
                    </param>
                    <param type=\"string\">
                        <name>Material</name>
                        <cTag>material</cTag>
                        <desc>material</desc>
                        <value>P</value>
                    </param>
                    <param type=\"string\">
                        <name>comp1_mat</name>
                        <cTag>comp1_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>x concentration</name>
                        <cTag>xval</cTag>
                        <desc>x concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"string\">
                        <name>comp2_mat</name>
                        <cTag>comp2_mat</cTag>
                        <desc>not sure about this</desc>
                        <value>NULL</value>
                    </param>
                    <param type=\"real\">
                        <name>y concentration</name>
                        <cTag>yval</cTag>
                        <desc>y concentration</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>temperature</name>
                        <cTag>tkel</cTag>
                        <desc>temperature</desc>
                        <value>300</value>
                    </param>
                     <param type=\"string\">
                        <name>substrate</name>
                        <cTag>substrate</cTag>
                        <desc>substrate</desc>
                        <value>GaAs</value>
                    </param>
                    <param type=\"real\">
                        <name>Ec</name>
                        <cTag>Ec</cTag>
                        <desc>NONE</desc>
                        <value>1.42116</value>
                    </param>
                    <param type=\"real\">
                        <name>Ev</name>
                        <cTag>Ev</cTag>
                        <desc>NONE</desc>
                        <value>-0.00305492</value>
                    </param>
                    <param type=\"real\">
                        <name>electron dos mass</name>
                        <cTag>mstar_c_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.0665337</value>
                    </param>
                    <param type=\"real\">
                        <name>hole dos mass</name>
                        <cTag>mstar_v_dos</cTag>
                        <desc>NONE</desc>
                        <value>0.402511</value>
                    </param>
                    <param type=\"real\">
                        <name>dielectric constant</name>
                        <cTag>epsilon</cTag>
                        <desc>NONE</desc>
                        <value>12.8464</value>
                    </param>
                    <param type=\"real\">
                        <name>unstrained lattice constant</name>
                        <cTag>unstrnd_cubic_cell_length</cTag>
                        <desc>NONE</desc>
                        <value>0.543095</value>
                    </param>
                    <param type=\"real\">
                        <name>alpha (Keating)</name>
                        <cTag>strain_alpha</cTag>
                         <desc>The set of values for alpha and beta
                           chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>48.5</value>
                    </param>
                    <param type=\"real\">
                        <name>beta (Keating)</name>
                        <cTag>strain_beta</cTag>
                        <desc>The set of values for alpha and beta
                            chosen because they correspond most closely
                            to the elastic constants used in nemo [see
                            Chen/Sher p.47]</desc>
                        <value>13.8</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Strain dependence of VFF</name>
                        <cTag>strain_VFF</cTag>
                        <desc>Parameters describing strain dependence of
                            VFF constants</desc>
                        <value>{7.67 4.67 4.53  4.77 4.67 5.34 }</value>
                    </param>
                    <param type=\"array_real\">
                        <name>M_atom</name>
                        <cTag>M_atom</cTag>
                        <desc>Mass of atomic species in 10^-26 kg</desc>
                        <value>{4.6637}</value>
                    </param>
                    <param type=\"string\">
                        <name>band model</name>
                        <cTag>band_model</cTag>
                        <desc>NONE</desc>
                        <value>Bands_20_sp3d5ss_spin</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_pv</name>
                        <cTag>mb_pv</cTag>
                        <desc>NONE</desc>
                        <value>{-5.3367 3.3462 -5.3367 3.3462
                            20.49575 20.49575 14.1074 14.1074 0.02485
                            0.02485  -1.83493 -4.36635 -1.16274 -1.16274
                            2.90562 2.90562 2.79465 2.79465 -2.42485
                            -2.42485 -1.02993 -1.02993 4.10364 -1.51801
                            -1.35554 -1.35554 2.38479 2.38479 -1.68136
                            2.58880 -1.81400 0.543095  1 1 1 }</value>
                    </param>
                      <param type=\"real\">
                        <name>VB offset (eV)</name>
                        <cTag>Ev_offset</cTag>
                        <desc>NONE</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>mb_strain</name>
                        <cTag>mb_strain</cTag>
                        <desc>NONE</desc>
                        <value>{ 0.13203 0.19237 0.56247 2.36548 0.49435
                            1.84385 2.5672 0.34492 2.23636 4.5125
                            2.00203 4.66836 2.30238 0.92391 1.08601 27
                            27 0.83942 0.38885 0.89216 0.89216 0.24004
                            0.24004 1.7944 1.7944 0 0 0.17104 0.17104
                            2.0332 0 0 2.64866 }</value>
                    </param>
                </group>
";
 my($new_GaAs_param) =
"                        <value>{-5.50042 4.15107 -0.24119 6.70776
                            19.71059 22.66352 13.03169 12.74846 0.17234
                            0.02179 -1.64508 -3.704550 -2.20777 -1.31491
                            2.66493 2.96032 1.97650 1.02755 -2.609400
                            -2.32059 -0.62820 0.13324 4.15080 -1.42744
                            -1.87428 -1.88964 2.52926 2.54913 -1.26996
                            2.50536 -0.85174 0.56532 1 1 1 }</value>
";
 my($new_GaAs_mb_strain) =
"                        <value>{ 0.00000 0.21266 2.06001 1.38498 2.68497
                            1.31405 1.89889 1.39930 1.81235 2.37964
                            2.93686 1.72443 1.97253 1.89672 1.78540 27.0
                            27.0                             0.58696
                            0.48609 0.88921 0.77095 0.75979
                            1.45891 0.81079 1.21202 1.07015 0.38053
                            1.03256 1.32385 2.00000 1.61350 1.50000
                            1.26262 }</value>
";
 my($new_InAs_param) = 
"                        <value>{ -5.50042 4.15107 -0.581930 6.97163
                            19.71059 19.94138 13.03169 13.30709 0.17234
                            0.139670 -1.69435 -4.21045 -2.42674 -1.159870
                            2.598230 2.809360 2.06766 0.93734 -2.26837
                            -2.293090 -0.899370 -0.488990 4.310640 -1.288950
                            -1.731410 -1.978420 2.188860 2.456020 -1.584610
                            2.71793 -0.505090 0.60583 1 1 1 }</value>
";
my($new_InAs_mb_strain) =
"                        <value>{ 0.060800 0.000810 1.92494 1.570030 2.06151
                            1.602470 1.765660 1.949370 2.383820 2.455600
                            1.94235 2.32291 1.84572 2.3296 2.02387 27.0
                            27.0                             1.258286
                            2.481447 1.086223 4.557774 4.367575
                            7.02966 3.298598 7.029496 0.0000 0.187036
                            1.195042 1.769483 4.624483 0.0000 0.0000
                            0.246999 }</value>
";
 my($new_Si_param) = 
"                        <value>{ -2.15168 4.22925 -2.15168 4.22925 19.11650
                                 19.11650 13.78950 13.78950 0.01989 0.01989
                                 -1.95933 -4.24135 -1.52230 -1.52230 3.02562
                                 3.02562 3.15565 3.15565 -2.28485 -2.28485
                                 -0.80993 -0.80993 4.10364 -1.51801 -1.35554
                                 -1.35554 2.38479 2.38479 -1.68136 2.58880
                                 -1.81400 0.543095 1 1 1}</value>
";
 my($new_solver_options)=
"                            <option>PQTraceminSolver</option>
                            <option>PCTraceminSolver</option>
                            <option>BlockLanczos</option>
";
 my($tracemin_options)=
"                        <group
                            enabledif=\"/application/input/top/Execution parameters/Electronic structure calculation/Algorithms/Electronic eigensolver==PQTraceminSolver\" type=\"struct\">
                            <name>PQTracemin</name>
                            <cTag>PQTracemin</cTag>
                            <desc>PQTracemin parameters</desc>
                            <param type=\"integer\">
                                <name>number of eigenvalues requested</name>
                                <cTag>NumEigVal</cTag>
                                <desc>number of eigenvalues requested</desc>
                                <value>6</value>
                            </param>
                            <param type=\"real\">
                                <name>reference eigenvalue</name>
                                <cTag>RefEigVal</cTag>
                                <desc>reference eigenvalue</desc>
                                <value>1.12</value>
                            </param>
                            <param type=\"real\">
                                <name>tolerance</name>
                                <cTag>Tolerance</cTag>
                                <desc>tolerance</desc>
                                <value>5e-07</value>
                            </param>
                            <param type=\"real\">
                                <name>initial tol of GC</name>
                                <cTag>cgtol</cTag>
                                <desc>cgtol</desc>
                                <value>1.0</value>
                            </param>
                            <param type=\"integer\">
                                <name>Max. # of iterations</name>
                                <cTag>maxi</cTag>
                                <desc>maxi</desc>
                                <value>40</value>
                            </param>
                            <param type=\"integer\">
                                <name>Number of iterations of CG</name>
                                <cTag>cgmaxim</cTag>
                                <desc>cgmaxim</desc>
                                <value>100</value>
                            </param>
                            <param type=\"boolean\">
                                <name>calculate eigenvectors</name>
                                <cTag>CalEigVec</cTag>
                                <desc>calculate eigenvectors?</desc>
                                <value>true</value>
                            </param>
                        </group>
                        <group
                            enabledif=\"/application/input/top/Execution parameters/Electronic structure calculat
ion/Algorithms/Electronic eigensolver==PCTraceminSolver\" type=\"struct\">
                            <name>PCTracemin</name>
                            <cTag>PCTracemin</cTag>
                            <desc>PCTracemin parameters</desc>
                            <param type=\"integer\">
                                <name>number of eigenvalues requested</name>
                                <cTag>NumEigVal</cTag>
                                <desc>number of eigenvalues requested</desc>
                                <value>6</value>
                            </param>
                            <param type=\"real\">
                                <name>reference eigenvalue</name>
                                <cTag>RefEigVal</cTag>
                                <desc>reference eigenvalue</desc>
                                <value>0.0</value>
                            </param>
                            <param type=\"real\">
                                <name>tolerance</name>
                                <cTag>Tolerance</cTag>
                                <desc>tolerance</desc>
                                <value>5e-07</value>
                            </param>
                            <param type=\"array_real\">
                                <name>Interval in which the eigenvalues are wanted</name>
                                <cTag>WantEigRange</cTag>
                                <desc>Interval [a,b] in which the eigenvalues are wanted</desc>
                                <value>{ -1.0 1.6 }</value>
                            </param>
                            <param type=\"array_real\">
                                <name>Interval that includes all eigenvalues</name>
                                <cTag>AllEigRange</cTag>
                                <desc>Interval [c,d], which includes all eigenvalues in it</desc>
                                <value>{ -14 40 }</value>
                            </param>
                            <param type=\"integer\">
                                <name>Max. # of iterations</name>
                                <cTag>maxi</cTag>
                                <desc>maxi</desc>
                                <value>40</value>
                            </param>
                            <param type=\"integer\">
                                <name>degree of Chebyshev polynomial to be used</name>
                                <cTag>ChebyDeg</cTag>
                                <desc>ChebyDeg</desc>
                                <value>240</value>
                            </param>
                            <param type=\"boolean\">
                                <name>calculate eigenvectors</name>
                                <cTag>CalEigVec</cTag>
                                <desc>calculate eigenvectors?</desc>
                                <value>true</value>
                            </param>
                        </group>
";
 my($tracemin_options2)=
"                    <group type=\"struct\">
                        <name>Tracemin</name>
                        <cTag>PhonTracemin</cTag>
                        <desc>Tracemin parameters</desc>
                        <param type=\"integer\">
                            <name>maximum iteration number</name>
                            <cTag>MaxIter</cTag>
                            <desc>maximum iteration number</desc>
                            <value>3600</value>
                        </param>
                        <param type=\"integer\">
                            <name>number of Arnoldi vectors</name>
                            <cTag>NumArnoldiVector</cTag>
                            <desc>How many times of the number of
                                eigenvalues
                                is the number of Arnoldi vectors?</desc>
                            <value>3</value>
                        </param>
                        <param type=\"real\">
                            <name>reference eigenvalue</name>
                            <cTag>RefEigVal</cTag>
                            <desc>reference eigenvalue, eV</desc>
                            <value>0.0</value>
                        </param>
                        <param type=\"real\">
                            <name>tolerance</name>
                            <cTag>Tolerance</cTag>
                            <desc>tolerance</desc>
                            <value>1e-10</value>
                        </param>
                    </group>
";
my($OutforAtomicPos)=
"                    <param type=\"boolean\">
                        <name>Dump atomic positions before strain?</name>
                        <cTag>AtomPosBeforeStrain</cTag>
                        <desc>AtomPosBeforeStrain</desc>
                        <value>false</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Dump atomic positions at equil.? NOTE: need also PsiSqr set to true</name>
                        <cTag>AtomPosEquil</cTag>
                        <desc>AtomPosEquil</desc>
                        <value>false</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Dump atomic positions after strain?</name>
                        <cTag>AtomPosAfterStrain</cTag>
                        <desc>AtomPosAfterStrain</desc>
                        <value>false</value>
                    </param>
";
my($Bgauge)=
"                    <param type=\"string\">
                        <name>Gauge choice for Bfield</name>
                        <cTag>Bgauge</cTag>
                        <desc>Gauge for the magnetic field, 'Symmetric' gives A = (-By/2,Bx/2,0), 'Asymmetric_x' gives  A = (0,Bx,0) useful with Efield along x, 'Asymmetric_y' gives  A = (-By,0,0), useful with Efield along y</desc>
                        <value>NONE</value>
		        <option>NONE</option>
		        <option>Symmetric</option>
		        <option>Asymmetric_x</option>
		        <option>Asymmetric_y</option>
                    </param>
";
my($donor_opts)=
"                    <param type=\"integer\">
                        <name>Number of impurities?</name>
                        <cTag>ImpNumber</cTag>
                        <desc>Number of P+ donors in the system (max.3)</desc>
                        <value>0</value>
                    </param>
                    <param type=\"array_real\">
                        <name>Coords. of impurities. Needed x1,y1,z1, x2,y2,z2,...,xN,yN,zN with N=ImpNumber </name>
                        <cTag>ImpsXYZ</cTag>
                        <desc>ImpsXYZ</desc>
                        <value>{0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0}</value>
                    </param>
                    <param type=\"integer\">
                        <name>imp_model</name>
                        <cTag>imp_model</cTag>
                        <desc>Use 1 to include U0 only, 2 includes also Si:P parameters</desc>
                        <value>2</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Include nuclear interaction between impurities?</name>
                        <cTag>include_nuclear_term</cTag>
                        <desc>include_nuclear_term</desc>
                        <value>false</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Enable constant electric field?</name>
                        <cTag>EfieldON</cTag>
                        <desc>EfieldON</desc>
                        <value>false</value>
                    </param>
                    <param type=\"real\">
                        <name>Ex</name>
                        <cTag>Ex</cTag>
                        <desc>x component of electric field (eV)</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>Ey</name>
                        <cTag>Ey</cTag>
                        <desc>y component of electric field (eV)</desc>
                        <value>0</value>
                    </param>
                    <param type=\"real\">
                        <name>Ez</name>
                        <cTag>Ez</cTag>
                        <desc>z component of electric field (eV)</desc>
                        <value>0</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Dump donor file?</name>
                        <cTag>donor_file</cTag>
                        <desc>donor_file</desc>
                        <value>false</value>
                    </param>
                    <param type=\"boolean\">
                        <name>Dump Efield file?</name>
                        <cTag>electric_file</cTag>
                        <desc>electric_file</desc>
                        <value>false</value>
                    </param>
";
 my($dx_option)= 
"                        <param type=\"boolean\">
                            <name>StructureDX</name>
                            <cTag>StructureDX</cTag>
                            <desc>Dump out .dx file for eigenvectors, for volqd visualization</desc>
                            <value>false</value>
                        </param>
";
 my($New_Lancz_opts)= 
"                            <param type=\"integer\">
                                <name>conv_method</name>
                                <cTag>conv_method</cTag>
                                <desc>1 for eigenvalues and 2 for eigenvectors</desc>
                                <value>1</value>
                            </param>
                            <param type=\"integer\">
                                <name>conv_vec_method</name>
                                <cTag>conv_vec_method</cTag>
                                <desc>If conv_method=2: 1 for ALL components and 2 for selected comps. </desc>
                                <value>2</value>
                            </param>
";
 my($New_Lancz_opts2)=
"                           <param type=\"array_real\">
                                <name>Valence band energy range</name>
                                <cTag>Erange_vb</cTag>
                                <desc>Search range for valence band eigenvalues</desc>
                                <value>{ -0.6 0.2 }</value>
                            </param>
                            <param type=\"array_real\">
                                <name>Conduction band energy range</name>
                                <cTag>Erange_cb</cTag>
                                <desc>Search range for conduction band eigenvalues</desc>
                                <value>{ 1.4 2.0}</value>
                            </param>
";
 my($BL_param)=
"
                        <group
                            enabledif=\"/application/input/top/Execution parameters/Electronic structure calculation/Algorithms/Electronic eigensolver==BlockLanczos\" type=\"struct\">
                            <name>BLanczos</name>
                            <cTag>BLanczos</cTag>
                            <desc>block lanczos parameters</desc>
                            <param type=\"integer\">
                                <name>BL_size</name>
                                <cTag>BL_size</cTag>
                                <desc>number of vectors in a block for each iteration:gives deg.</desc>
                                <value>6</value>
                            </param>
                            <param type=\"integer\">
                                <name>Deflation</name>
                                <cTag>DeflationOn</cTag>
                                <desc>check for deflation at each iteration step</desc>
                                <value>0</value>
                            </param>
                            <param type=\"real\">
                                <name>DefTolerance</name>
                                <cTag>DefTolerance</cTag>
                                <desc>vector discarded if its norm is less than DefTolerance</desc>
                                <value>1e-07</value>
                            </param>
                            <param type=\"integer\">
                                <name>Max_no_success</name>
                                <cTag>Max_no_success</cTag>
                                <desc>stop if no eigenvalues found in Max_no_success iterations</desc>
                                <value>1000</value>
                            </param>
                            <param type=\"integer\">
                                <name>Eig_dump</name>
                                <cTag>Eig_dump</cTag>
                                <desc>directly dump out the eigenvalues</desc>
                                <value>0</value>
                            </param>
                            <param type=\"integer\">
                                <name>num_cb</name>
                                <cTag>num_cb</cTag>
                                <desc>if Eig_dump=true, number of CB eigenvalues to be dumped</desc>
                                <value>15</value>
                            </param>
                            <param type=\"integer\">
                                <name>num_vb</name>
                                <cTag>num_vb</cTag>
                                <desc>if Eig_dump=true, number of VB eigenvalues to be dumped</desc>
                                <value>0</value>
                            </param>
                            <param type=\"integer\">
                                <name>Dump_blockTridiagonal_Sparse_Matrix_Format</name>
                                <cTag>Dump_blockTridiagonal_Sparse_Matrix_Format</cTag>
                                <desc>Dump_blockTridiagonal_Sparse_Matrix_Format</desc>
                                <value>0</value>
                            </param>
                            <param type=\"integer\">
                                <name>Dump_blockTridiagonal_Band_Matrix_Format</name>
                                <cTag>Dump_blockTridiagonal_Band_Matrix_Format</cTag>
                                <desc>Dump_blockTridiagonal_Band_Matrix_Format</desc>
                                <value>0</value>
                            </param>
                            <param type=\"integer\">
                                <name>Read_blockTridiagonal_Band_Matrix_Format</name>
                                <cTag>Read_blockTridiagonal_Band_Matrix_Format</cTag>
                                <desc>Read_blockTridiagonal_Band_Matrix_Format</desc>
                                <value>0</value>
                            </param>
                            <param type=\"boolean\">
                                <name>ReadInTridiag</name>
                                <cTag>ReadInTridiag</cTag>
                                <desc>ReadInTridiag</desc>
                                <value>false</value>
                            </param>
                            <param type=\"real\">
                                <name>Convergence Tolerance</name>
                                <cTag>ConvTolerance</cTag>
                                <desc>NONE</desc>
                                <value>1e-07</value>
                            </param>
                            <param type=\"integer\">
                                <name>Maximum Iteration Number</name>
                                <cTag>MaxIter</cTag>
                                <desc>NONE</desc>
                                <value>2400</value>
                            </param>
                            <param type=\"integer\">
                                <name>conv_method</name>
                                <cTag>conv_method</cTag>
                                <desc>1 for eigenvalues and 2 for eigenvectors</desc>
                                <value>1</value>
                            </param>
                            <param type=\"integer\">
                                <name>conv_vec_method</name>
                                <cTag>conv_vec_method</cTag>
                                <desc>If conv_method=2: 1 for ALL components and 2 for selected comps. </desc>
                                <value>2</value>
                            </param>
                            <param type=\"integer\">
                                <name>Number of requested VB eigenvalues</name>
                                <cTag>NumEigReq_vb</cTag>
                                <desc>NONE</desc>
                                <value>1</value>
                            </param>
                           <param type=\"array_real\">
                                <name>Valence band energy range</name>
                                <cTag>Erange_vb</cTag>
                                <desc>Search range for valence band eigenvalues</desc>
                                <value>{ -1.0 0.0 }</value>
                            </param>
                            <param type=\"array_real\">
                                <name>Conduction band energy range</name>
                                <cTag>Erange_cb</cTag>
                                <desc>Search range for conduction band eigenvalues</desc>
                                <value>{ 0.9 1.6 }</value>
                            </param>
                            <param type=\"integer\">
                                <name>Number of requested CB eigenvalues</name>
                                <cTag>NumEigReq_cb</cTag>
                                <desc>NONE</desc>
                                <value>1</value>
                            </param>
                            <param type=\"integer\">
                                <name>Start Convergence Check</name>
                                <cTag>ConvCheckStartIter</cTag>
                                <desc>first iteration number to check
                                    for convergence</desc>
                                <value>200</value>
                            </param>
                            <param type=\"integer\">
                                <name>Convergence Check Interval</name>
                                <cTag>ConvCheckSkipRate</cTag>
                                <desc>convergence check skip rate</desc>
                                <value>50</value>
                            </param>
                            <param type=\"option\">
                                <name>Initial condition</name>
                                <cTag>Start</cTag>
                                <desc>inital guess</desc>
                                <value>CyclicRandom</value>
                                <option>CyclicRandom</option>
                                <option>Random</option>
                                <option>Const1</option>
                                <option>ReadFile</option>
                            </param>
                            <param type=\"string\">
                                <name>File for Initial Guess</name>
                                <cTag>StartFileRead</cTag>
                                <desc>NONE</desc>
                                <value>NULL</value>
                            </param>
                            <param type=\"integer\">
                                <name>Max number of iterations w/o eigval</name>
                                <cTag>ConvSingleValue</cTag>
                                <desc>maximum number of iterations
                                    without a
                                    converged eigenvalue</desc>
                                <value>1000</value>
                            </param>
                            <param type=\"array_real\">
                                <name>CyclicString</name>
                                <cTag>CyclicString</cTag>
                                <desc>CyclicString</desc>
                                <value>{ 1.2 0.8 -1.4 0.1 -0.9 -0.2 0.6 }</value>
                            </param>
                        </group>
";
my($tiltbox) = 
"                            <option>TiltedBox</option>
";
my($tiltbox2) = 
"                                    <option>TiltedBox</option>
";
my($newCylinderOpt)=
"                            <option>Cylinder_x</option>
";
my($newCylinderOpt2)=
"                                    <option>Cylinder_x</option>
";
my($DXOption) = 
"                        <param type=\"boolean\">
                            <name>Visualization_3D</name>
                            <cTag>Visualization_3D</cTag>
                            <desc>write output files of dx-format for 3D visualization</desc>
                            <value>false</value>
                        </param>
";
my($TensionOutpt)=
"                        <param type=\"boolean\">
                            <name>output strain from the electronic domain only?</name>
                            <cTag>StrainElOnly</cTag>
                            <desc>StrainElOnly</desc>
                            <value>false</value>
                        </param>
                        <param type=\"boolean\">
                            <name>output file for tension</name>
                            <cTag>TensionOutput</cTag>
                            <desc>write output files of six tension components for visualization</desc>
                            <value>false</value>
                        </param>
";

 update_version_ID($fil,$fil_out);

  @strList = ("ShapeName", "Ellipsoid");
  insert_txt($fil_out, \@strList, $tiltbox2);

  @strList = ("ShapeName_1", "Ellipsoid");
  insert_txt($fil_out, \@strList, $tiltbox);

  @strList = ("ShapeName", "Cylinder");
  insert_txt($fil_out, \@strList, $newCylinderOpt2);

  @strList = ("ShapeName_1", "Cylinder");
  insert_txt($fil_out, \@strList, $newCylinderOpt);

  @strList = ("a_z", "a_lattice_z", "/param>");
  insert_txt($fil_out, \@strList, $tilted_hkl);

  @strList = ("material_13<", "mb_strain", "</group>");
  insert_txt($fil_out, \@strList, $new_mats);
  
  @start = ("material_1<", "GaAs", "mb_pv","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1<", "GaAs", "mb_pv","</desc>");
  insert_txt($fil_out, \@strList, $new_GaAs_param);

  @start = ("material_1<", "GaAs", "mb_strain","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1<", "GaAs", "mb_strain","</desc>");
  insert_txt($fil_out, \@strList, $new_GaAs_mb_strain);

  @start = ("material_2<", "InAs", "mb_pv","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2<", "InAs", "mb_pv","</desc>");
  insert_txt($fil_out, \@strList, $new_InAs_param);

  @start = ("material_2<", "InAs", "mb_strain","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2<", "InAs", "mb_strain","</desc>");
  insert_txt($fil_out, \@strList, $new_InAs_mb_strain);

  @start = ("material_4<", "Si", "mb_pv","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_4<","Si", "mb_pv","</desc>");
  insert_txt($fil_out, \@strList, $new_Si_param);

  @strList = ("ElAlg","<param","Readstates","</option>");
  insert_txt($fil_out, \@strList, $new_solver_options);
  
  @strList = ("enabledif", "Parpack","</param>", "</group>");
  insert_txt($fil_out, \@strList, $tracemin_options);

  @strList = ("PhonParpack", "Tolerance","</param>", "</group>");
  insert_txt($fil_out, \@strList, $tracemin_options2);

  @strList = ("enabledif", "Lanczos","MaxIter","</param>");
  insert_txt($fil_out, \@strList, $New_Lancz_opts);

  @strList = ("enabledif", "Lanczos","NumEigReq_cb","</param>");
  insert_txt($fil_out, \@strList, $New_Lancz_opts2);

  @strList = ("enabledif", "Lanczos","</param>", "</group>");
  insert_txt($fil_out, \@strList, $BL_param);

  @start = ("ResFind", "method", "</param>", "<param");
  @end = ("</param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("ResFind", "method", "</param>", "<param");
  @end = ("</param>");
  delete_txt($fil_out, \@start, \@end);

  @start = ("Bz", "field", "</param>", "<param");
  @end = ("</param>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("Bz","field","</param>");
  insert_txt($fil_out, \@strList, $donor_opts);

  @strList = ("Bz","field","</param>");
  insert_txt($fil_out, \@strList, $Bgauge);

  @strList = ("PsiSqr_cell", "modulus", "</param>");
  insert_txt($fil_out, \@strList, $DXOption);
  #@strList = ("StructurePdb", "StructurePDB", "</param>");
  #insert_txt($fil_out, \@strList, $dx_option);

  @start = ("AtomInfo", "NbrIndx", "</param>", "<param");
  @end = ("</param>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("AtomInfo", "NbrIndx", "</param>");
  insert_txt($fil_out, \@strList, $OutforAtomicPos);

  @strList = ("Bondlength", "value", "</param>");
  insert_txt($fil_out, \@strList, $TensionOutpt);
}

sub update_version_0_18 {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);
  my($newer_GaAs_param) =
"                        <value>{-5.50042 4.15107 -0.24119 6.70776
                            19.71059 22.66352 13.03169 12.74846 0.17234
                            0.02179 -1.64508 -3.704550 -2.20777 -1.31491
                            2.66493 2.96032 1.97650 1.02755 -2.609400
                            -2.32059 -0.62820 -0.13324 4.15080 -1.42744
                            -1.87428 -1.88964 2.52926 2.54913 -1.26996
                            2.50536 -0.85174 0.56532 1 1 1 }</value>
";
  my($newer_GaAs_mb_strain) =
"                        <value>{ 0.00000 0.21266 2.06001 1.38498 2.68497
                            1.31405 1.89889 1.39930 1.81235 2.37964
                            1.0000 1.72443 1.97253 1.89672 1.78540 27.0
                            27.0                             0.58696
                            0.48609 0.88921 0.77095 0.75979
                            1.45891 0.81079 1.21202 1.07015 0.58053
                            1.03256 1.32385 2.00000 1.61350 1.50000
                            1.26262 }</value>
";
 my($newer_InAs_param) =
"                        <value>{ -5.50042 4.15107 -0.581930 6.97163
                            19.71059 19.94138 13.03169 13.30709 0.17234
                            0.131200 -1.69435 -4.21045 -2.42674 -1.159870
                            2.598230 2.809360 2.06766 0.93734 -2.26837
                            -2.293090 -0.899370 -0.488990 4.310640 -1.288950
                            -1.731410 -1.978420 2.188860 2.456020 -1.584610
                            2.71793 -0.505090 0.60583 1 1 1 }</value>
";
my($newer_InAs_mb_strain) =
"                        <value>{ 0.060800 0.000810 1.92494 1.570030 2.06151
                            1.602470 1.765660 1.949370 2.383820 2.455600
                            1.0000 2.32291 1.615890 2.3296 2.02387 27.0
                            27.0                             1.258286
                            2.481447 1.086223 4.557774 4.367575
                            7.02966 3.298598 7.029496 0.0000 0.187036
                            1.195042 1.769483 4.624483 0.0000 0.0000
                            0.246999 }</value>
";

 update_version_ID($fil,$fil_out);

  @start = ("material_1<", "GaAs", "mb_pv","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1<", "GaAs", "mb_pv","</desc>");
  insert_txt($fil_out, \@strList, $newer_GaAs_param);

  @start = ("material_1<", "GaAs", "mb_strain","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_1<", "GaAs", "mb_strain","</desc>");
  insert_txt($fil_out, \@strList, $newer_GaAs_mb_strain);

  @start = ("material_2<", "InAs", "mb_pv","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2<", "InAs", "mb_pv","</desc>");
  insert_txt($fil_out, \@strList, $newer_InAs_param);

  @start = ("material_2<", "InAs", "mb_strain","<value>");
  @end = ("</value>");
  delete_txt($fil_out, \@start, \@end);

  @strList = ("material_2<", "InAs", "mb_strain","</desc>");
  insert_txt($fil_out, \@strList, $newer_InAs_mb_strain);
}
#-------------------------------------------------------------------------------
# probably don't need to touch these functions
#-------------------------------------------------------------------------------

sub get_next_version_ID {
  my($version) = @_;
  my($version_next) = "";

  $i_prev="";
  foreach $i (@versionList) {
    if ($i_prev eq $version) {
      return $i;
    }
    $i_prev = $i;
  }

  return "";
}


sub get_current_version_ID {
  my($fil) = @_;
  my($version) = "";
  my $line;

  open(FP, $fil) || die("Could not open file:  '", $fil,"'");
  while (defined($line = <FP>) && $version eq "") {
    ($y) = $line =~ m/(version of XML file)/;
    if ($y ne "") {
      $line = <FP>;
      ($version) = $line =~ m/>(.+)</;
    }
  }
  close(FP);

  $version;
}


sub update_version_ID {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);
  local $in_and_out_are_same = 0;

  my $version=get_current_version_ID($fil);
  my $version_next=get_next_version_ID($version);

  $in_and_out_are_same=1 if ($fil eq $fil_out);

  if ($in_and_out_are_same==1) {
    $fil_out="tmp_lkjxvizpwr";
  }

  open(FP, $fil) || die("Could not open file:  '", $fil,"'");
  open(FP_OUT, ">".$fil_out);
  while (defined($line = <FP>)) {
    my $y;

    ($y) = $line =~ m/(version of XML file)/;
    if ($y ne "") {
      print FP_OUT $line;
      $line = <FP>;
      $line =~ s?$version?$version_next?;
    }

    print FP_OUT $line;

  }
  close(FP_OUT);
  close(FP);

  if ($in_and_out_are_same==1) {
    system("mv -f tmp_lkjxvizpwr " . $fil);
  }
}

#delete text between the start and end lines (including the boundary lines)
sub delete_txt {
  my($end) = pop(@_);
  my($start) = pop(@_);
  my($fil) = pop(@_);

  my($output_string)="";
  my($deleting) = 0;

  open(FP, $fil) || die("Could not open file:  '", $fil,"'");
  while (defined($line = <FP>)) {
    my($s) = $line =~ m/(@$start[0])/;
    my($e) = $line =~ m/(@$end[0])/;

    if ($s ne "" && $deleting==0) {

      shift(@$start);

      if (@$start[0] eq "") {
	$deleting=1;
      }
    }

    if ($deleting==0) {
      $output_string = $output_string . $line;
    }

    if ($e ne "" && $deleting==1) {
      shift(@$end);

      if (@$end[0] eq "") {
	$deleting=0;
      }
    }

  }
  close(FP);

  open(FP_OUT, ">".$fil);
  print FP_OUT $output_string;
  close(FP_OUT);
}

#insert text below the line pointer
sub insert_txt {
  my($str) = pop(@_);
  my($start) = pop(@_);
  my($fil) = pop(@_);

  my($output_string)="";
  my($no_insertion_yet)=1;

  open(FP, $fil) || die("Could not open file:  '", $fil,"'");
  while (defined($line = <FP>)) {
    my($s) = $line =~ m/(@$start[0])/;

    $output_string = $output_string . $line;

    if ($s ne "" && $no_insertion_yet) {
      shift(@$start);

      if (@$start[0] eq "") {
	$output_string = $output_string . $str;
	$no_insertion_yet = 0;
      }
    }
  }
  close(FP);

  open(FP_OUT, ">".$fil);
  print FP_OUT $output_string;
  close(FP_OUT);
}


sub update_version {
  my($fil_out) = pop(@_);
  my($fil) = pop(@_);
  local $fil_tmp_in = $fil;
  local $fil_tmp_out;
  my $version, $version_next, $cmd;

  $version=get_current_version_ID($fil);
  while ($version ne "") {
    $version_next = get_next_version_ID($version);
    if ($version_next ne "") {
      $fil_tmp_out = $fil . "_version=" . $version_next;
      print "Updating:  " . $fil_tmp_in .  "\n"; #"  ==> " . $fil_tmp_out . "\n";

      ($str = $version) =~ s/\./_/;
#line commented by Marta, so the input and output file is THE SAME
      eval("update_version_" . $str . '($fil_tmp_in,$fil_tmp_in);');
    #  eval("update_version_" . $str . '($fil_tmp_in,$fil_tmp_out);');
    }

    $fil_tmp_in = $fil_tmp_out;
    $version = $version_next;
  }

  if ($fil_out ne "") {
    system("cp -f " . $fil_tmp_out . " " . $fil);
  }

}


#-------------------------------------------------------------------------------
#  main
#-------------------------------------------------------------------------------

$fil_in = $ARGV[0];

if ($#ARGV==0) {
  $fil_out = "";
}
elsif ($#ARGV==1) {
  $fil_out = $ARGV[1];
}
else {
  die("Incorrect number of arguments.");
}

update_version($fil_in,$fil_out);
