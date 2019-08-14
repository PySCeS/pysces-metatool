pysces-metatool
===============

Metatool and bindings for use with PySCeS - adds elementary mode analysis. The MetaTool source code has been slightly altered to allow building under Linux and for use with PySCeS

Brett G. Olivier, Amsterdam 2017.

PySCeS
------

The official code repository (http://github.com/PySCeS) of The Python Simulator for Cellular Systems: PySCeS project.

[![Research software impact](http://depsy.org/api/package/pypi/PySCeS/badge.svg)](http://depsy.org/package/python/PySCeS)

Copyright (c) 2004 - 2017, Brett G. Olivier, Johann M. Rohwer and Jan-Hendrik S. Hofmeyr
All rights reserved.

Metatool
--------

METATOOL is a C program developed from 1998 to 2000 by Thomas Pfeiffer
(Berlin) in cooperation with Stefan Schuster and Ferdinand
Moldenhauer (Berlin) and Juan Carlos Nuno (Madrid).
It serves to derive conclusions about the pathway
structure of metabolic networks from the stoichiometric reaction
equations and information about reversibility  and irreversibility of
enzymes.
It should preferably be compiled with the GNU compiler.
For DOS and Win32 console applications, comment out the two lines
#include<conio.h> and #include<malloc.h>.
The program requires the two names of the input and output files at
the command line.
To explain the format of the input file, we give an example file
(Example.dat), which codifies a reaction scheme comprising the
tricarboxylic acid cycle, glyoxylate shunt and adjacent reactions of
amino acid synthesis in E. coli (cf. Ref. 1).

-ENZREV
Eno Acn SucCD Sdh Fum Mdh AspC Gdh IlvEAvtA

-ENZIRREV
Pyk AceEF GltA Icd SucAB Icl Mas AspCon AspA Pck Ppc Pps GluCon
AlaCon SucCoACon

-METINT
Ala Asp Glu Gly Mal Fum Succ SucCoA OG IsoCit Cit OAA AcCoA CoA
Pyr PEP

-METEXT
Sucex Alaex Gluex ADP ATP AMP NH3 Aspex FADH2 FAD GTP
GDP NADPH NADP NADH CO2 NAD PG

-CAT
Eno : PG = PEP .
Pyk : PEP + ADP = Pyr + ATP .
AceEF : Pyr + NAD + CoA = AcCoA + CO2 + NADH .
GltA : OAA + AcCoA = Cit + CoA .
Acn : Cit = IsoCit .
Icd : IsoCit + NADP = OG + CO2 + NADPH .
SucAB : OG + NAD + CoA = SucCoA + CO2 + NADH .
SucCD : SucCoA + ADP = Succ + ATP + CoA .
Sdh : Succ + FAD = Fum + FADH2 .
Fum : Fum = Mal .
Mdh : Mal + NAD = OAA + NADH .
Icl : IsoCit = Succ + Gly .
Mas : Gly + AcCoA = Mal + CoA .
AspC : OAA + Glu = Asp + OG .
AspCon : Asp = Aspex .
AspA : Asp = Fum + NH3 .
Gdh : OG + NH3 + NADPH = Glu + NADP .
Pck : OAA + ATP = PEP + ADP + CO2 .
Ppc : PEP + CO2 = OAA .
Pps : Pyr + ATP = PEP + AMP .
GluCon : Glu = Gluex .
IlvEAvtA : Pyr + Glu = Ala + OG .
AlaCon : Ala = Alaex .
SucCoACon : SucCoA = Sucex + CoA .

[Explanation:
-ENZREV, -ENZIRREV
After the key words -ENZREV and -ENZIRREV, names or abbreviations of
the reversible and irreversible enzymes, respectively, have to be
written.
-METINT, -METEXT
After the key word -METINT, names or abbreviations of the internal
metabolites have to be written. These are the substances which have to
fulfil a steady-state condition (production = consumption).
After the key word -METEXT names or abbreviations of the external
metabolites have to be written. External metabolites (sources and
sinks) need not be balanced in the scheme under consideration.
The order of these four fields is important. All internal and external
metabolites must have an underscore or a letter (no number) as the first
character and must not include a white space.
-CAT
After the key word -CAT, the reaction equations are listed in any
order. The raction name is written first just as after the key words
-ENZREV and -ENZREV. The reaction name is followed by a white space
(space or tab), a colon and a white space. Then the
stoichiometric reaction equation follows. Stoichiometric coefficients are
integers separated by a white space from the metabolites. After the
metabolites, a white space and a plus or a white space and an
equal sign follow. The end of each reaction is formed by a white space and a
full stop. Metabolites are written in the same way as after the key words
-METINT and -METEXT.
The order of metabolites in the reaction equations makes no
difference. However, the sides of the reaction equations are exchangeable
only
in the reversible reactions. The metabolites that are formed by the
irreversible reactions have to be written on the right side of the reaction
equations.]

----------------------------------------------------------------------
The program writes the results in an output file. For our example,
this file reads as follows:

METATOOL OUTPUT Version 2.0 [your path to]\meta2.exe

INPUT FILE: Example.dat

INTERNAL METABOLITES: 19
REACTIONS: 24

STOICHIOMETRIC MATRIX:

 matrix dimension 16 x 24
 -1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -1 -1  0  0  0  0  0  0
  0  0  0  0  0  0 -1  1 -1  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0
  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0
  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
  0  0 -1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0 -1
  0  0  0  0  0  0  1 -1  1  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0
  0  1  0  0  0  0  0  0  0  0  0  0 -1  0 -1  0  0  0  0  0  0  0  0  0
  0 -1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  1 -1  0  0  0  0 -1  0  0  0  0  0  0 -1  1  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  1 -1  0  0  0 -1  0  0  0  0  0  0  0  0
  0  0  1  0  0  0  0  0  0  0 -1  1  0 -1  0  1  0  0  0  0  0  0  0  1
  0  0  0  0  0  0  0  0 -1  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0
  1  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  1 -1  1  0  0  0
following line indicates reversible (0) and irreversible reactions (1)
  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1

[Explanation: The program gives the numbers of internal metabolites and
reactions. It also parses the reaction equations and translates them into
a stoichiometric matrix. This matrix includes the stoichiometric
coefficients (molecularities) of the internal metabolites in all the
reaction
equations, with the rows corresponding to internal metabolites and the
columns corresponding to reactions.
The line following the stoichiometric matrix indicates the reversible
and irreversible reactions in the same order as after the key words
-ENZREV and -ENZIRREV.]


KERNEL

matrix dimension 9 x 24
 -1 -1 -1 -1 -1 -1  0  0  0 -1 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0  0
 -2 -1  0 -1 -1 -2 -1 -1  0 -2 -2 -1  0  0 -1 -1 -1  0  0  0  0  0  0  0
  0  0  0  0 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0
 -1 -1  0 -1 -1 -2  0  0  0 -2 -2 -1  0  0 -1 -1  0  0 -1  0  0  0  0  0
  1  1  0  1  1  2  0  0  0  2  2  1  0  0  1  1  0  0  0 -1  0  0  0  0
  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  0  0  0
 -3 -2  0 -1 -1 -2  0 -1  0 -3 -3 -2 -1  0 -1 -1  0  0  0  0  0 -1  0  0
  1  0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  0
  2  1 -1  0  0  1  0  0  0  2  2  1  0  0  1  1  0  0  0  0  0  0  0  1

[Explanation: The kernel or nullspace is the subspace of all flux
vectors V that satisfy the equation Stoich. matrix times V = 0
(see Ref. 2). The rows of the above matrix span this subspace.]

enzymes

 1:      -Eno -Acn -SucCD -Sdh -Fum -Mdh -Pyk -AceEF -GltA -Icd -SucAB
irreversible
 2:      (-2 Eno) -Acn -Sdh -Fum (-2 Mdh) -AspC -Gdh (-2 Pyk) (-2 AceEF)
-GltA -Icl
         -Mas -AspCon irreversible
 3:      -Fum -Mdh -AspC -Gdh -AspA irreversible
 4:      -Eno -Acn -Sdh -Fum (-2 Mdh) (-2 Pyk) (-2 AceEF) -GltA -Icl -Mas
-Pck
         irreversible
 5:      Eno Acn Sdh Fum (2 Mdh) (2 Pyk) (2 AceEF) GltA Icl Mas -Ppc
irreversible
 6:      Pyk Pps irreversible
 7:      (-3 Eno) (-2 Acn) -Sdh -Fum (-2 Mdh) -Gdh (-3 Pyk) (-3 AceEF) (-2
GltA) -Icd
         -Icl -Mas -GluCon irreversible
 8:      Eno Gdh IlvEAvtA Pyk AlaCon irreversible
 9:      (2 Eno) Acn -SucCD Mdh (2 Pyk) (2 AceEF) GltA Icl Mas SucCoACon
irreversible

[Explanation: This list contains the enzymes that correspond to the
rows of the kernel matrix. The coefficients indicate relative fluxes
carried by the enzymes. A minus sign before an enzyme name stands for -1.
The following nine lines contain the
sum of metabolites which are involved in these enzyme reactions. E.g.
in reaction 6, Pyk and Pps catalyse PEP + ADP = Pyr + ATP and Pyr +
ATP = PEP + AMP, respectively, which gives, as the overall reaction:
ADP = AMP]

 overall reaction

 1:     2 ATP + FADH2 + NADPH + 3 NADH + 3 CO2 = 2 ADP + FAD + NADP + 3 NAD
+ PG
 2:     2 ATP + Aspex + FADH2 + NADP + 4 NADH + 2 CO2 = 2 ADP + NH3 + FAD +
NADPH + 4
        NAD + 2 PG
 3:     NADP + NADH = NADPH + NAD
 4:     ATP + FADH2 + 4 NADH + 3 CO2 = ADP + FAD + 4 NAD + PG
 5:     2 ADP + FAD + 4 NAD + PG = 2 ATP + FADH2 + 4 NADH + 3 CO2
 6:     ADP = AMP
 7:     Gluex + 3 ATP + FADH2 + 5 NADH + 4 CO2 = 3 ADP + NH3 + FAD + 5 NAD +
3 PG
 8:     ADP + NH3 + NADPH + PG = Alaex + ATP + NADP
 9:     ADP + 3 NAD + 2 PG = Sucex + ATP + 3 NADH + 2 CO2

SUBSETS of reactions (21 rows)

 matrix dimension 21 x 24
 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0
  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1

[Explanation: Enzyme subsets are sets of enzymes that always operate
together in fixed flux ratios. For example, if aconitase (Acn) is operative,
then also citrate synthase (GltA) is operative. This information can be
written in the form of a matrix (see above). For example, the second
row contains ones at positions 2 and 12, which correspond to Acn and
GltA. Below, this information is given in more detailed form, together
with the overall reactions of the subsets.]

enzymes

 1:      -Eno reversible
 2:      Acn GltA irreversible
 3:      -SucCD reversible
 4:      -Sdh reversible
 5:      -Fum reversible
 6:      -Mdh reversible
 7:      -AspC reversible
 8:      -Gdh reversible
 9:      IlvEAvtA AlaCon irreversible
 10:     Pyk irreversible
 11:     AceEF irreversible
 12:     Icd irreversible
 13:     SucAB irreversible
 14:     Icl Mas irreversible
 15:     AspCon irreversible
 16:     AspA irreversible
 17:     Pck irreversible
 18:     Ppc irreversible
 19:     Pps irreversible
 20:     GluCon irreversible
 21:     SucCoACon irreversible

 overall reaction

 1:     PEP = PG
 2:     OAA + AcCoA = IsoCit + CoA
 3:     Succ + CoA + ATP = SucCoA + ADP
 4:     Fum + FADH2 = Succ + FAD
 5:     Mal = Fum
 6:     OAA + NADH = Mal + NAD
 7:     Asp + OG = Glu + OAA
 8:     Glu + NADP = OG + NH3 + NADPH
 9:     Glu + Pyr = OG + Alaex
 10:    PEP + ADP = Pyr + ATP
 11:    CoA + Pyr + NAD = AcCoA + NADH + CO2
 12:    IsoCit + NADP = OG + NADPH + CO2
 13:    OG + CoA + NAD = SucCoA + NADH + CO2
 14:    IsoCit + AcCoA = Mal + Succ + CoA
 15:    Asp = Aspex
 16:    Asp = Fum + NH3
 17:    OAA + ATP = PEP + ADP + CO2
 18:    PEP + CO2 = OAA
 19:    Pyr + ATP = PEP + AMP
 20:    Glu = Gluex
 21:    SucCoA = CoA + Sucex

[Explanation: Enzymes belonging to the same subset can be
lumped. This gives rise to the following reduced reaction system.]


REDUCED SYSTEM with 13 branch point metabolites in 21 reactions
(columns)

 matrix dimension 13 x 21
  0  0  0  0  0  0 -1  0  0  0  0  0  0  0 -1 -1  0  0  0  0  0
  0  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0 -1  0
  0  0  0  0 -1  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
  0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0
  0  0 -1  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
  0  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1
  0  0  0  0  0  0 -1  1  1  0  0  1 -1  0  0  0  0  0  0  0  0
  0  1  0  0  0  0  0  0  0  0  0 -1  0 -1  0  0  0  0  0  0  0
  0 -1  0  0  0 -1  1  0  0  0  0  0  0  0  0  0 -1  1  0  0  0
  0 -1  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0  0  0  0  0
  0  1 -1  0  0  0  0  0  0  0 -1  0 -1  1  0  0  0  0  0  0  1
  0  0  0  0  0  0  0  0 -1  1 -1  0  0  0  0  0  0  0 -1  0  0
 -1  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  1 -1  1  0  0
following line indicates reversible (0) and irreversible reactions (1)
  0  1  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1

[Explanation: "The simplified system is a kind of skeleton model of
the original system. It contains only metabolites at branch points.
Skeleton models are often used in metabolic modeling to reduce the
number of variables" (see Refs. 5-7).

CONVEX BASIS

 matrix dimension 12 x 21
  0  0  0  0 -1 -1 -1 -1  0  0  0  0  0  0  0  1  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0
  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0
 -1  0  0  0  0  0 -1 -1  0  0  0  0  0  0  1  0  0  1  0  0  0
 -1  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1
 -1  0  0  0  0  0  0 -1  1  1  0  0  0  0  0  0  0  0  0  0  0
 -1  1 -1 -1 -1 -1  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0
 -2  1  0  0  0  0  0 -1  0  1  1  1  0  0  0  0  0  1  0  1  0
 -1  1  0 -1 -1 -2  0  0  0  2  2  0  0  1  0  0  1  0  0  0  0
 -2  1  0 -1 -1 -2 -1 -1  0  2  2  0  0  1  1  0  0  0  0  0  0
 -2  1  1  0  0 -1  0  0  0  2  2  0  0  1  0  0  0  0  0  0  1
 -3  2  0 -1 -1 -2  0 -1  0  3  3  1  0  1  0  0  0  0  0  1  0

 enzymes

 1:      Fum Mdh AspC Gdh AspA irreversible
 2:      Pck Ppc irreversible
 3:      Pyk Pps irreversible
 4:      Eno AspC Gdh AspCon Ppc irreversible
 5:      Eno -SucCD -Sdh -Fum -Mdh Ppc SucCoACon irreversible
 6:      Eno Gdh IlvEAvtA Pyk AlaCon irreversible
 7:      Eno Acn SucCD Sdh Fum Mdh Pyk AceEF GltA Icd SucAB irreversible
 8:      (2 Eno) Acn Gdh Pyk AceEF GltA Icd Ppc GluCon irreversible
 9:      Eno Acn Sdh Fum (2 Mdh) (2 Pyk) (2 AceEF) GltA Icl Mas Pck
irreversible
 10:     (2 Eno) Acn Sdh Fum (2 Mdh) AspC Gdh (2 Pyk) (2 AceEF) GltA Icl Mas
AspCon
         irreversible
 11:     (2 Eno) Acn -SucCD Mdh (2 Pyk) (2 AceEF) GltA Icl Mas SucCoACon
irreversible
 12:     (3 Eno) (2 Acn) Sdh Fum (2 Mdh) Gdh (3 Pyk) (3 AceEF) (2 GltA) Icd
Icl Mas
         GluCon irreversible

 overall reaction

 1:     NADPH + NAD = NADP + NADH
 2:     ATP = ADP
 3:     ADP = AMP
 4:     NH3 + NADPH + CO2 + PG = Aspex + NADP
 5:     ATP + FADH2 + NADH + CO2 + PG = Sucex + ADP + FAD + NAD
 6:     ADP + NH3 + NADPH + PG = Alaex + ATP + NADP
 7:     2 ADP + FAD + NADP + 3 NAD + PG = 2 ATP + FADH2 + NADPH + 3 NADH + 3
CO2
 8:     ADP + NH3 + NAD + 2 PG = Gluex + ATP + NADH + CO2
 9:     ADP + FAD + 4 NAD + PG = ATP + FADH2 + 4 NADH + 3 CO2
 10:    2 ADP + NH3 + FAD + NADPH + 4 NAD + 2 PG = 2 ATP + Aspex + FADH2 +
NADP + 4
        NADH + 2 CO2
 11:    ADP + 3 NAD + 2 PG = Sucex + ATP + 3 NADH + 2 CO2
 12:    3 ADP + NH3 + FAD + 5 NAD + 3 PG = Gluex + 3 ATP + FADH2 + 5 NADH +
4 CO2

[Explanation: The convex basis is the minimum number of elementary
modes to reconstruct the whole reaction system. Any admissible flux
distribution in the system (i.e. any distribution that is compatible
with the steady-state condition and the directionality of the irreversible
reactions) can be written as a non-negative linear
combination of the vectors forming the convex basis (Ref. 5). These vectors
form
the rows of the above matrix. These rows are then translated into lists
of enzymes in the same way as have been translated above the rows of the
null-space matrix. A basis vector is reversible if its negative is an
admissible flux distribution as well, otherwise it is irreversible.]

CONSERVATON RELATIONS

 matrix dimension 1 x 16
  0  0  0  0  0  0  0 -1  0  0  0  0 -1 -1  0  0

[Explanation - Conservation relations indicate that linear combinations
(e.g. the sum) of several internal
metabolites are constant. The metabolites are in the same order as
after the keyword -METINT. The above row means that
SucCoA + AcCoA + CoA = const. (The minus sign in the above row is irrelevant
because we can multiply the equation by -1).]

ELEMENTARY MODES

 matrix dimension 16 x 21
  0  0  0  0 -1 -1 -1 -1  0  0  0  0  0  0  0  1  0  0  0  0  0
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0
  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0
 -1  0  0  0  0  0 -1 -1  0  0  0  0  0  0  1  0  0  1  0  0  0
 -1  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1
 -1  0  1  1  0  0 -1 -1  0  0  0  0  0  0  0  1  0  1  0  0  1
 -1  0  0  0  0  0  0 -1  1  1  0  0  0  0  0  0  0  0  0  0  0
 -2  1  1  0  0 -1  0  0  0  2  2  0  0  1  0  0  0  0  0  0  1
 -1  1 -1 -1 -1 -1  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0
 -3  1  2  1  1  0  0  0  0  2  2  0  0  1  0  0  0  1  0  0  2
 -2  1  0  0  0  0  0 -1  0  1  1  1  0  0  0  0  0  1  0  1  0
 -2  1  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  1  0  0  1
 -1  1  0 -1 -1 -2  0  0  0  2  2  0  0  1  0  0  1  0  0  0  0
 -2  1  0 -1 -1 -2 -1 -1  0  2  2  0  0  1  1  0  0  0  0  0  0
 -3  2  0 -1 -1 -2  0 -1  0  3  3  1  0  1  0  0  0  0  0  1  0
 -3  2  0 -1 -1 -2  0  0  0  3  3  1  1  1  0  0  0  0  0  0  1

[Explanation: The choice of the basis vectors of the kernel (or
nullspace) is not unique. Therefore, it was proposed (Refs. 1,3-5,7) to
take a complete set of the simplest basis vectors compatible with the
directionality of the irreversible reactions. These are called
elementary modes. There may be more of them then actually needed to
span the admissible region in flux space, but they have the favourable
property to be uniquely determined (up to scalar multiples). These
modes can be brought in relation with the biochemical pathways in the
system. The rows of the elementary modes matrix give the elementary
modes for our example system.]

[Explanation: Below goes the verbal listing of the elementary modes and
of the overall reactions in terms of the external metabolites:]

enzymes

 1:      Fum Mdh AspC Gdh AspA irreversible
 2:      Pck Ppc irreversible
 3:      Pyk Pps irreversible
 4:      Eno AspC Gdh AspCon Ppc irreversible
 5:      Eno -SucCD -Sdh -Fum -Mdh Ppc SucCoACon irreversible
 6:      Eno -SucCD -Sdh AspC Gdh AspA Ppc SucCoACon irreversible
 7:      Eno Gdh IlvEAvtA Pyk AlaCon irreversible
 8:      (2 Eno) Acn -SucCD Mdh (2 Pyk) (2 AceEF) GltA Icl Mas SucCoACon
irreversible
 9:      Eno Acn SucCD Sdh Fum Mdh Pyk AceEF GltA Icd SucAB irreversible
 10:     (3 Eno) Acn (-2 SucCD) -Sdh -Fum (2 Pyk) (2 AceEF) GltA Icl Mas Ppc
         (2 SucCoACon) irreversible
 11:     (2 Eno) Acn Gdh Pyk AceEF GltA Icd Ppc GluCon irreversible
 12:     (2 Eno) Acn Pyk AceEF GltA Icd SucAB Ppc SucCoACon irreversible
 13:     Eno Acn Sdh Fum (2 Mdh) (2 Pyk) (2 AceEF) GltA Icl Mas Pck
irreversible
 14:     (2 Eno) Acn Sdh Fum (2 Mdh) AspC Gdh (2 Pyk) (2 AceEF) GltA Icl Mas
AspCon
         irreversible
 15:     (3 Eno) (2 Acn) Sdh Fum (2 Mdh) Gdh (3 Pyk) (3 AceEF) (2 GltA) Icd
Icl Mas
         GluCon irreversible
 16:     (3 Eno) (2 Acn) Sdh Fum (2 Mdh) (3 Pyk) (3 AceEF) (2 GltA) Icd
SucAB Icl Mas
         SucCoACon irreversible

 overall reaction

 1:     NADPH + NAD = NADP + NADH
 2:     ATP = ADP
 3:     ADP = AMP
 4:     NH3 + NADPH + CO2 + PG = Aspex + NADP
 5:     ATP + FADH2 + NADH + CO2 + PG = Sucex + ADP + FAD + NAD
 6:     ATP + FADH2 + NADPH + CO2 + PG = Sucex + ADP + FAD + NADP
 7:     ADP + NH3 + NADPH + PG = Alaex + ATP + NADP
 8:     ADP + 3 NAD + 2 PG = Sucex + ATP + 3 NADH + 2 CO2
 9:     2 ADP + FAD + NADP + 3 NAD + PG = 2 ATP + FADH2 + NADPH + 3 NADH + 3
CO2
 10:    FADH2 + 2 NAD + 3 PG = 2 Sucex + FAD + 2 NADH + CO2
 11:    ADP + NH3 + NAD + 2 PG = Gluex + ATP + NADH + CO2
 12:    ADP + NADP + 2 NAD + 2 PG = Sucex + ATP + NADPH + 2 NADH + 2 CO2
 13:    ADP + FAD + 4 NAD + PG = ATP + FADH2 + 4 NADH + 3 CO2
 14:    2 ADP + NH3 + FAD + NADPH + 4 NAD + 2 PG = 2 ATP + Aspex + FADH2 +
NADP + 4
        NADH + 2 CO2
 15:    3 ADP + NH3 + FAD + 5 NAD + 3 PG = Gluex + 3 ATP + FADH2 + 5 NADH +
4 CO2
 16:    3 ADP + FAD + NADP + 6 NAD + 3 PG = Sucex + 3 ATP + FADH2 + NADPH +
6 NADH
        + 5 CO2


References

1. Schuster, S., Dandekar, T and Fell, D. (1999) Detection of
elementary flux modes in biochemical networks: a promising tool for
pathway analysis and metabolic engineering, TIBTECH, 17, 53-60.

2. Reder, C. (1988) Metabolic control theory: a structural
approach. J. theor. Biol. 135, 175-201.

3. Schuster, S. and Hilgetag, C. (1994) On elementary flux
modes in biochemical reaction systems at steady state.
J. Biol. Syst. 2, 165-182.

4. Schuster, S., Hilgetag, C., Woods, J. H. and Fell, D. A.
(1996) Elementary modes of functioning in biochemical networks.
In: Computation in Cellular and Molecular Biological Systems
(Cuthbertson, R., Holcombe, M. and Paton, R., eds), pp. 151-165, World
Scientific, Singapore.

5. T. Pfeiffer, I. Sanchez-Valdenebro, J. C. Nuno, F. Montero and S.
Schuster: METATOOL: For Studying Metabolic Networks, Bioinformatics 15
(1999) 251-257.

6. R. Heinrich, H.-G. Holzhuetter and S. Schuster (1987) A theoretical
approach to the evolution and structural design of enzymatic networks;
linear enzymatic chains, branched pathways and glycolysis of
erythrocytes, Bull. Math. Biol. 49, 539-595.

7. R. Heinrich and S. Schuster (1996) The Regulation of Cellular
Systems, Chapman & Hall, New York.

See also
http://www.biologie.hu-berlin.de/biophysics/Theory/tpfeiffer/metatool.html.
