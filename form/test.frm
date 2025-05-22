#-
#: IncDir procedures
#: SmallExtension 100M
#: WorkSpace      2G
#: MaxTermSize 4M
Off Statistics;

*#include declarations.h
cf log(symmetric);
v q;
i i1, i2, i3;
i c11, c2, c12;
auto s marker;
cf func;
s sym;
cf summe, summeP;


l test = log(i1, q);
.sort
if(occurs(i2)) mul marker;

print test;
.end

