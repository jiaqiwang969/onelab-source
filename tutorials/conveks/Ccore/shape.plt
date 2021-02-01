#set terminal pdf

set grid
set ytics nomirror
set y2tics
set style data lines

set xlabel "tau"
set y2label "w"
set ylabel "dw/dtau"

set output "aaa.pdf"
plot "aaa" u 2:3 axis x1y2 t"w",\
     "aaa" u 2:4 t"dwdt",\
     "aaa" u 2:5 w linesp t"dwdt FD"

pause -1

set output "optishape1.pdf"
plot "Torque_linTH.txt" u 2:3 axis x1y2 t"w",\
     "Torque_linTH.txt" u 2:4 t"dwdt",\
     "Torque_linTH.txt" u 2:5 w linesp t"dwdt FD"

pause -1

set output "optishape2.pdf"
plot "LieDevSmoother.txt" u 2:3 axis x1y2 t"w",\
     "LieDevSmoother.txt" u 2:4 t"dwdt",\
     "LieDevSmoother.txt" u 2:5 w linesp t"dwdt FD"

pause -1

set output "optishape3.pdf"
plot "Torque_linTH2.txt" u 2:3 axis x1y2 t"w",\
     "Torque_linTH2.txt" u 2:4 t"dwdt",\
     "Torque_linTH2.txt" u 2:5 w linesp t"dwdt FD"

pause -1

set output "optishape4.pdf"
plot "LieDevSmoother2.txt" u 2:3 axis x1y2 t"w",\
     "LieDevSmoother2.txt" u 2:4 t"dwdt",\
     "LieDevSmoother2.txt" u 2:5 w linesp t"dwdt FD"

pause -1

set output "optishape5.pdf"
plot "Torque_TH.txt" u 2:3 axis x1y2 t"w",\
     "Torque_TH.txt" u 2:4 t"dwdt",\
     "Torque_TH.txt" u 2:5 w linesp t"dwdt FD"

pause -1

