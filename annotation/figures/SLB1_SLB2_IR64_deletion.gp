set terminal png tiny size 1400,1400
set output "0_annotation/figures/SLB1_SLB2_IR64_deletion.png"
set size 1,1
set grid
unset key
set border 15
set tics scale 0
set xlabel "Nipponbare#1#chr1:28966552-29064703"
set ylabel "IR64#1#chr1:30213095-30306203"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:98152]
set yrange [1:93109]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "0_annotation/figures/SLB1_SLB2_IR64_deletion.fplot" title "FWD" w lp ls 1, \
 "0_annotation/figures/SLB1_SLB2_IR64_deletion.rplot" title "REV" w lp ls 2
