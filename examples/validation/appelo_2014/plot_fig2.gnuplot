set terminal pngcairo size 1100,800 enhanced font "sans,11"
set output "/tmp/appelo_fig2.png"
set datafile separator "\t"

set title "Appelo et al. (2014), Figure 2 SUPCRT comparison"
set xlabel "Temperature [C]"
set ylabel "Intrinsic volume [cm^3/mol]"
set xrange [0:200]
set yrange [-42:30]
set grid
set key outside right center

species = "Na+ K+ Mg+2 Ca+2 Cl- HCO3- SO4-2"
current = "/tmp/appelo_fig2.tsv"
reference = "examples/validation/appelo_2014/fig2_supcrt92_reference.tsv"

plot for [i=1:words(species)] current \
       using (strcol(1) eq word(species, i) ? $2 : 1/0):4 \
       with lines linewidth 2 linecolor i title word(species, i), \
     for [i=1:words(species)] reference \
       using (strcol(1) eq word(species, i) ? $2 : 1/0):3 \
       with points pointtype 7 pointsize 0.7 linecolor i notitle
