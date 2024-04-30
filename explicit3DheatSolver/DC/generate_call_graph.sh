python gprof2dot.py -help
python gprof2dot.py -e 0.0 prof.outputPar > prog.dotPar
dot -Tpng -o progPar.png prog.dotPar
display progPar.png
