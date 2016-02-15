import math, os, sys

part = sys.argv[1]
pdb_file  = sys.argv[2]


vmd_input = open(part + '.vmd', 'w')
part = open(part, 'r')


vmd_input.write( '''#!/usr/local/bin/vmd
# VMD script written by graph2py by M. Wolter
# VMD version: 1.9.2
set viewplist {}
set fixedlist {}
mol new %s type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
''' % pdb_file)

res = 1

for line in part:
    frag_num = int(line)
    vmd_input.write( '''mol representation NewRibbons 0.350000 12.000000 5.500000 0
mol color ColorID %d
mol selection {resid %d to %d}
mol material Opaque
mol addrep top
''' % (frag_num, res, res))
    res = res+1

#works only for dp
#for line in part:
#    frag_num = int(line)
#    if frag_num < frag_num_last and frag_num_last >= 0:
#        frag_end = res-1
#        vmd_input.write( '''mol representation NewRibbons 0.350000 12.000000 5.500000 0
#mol color ColorID %d
#mol selection {resid %d to %d}
#mol material Opaque
#mol addrep top
#''' % (colorid, frag_start, frag_end))
#        colorid = colorid + 1
#        frag_start = res
#    res = res+1
#    frag_num_last = frag_num

vmd_input.write( '''unset viewplist
unset fixedlist
proc vmdrestoremycolors {} {
  color scale method RWB
  set colorcmds {
    {color Display {Background} white}
  }
  foreach colcmd $colorcmds {
    set val [catch {eval $colcmd}]
  }
}
vmdrestoremycolors
''')

vmd_input.close
part.close
