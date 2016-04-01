load 3BO5.pdb

# UniProt: SETMR_HUMAN

select preset=(i;57-122)   # UniProt: 73-136
select set=(i;125-249)   # UniProt: 139-263
select postset=(i;269-285)   # UniProt: 283-299
select zinc=(name; ZN)

hide all
show cartoon
show spheres, zinc
color gray, all
color yellow, preset
color orange, set
color red, postset

#background white
set seq_view, 1


