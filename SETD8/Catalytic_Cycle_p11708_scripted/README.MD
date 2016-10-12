1) run master.sh (bash master.sh in interactive job) --> this will generate solvated structures
2) run equilibration scripts in equil_scripts individually per structure (qsub X.sh to run GPU batch jobs) - this will generate packaged (in OpenMM 7.0) structures
3) run package.sh (bash package.sh in interactive job) --> this will repackage all structures in OpenMM 6.3
4) Ready for FAH!
