# This script was generated by OpenMMDL-Setup on 2023-11-12.


        #       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
        #     .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      
        #    / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      
        #   ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    
        #   |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
        #   : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
        #    \ `"/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  
        #     '. \_/``".'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ 
        #       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` 
                                                                                                      
                                                                                                      
    
#!/bin/bash

################################## Receptor ######################################
rcp_nm=8EFO_protein # the file name of ligand without suffix `pdb`
rcp_ff=leaprc.protein.ff19SB


## Clean the PDB file by pdb4amber
pdb4amber -i ${rcp_nm}.pdb -o ${rcp_nm}_amber.pdb

## `tleap` requires that all residues and atoms have appropriate types to ensure compatibility with the specified force field.
## To avoid `tleap` failing, we delete non-essential atoms, such as hydrogens, but preserve important atoms like carbon and nitrogen within the caps residues.
## Don' worry about the missing atoms as tleap has the capability to reconstruct them automatically. 
awk '! ($2 ~ "(CH3|HH31|HH32|HH33)" || $3 ~ "(CH3|HH31|HH32|HH33)" )' ${rcp_nm}_amber.pdb > ${rcp_nm}_amber_f.pdb 
grep -v '^CONECT' ${rcp_nm}_amber_f.pdb > ${rcp_nm}_cnt_rmv.pdb

################################## Ligand ######################################
# Normal Ligand
nmLigFile=8QY # the file name of ligand without suffix `pdb`
charge_method=bcc # refers to the charge method that antechamber will adopt
charge_value=1 # Enter the net molecular charge of the ligand as integer (e.g. 1 or -2)
lig_ff=gaff2 # Ligand force field 

## Clean the PDB file by pdb4amber
pdb4amber -i ${nmLigFile}.pdb -o ${nmLigFile}_amber.pdb

## Generate a prepc file and an additional frcmod file by `antechamber`
antechamber -fi pdb -fo prepc -i ${nmLigFile}_amber.pdb -o ${nmLigFile}.prepc -c ${charge_method} -at ${lig_ff} -nc ${charge_value} -pf y
parmchk2 -f prepc -i ${nmLigFile}.prepc -o ${nmLigFile}.frcmod

## Rename ligand pdb
antechamber -i ${nmLigFile}.prepc -fi prepc -o rename_${nmLigFile}.pdb -fo pdb

######################  Combine All Components to Be Modelled ####################
cat > tleap.combine.in <<EOF

source ${rcp_ff}
source leaprc.${lig_ff}

loadamberprep ${nmLigFile}.prepc
loadamberparams ${nmLigFile}.frcmod

rcp = loadpdb ${rcp_nm}_cnt_rmv.pdb
nmLig = loadpdb rename_${nmLigFile}.pdb 
comp = combine{rcp nmLig}
savepdb comp comp.pdb

quit

EOF

tleap -s -f tleap.combine.in > tleap.combine.out
grep -v '^CONECT' comp.pdb > comp_cnt_rmv.pdb

################################ Add Water/Membrane ##############################
lipid_tp=POPC
lipid_ratio=1
lipid_ff=lipid21
dist2Border=15  # The minimum distance between the maxmin values for x y and z to the box boundaries. Flag --dist
padDist=17  # The width of the water layer over the membrane or protein in the z axis. Flag --dist_wat
water_ff=opc
pos_ion=Na+
neg_ion=Cl-
ionConc=0.15


## Build the membrane
packmol-memgen --pdb comp.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate

## Clean the complex pdb by `pdb4amber` for further `tleap` process
pdb4amber -i bilayer_comp.pdb -o clean_bilayer_comp.pdb
grep -v '^CONECT' clean_bilayer_comp.pdb > clean_bilayer_comp_cnt_rmv.pdb


##################### Generate Prmtop and Frcmod File for the Complex ###################### 
cat > tleap.in <<EOF

source ${rcp_ff}
source leaprc.water.${water_ff}
source leaprc.${lig_ff}
source leaprc.${lipid_ff}

loadamberprep ${nmLigFile}.prepc
loadamberparams ${nmLigFile}.frcmod

system = loadpdb clean_bilayer_comp_cnt_rmv.pdb

setBox system "vdw"
check system 
charge system

savepdb system system.${water_ff}.pdb
saveamberparm system system.${water_ff}.prmtop system.${water_ff}.inpcrd

quit

EOF

tleap -s -f tleap.in > tleap.out