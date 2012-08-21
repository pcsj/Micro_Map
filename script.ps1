$GRAD_DEF=8.0
$GRAD_FOC=8.0
$LENGTH_TOTAL=1.5
$LENGTH_F=0.1
$NUMBER_OF_STEPS_IN_SCAN=10
$NUMBER_OF_STEPS_DRIFT=10
$NUMBER_OF_STEPS_FOC=10
$INCREASE_STEP_GRAD=2
$INCREASE_STEP_LENG=0.05
$NUMBER_OF_STEP_PER_SIM=100
$ENERGY=30
$X_MAX=1.5
$Y_MAX_POS=0.1
$Y_MAX_OPT=5
$vero=$FALSE
$P_RIF=1
$X_RIF=2.5
$PERC=0.2

#######
$grad_d=${GRAD_DEF}
$grad_f=${GRAD_FOC}
$lung_elem=${LENGTH_F}
$lung_drift_i=0.10
$energy=${ENERGY}
$a=1
$i=1
$j=1
$k=1
if (!(Test-Path "Funzioni_Ottiche")){
new-Item -type directory ("Funzioni_Ottiche")}

if (!(Test-Path "Posizione_Particelle")){
new-Item -type directory ("Posizione_Particelle")}


while ($a -le $NUMBER_OF_STEPS_FOC)
{
$grad_d=${GRAD_DEF}
while ($i -le $NUMBER_OF_STEPS_FOC)
{
$lung_elem=${LENGTH_F}
while ($j -le $NUMBER_OF_STEPS_IN_SCAN)
{
$lung_drift_m=0.5
while ($k -le $NUMBER_OF_STEPS_DRIFT)
{
	echo $grad_d
	echo $grad_f
	echo $lung_elem
	echo $lung_drift_m

if (!(Test-Path ".\Posizione_Particelle\FOC_${grad_f}_e_DEFOC_${grad_d}")){
new-Item -type directory (".\Posizione_Particelle\FOC_${grad_f}_e_DEFOC_${grad_d}")}

if (!(Test-Path "Funzioni_Ottiche\FOC_${grad_f}_e_DEFOC_${grad_d}")){
new-Item -type directory (".\Funzioni_Ottiche\FOC_${grad_f}_e_DEFOC_${grad_d}")}

if (!(Test-Path ".\Posizione_Particelle\FOC_${grad_f}_e_DEFOC_${grad_d}\LUNG_CELL_${lung_elem}")){
new-Item -type directory (".\Posizione_Particelle\FOC_${grad_f}_e_DEFOC_${grad_d}\LUNG_CELL_${lung_elem}")}

if (($lung_drift_m+$lung_drift_i) -gt (${LENGTH_TOTAL}+$lung_elem+$lung_elem))
{
$lung_drift_m=0.5
$vero=$TRUE
} 

	$riga="O 0.0 0.1"
	out-file -filepath .\parametri.txt -inputobject $riga -encoding ASCII
    
    $riga="F $grad_f $lung_elem"
	out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	$riga="O 0.0 $lung_drift_m"
	out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	$riga="D $grad_d $lung_elem"
	out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
    
    $lung_drift_i=${LENGTH_TOTAL}-$lung_elem-$lung_elem-$lung_drift_m
    
	$riga="O 0.0 $lung_drift_i"
	out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
	.\FODO_sing_part_1.7.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM -xmax_opt $X_MAX -xmax_pos $X_MAX -compare_X $X_RIF -compare_P $P_RIF -% $PERC
    
	if (Test-Path ".\graph_Funzioni_Ottiche.png")
	{
	Move-Item -Force graph_Funzioni_Ottiche.png  .\Funzioni_Ottiche\FOC_${grad_f}_e_DEFOC_${grad_d}\lung_cel_${lung_elem}__drift_${lung_drift_m}.png
	}

    Move-Item -Force graph_Posizione_Particelle.png  .\Posizione_Particelle\FOC_${grad_f}_e_DEFOC_${grad_d}\LUNG_CELL_${lung_elem}\drift_${lung_drift_m}.png


echo "Sono l'iteratore del primo ciclo "$k
if ($vero)
{$lung_drift_m=$lung_drift_m-$INCREASE_STEP_LENG}
$lung_drift_m=$lung_drift_m+$INCREASE_STEP_LENG
$k=$k+1
}
$vero=$FALSE
$k=1    
echo "Sono l'iteratore del secondo ciclo "$j
$lung_elem=$lung_elem+$INCREASE_STEP_LENG
$j=$j+1
}
$j=1
echo "Sono l'iteratore del terzo ciclo "$i
$grad_d=$grad_d+$INCREASE_STEP_GRAD
$i=$i+1
}
$i=1
echo "Sono l'iteratore del quarto ciclo "$a
$grad_f=$grad_f+$INCREASE_STEP_GRAD
$a=$a+1
}
