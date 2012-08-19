$INITIAL_STEP=8
$NUMBER_OF_STEPS_IN_SCAN=100
$NUMBER_OF_STEP_PER_SIM=100
$INCREASE_STEP=0.01
#######
$i=${INITIAL_STEP}
$a=1
while ($a -le $NUMBER_OF_STEPS_IN_SCAN){
echo $i
$riga="O 0.0 0.3"
out-file -filepath .\parametri.txt -inputobject $riga -encoding ASCII
$riga="D $i 0.2"
out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
$riga="O 0.0 0.3"
out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
$riga="F 8.0 0.2"
out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
.\Micro_Map.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM > check.txt
mv graph_Funzioni_Ottiche.png graph_Funzioni_Ottiche_${i}.png
mv graph_Posizione_Particelle.png graph_Posizione_Particelle_${i}.png
$i=${i}+$INCREASE_STEP
$a=$a+1
}
