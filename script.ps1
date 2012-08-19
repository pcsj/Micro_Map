$INITIAL_STEP=8
$NUMBER_OF_STEPS_IN_SCAN=100
$NUMBER_OF_STEP_PER_SIM=100
$INCREASE_STEP=0.10
#######
$i=${INITIAL_STEP}
$a=1
if (!(Test-Path "Funzioni_Ottiche")){
new-Item -type directory ("Funzioni_Ottiche")}

if (!(Test-Path "Posizione_Particelle")){
new-Item -type directory ("Posizione_Particelle")}

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
.\FODO_sing_part_1.5.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM > check.txt

if (Test-Path "Funzioni_Ottiche\graph_Funzioni_Ottiche_${i}.png")
{
Remove-Item "Funzioni_Ottiche\graph_Funzioni_Ottiche_${i}.png"
}

if (Test-Path "Posizione_Particelle\graph_Posizione_Particelle_${i}.png")
{
Remove-Item "Posizione_Particelle\graph_Posizione_Particelle_${i}.png"
}

if (Test-Path "graph_Funzioni_Ottiche.png")
{
mv graph_Funzioni_Ottiche.png  Funzioni_Ottiche\graph_Funzioni_Ottiche_${i}.png}

mv graph_Posizione_Particelle.png  Posizione_Particelle\graph_Posizione_Particelle_${i}.png

$i=${i}+$INCREASE_STEP
$a=$a+1
}