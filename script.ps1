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
$vero=0
$P_RIF=1
$X_RIF=2.5
$PERC=0.1
#######
$grad_d=${GRAD_DEF}
$grad_f=${GRAD_FOC}
$lung_elem=${LENGTH_F}
$energy=${ENERGY}
$a=0
$i=0
$j=0
$k=0

Invoke-BatchFile 'C:\Program Files\Microsoft Visual Studio 11.0\VC\vcvarsall.bat'
cd ..\FODO_sing_part_1.8\
cl /EHsc .\FODO_sing_part.cpp
mv -Force FODO_sing_part.exe ..\Release\
cd ..\Release


if (!(Test-Path "Funzioni_Ottiche")){
new-Item -type directory ("Funzioni_Ottiche")}

if (!(Test-Path "Posizione_Particelle")){
new-Item -type directory ("Posizione_Particelle")}

if (!(Test-Path "Math_rilevati")){
new-Item -type directory ("Math_rilevati")}


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

	        echo $k            
	        echo $j
	        echo $i
            echo $a
            

            if (($lung_drift_m+$lung_drift_i) -gt (${LENGTH_TOTAL}+$lung_elem+$lung_elem))
            {
                $lung_drift_m=0.5
                $vero=1
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
            .\FODO_sing_part.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM -compare_X $X_RIF -compare_P $P_RIF -perc $PERC
    	    #-xmax_opt $X_MAX -xmax_pos $X_MAX
        
            if (Test-Path ".\graph_Funzioni_Ottiche.png")
            {
	           Move-Item -Force graph_Funzioni_Ottiche.png  .\graph_Funz_${a}_${i}_${j}_${k}.png
	        }
            if (Test-Path ".\graph_Posizione_Particelle.png")
	        {
               Move-Item -Force graph_Posizione_Particelle.png  .\graph_Pos_${a}_${i}_${j}_${k}.png
            }
            if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche.png")
	        {
               Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche.png  .\graph_Ell_${a}_${i}_${j}_${k}.png
            }
            if (Test-Path ".\Math_rilevati.txt")
	        {            
            Move-Item -Force Math_rilevati.txt  .\Math_rilevati\Math_${a}_${i}_${j}_${k}.txt            
            }
            echo "Sono l'iteratore del primo ciclo "$k
            if ($vero)
            {
                $lung_drift_m=$lung_drift_m-$INCREASE_STEP_LENG
            }
            $lung_drift_m=$lung_drift_m+$INCREASE_STEP_LENG
            $k=$k+1
            
            }
            $vero=0
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

$a=0
$i=0
$j=0
$k=0
$grad_d=${GRAD_DEF}
$grad_f=${GRAD_FOC}
$lung_elem=${LENGTH_F}
$lung_drift_m=0.5

     
while ($a -le $NUMBER_OF_STEPS_FOC)
{
    while ($i -le $NUMBER_OF_STEPS_FOC)
    {
        while ($j -le $NUMBER_OF_STEPS_FOC)
        {
            while ($k -le $NUMBER_OF_STEPS_FOC)
            {
            	if (Test-Path ".\graph_Funz_${a}_${i}_${j}_${k}.png")
	            {
	                Move-Item -Force graph_Funz_${a}_${i}_${j}_${k}.png  .\Funzioni_Ottiche\FOC_${grad_f+a*$INCREASE_STEP_GRAD}__DEFOC_${grad_d+i*$INCREASE_STEP_GRAD}__L_cel_${lung_elem+j*$INCREASE_STEP_LENG}__drift_${lung_drift_m+k*$INCREASE_STEP_LENG}.png
                    #bisogna iterare anche i parametri grad_f ecc, basta moltiplicare il dato iniziale che deve essere riportato a quello iniziale per iteratore moltiplicato per incremento 
	            }
	            if (Test-Path ".\graph_Pos_${a}_${i}_${j}_${k}.png")
	            {
                    Move-Item -Force graph_Pos_${a}_${i}_${j}_${k}.png  .\Posizione_Particelle\FOC_${grad_f+a*$INCREASE_STEP_GRAD}__DEFOC_${grad_d+i*$INCREASE_STEP_GRAD}__L_cel_${lung_elem+j*$INCREASE_STEP_LENG}__drift_${lung_drift_m+k*$INCREASE_STEP_LENG}.png
                }
	            if (Test-Path ".\graph_Ell_${a}_${i}_${j}_${k}.png")
	            {
                    Move-Item -Force graph_Ell_${a}_${i}_${j}_${k}.png  .\Ellissi\FOC_${grad_f+a*$INCREASE_STEP_GRAD}__DEFOC_${grad_d+i*$INCREASE_STEP_GRAD}__L_cel_${lung_elem+j*$INCREASE_STEP_LENG}__drift_${lung_drift_m+k*$INCREASE_STEP_LENG}.png
                }
                $k=$k+1
            }
            $k=0
            $j=$j+1
        }
        $j=0
        $i=$i+1
    }
    $i=0
    $a=$a+1
}