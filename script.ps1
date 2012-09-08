$GRAD_DEF=8.0
$GRAD_FOC=8.0
$LENGTH_TOTAL=1.5
$LENGTH_F=0.1
$NUMBER_OF_STEPS_IN_SCAN=10
$NUMBER_OF_STEPS_DRIFT=10
$NUMBER_OF_STEPS_FOC=10
$INCREASE_STEP_GRAD=2
$INCREASE_STEP_LENG=0.05
$NUMBER_OF_STEP_PER_SIM=500
$ENERGY=30
$X_MAX=1.5
$Y_MAX_POS=0.1
$Y_MAX_OPT=5
$X_P_RIF=1
$X_X_RIF=2.5
$Y_P_RIF=1
$Y_X_RIF=2.5
$CONFRONTO_ASSI_X=0.5
$CONFRONTO_ASSI_P=0.1
$PERC=0.03
$NUMBER_OF_STEPS_ENERGY=10;
#######
$grad_d=${GRAD_DEF}
$grad_f=${GRAD_FOC}
$lung_elem=${LENGTH_F}
$energy=${ENERGY}
$a=0
$i=0
$j=0
$g=0
$k=0
$vero=$False
$prova=$false
$conto=0

Invoke-BatchFile 'C:\Program Files\Microsoft Visual Studio 11.0\VC\vcvarsall.bat'
cd ..\FODO_sing_part_1.8\
cl /EHsc .\FODO_sing_part.cpp
mv -Force FODO_sing_part.exe ..\Test\
cd ..\Test


if (!(Test-Path "Funzioni_Ottiche")){
new-Item -type directory ("Funzioni_Ottiche")}

if (!(Test-Path "Posizione_Particelle")){
new-Item -type directory ("Posizione_Particelle")}

if (!(Test-Path "Math_rilevati")){
new-Item -type directory ("Math_rilevati")}

if (!(Test-Path "Funzioni_Ottiche_T")){
new-Item -type directory ("Funzioni_Ottiche_T")}

if (!(Test-Path "Ellissi")){
new-Item -type directory ("Ellissi")}

if (!(Test-Path "Ellissi_T")){
new-Item -type directory ("Ellissi_T")}


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

            if ((($conto) -eq (1))-and($vero))
            {
                $lung_drift_m=0.5
                echo "io sono dentro l'if che mi fissa di nuovo la lunghezza del drift"
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
            .\FODO_sing_part.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM -compare_X $X_X_RIF $X_P_RIF -compare_Y $Y_X_RIF $Y_P_RIF -perc $PERC
    	    #-xmax_opt $X_MAX -xmax_pos $X_MAX

            if (Test-Path ".\graph_Funzioni_Ottiche.png")
            {
	           New-Item -Force -type directory ".\Funzioni_Ottiche\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}"              
	           Move-Item -Force graph_Funzioni_Ottiche.png  .\Funzioni_Ottiche\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	        }
            if (Test-Path ".\graph_Funzioni_Ottiche_T.png")
            {
	           New-Item -Force -type directory ".\Funzioni_Ottiche_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}"
	           Move-Item -Force graph_Funzioni_Ottiche_T.png  .\Funzioni_Ottiche_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	        }
            if (Test-Path ".\graph_Posizione_Particelle.png")
	        {
	           New-Item -Force -type directory ".\Posizione_Particelle\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}"
               Move-Item -Force graph_Posizione_Particelle.png  .\Posizione_Particelle\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
               $prova=$true            
            }
            if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche.png")
	        {
	           New-Item -Force -type directory ".\Ellissi\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}"
               Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche.png  .\Ellissi\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
            }
            if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche_T.png")
	        {
	           New-Item -Force -type directory ".\Ellissi_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}"
               Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche_T.png  .\Ellissi_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
            }
            if (Test-Path ".\Math_rilevati.txt")
	        {            
            Move-Item -Force Math_rilevati.txt  .\Math_rilevati\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}.txt           
            }
   
            if ($prova)
            {
                while ($g -le $NUMBER_OF_STEPS_ENERGY)
                {
                    $energy=$energy+2;
                    echo $energy
	                $riga="1"
	                out-file -filepath .\inputdata.txt -inputobject $riga -encoding ASCII
	                $riga="$energy"
	                out-file -filepath .\inputdata.txt -inputobject $riga -append -encoding ASCII
	                $riga="0.0 0.0 0.05 0.05"
	                out-file -filepath .\inputdata.txt -inputobject $riga -append -encoding ASCII

                    .\FODO_sing_part.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM -compare_X $X_X_RIF $X_P_RIF -compare_Y $Y_X_RIF $Y_P_RIF -perc 0.1
                    
                    if (Test-Path ".\graph_Posizione_Particelle.png")
	                {
                       Move-Item -Force graph_Posizione_Particelle.png  .\Posizione_Particelle\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
                    }
                    if (Test-Path ".\graph_Funzioni_Ottiche.png")
                    {            
	                   Move-Item -Force graph_Funzioni_Ottiche.png  .\Funzioni_Ottiche\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	                }
                    if (Test-Path ".\graph_Funzioni_Ottiche_T.png")
                    {
	                   Move-Item -Force graph_Funzioni_Ottiche_T.png  .\Funzioni_Ottiche_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	                }
                    if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche.png")
	                {
                       Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche.png  .\Ellissi\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
                    }
                    if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche_T.png")
	                {
                       Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche_T.png  .\Ellissi_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
                    }                    
                    
                    $g=$g+1
                }
                $energy=${ENERGY}
            }

            if ((($lung_drift_i) -lt (0))-or($vero))
            {          
                $lung_drift_m=$lung_drift_m-$INCREASE_STEP_LENG
                $vero=$True
                $conto=$conto+1
                echo "io sono dentro l'if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            }
            else
            {
                $lung_drift_m=$lung_drift_m+$INCREASE_STEP_LENG
            }
            $k=$k+1
            
            }
            $lung_drift_m=0.5
            $k=1    
            echo "Sono l'iteratore del secondo ciclo "$j
            $lung_elem=$lung_elem+$INCREASE_STEP_LENG
            $vero=$False
            $conto=0
            $prova=$false
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