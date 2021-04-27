
# compare MET and MHT together with MET + significance

jes="jecs/jecs_TTbar_PU200.210331.root"
jer="resolutions/TTbar_PU200_210331/resolutions.root"
minBias="data/perfNano_SingleNeutrino_PU200.210331.root"
ttbar="data/perfNano_TTbar_PU200.210331.root"
vbfInv="data/perfNano_VBF_HToInvisible_PU200.210331.root"

#outDir="webpage/test_210417"
outDir="webpage/test_newUnc_210421"
cfg="ch_cfg"
# mht defaults
eta="5.0"
# eta="2.4"
eta="1.5"
pt="30"
pt="10"

sigName=ttbar
#sigName=vbfInv

# for eta in 5 2.4 1.5; do
# for sigName in ttbar vbfInv; do
#for eta in 1.5 2.4 5; do
 
if [[ $sigName == "ttbar" ]]; then
    signal=$ttbar
elif [[ $sigName == "vbfInv" ]]; then
    signal=$vbfInv
fi

# for var in met mht mhtsig metCentral metBarrel mhtsig-mht80 mhtsig-mht60 mhtsig-mht40 mhtsig-mht20 mhtsig-mht0 mht-mhtsig0 mht-mhtsig2.5 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5; do
# for var in met mht mhtsig metCentral metBarrel; do
# # for var in mhtsig-mht80 mhtsig-mht60 mhtsig-mht40 mhtsig-mht20 mhtsig-mht0; do
# for var in mht-mhtsig0 mht-mhtsig2.5 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5; do
# for var in mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do #mhtcorr3
# # for var in mhtcorr-0 mhtcorr-0.5 mhtcorr-1 mhtcorr-1.5 mhtcorr-2 mhtcorr-2.5; do #mhtcorr3
# # for var in mht-mhtsig0 mht-mhtsig2 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5 mhtcorr-0 mhtcorr-0.5 mhtcorr-1 mhtcorr-1.5 mhtcorr-2 mhtcorr-2.5; do #mhtcorr3
# # # # # # var=mht
# for var in mhtsig mhtcorr0 mhtcorr1 mhtcorr2; do #mhtcorr3
# for var in mhtsig mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5 mhtcorr-0 mhtcorr-0.5 mhtcorr-1 mhtcorr-1.5 mhtcorr-2 mhtcorr-2.5; do #mhtcorr3


# for eta in 1.5 2.4 5; do
#var="mhtEvtCorr3"
# # #for var in mhttoy0 mhttoy-1 mhttoy-2 mhttoy1 mhttoy2; do
# # #for var in mht-mhtsig5 mht-mhtsig6 mht-mhtsig7 mht-mhtsig8; do
# # ##for var in mht-mhtsig8; do
# # #for var in mhtcorr0; do #mhtcorr-3
# # for var in mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do #mhtcorr-3
# # # # # # # for pt in 10 15 20 25 30 40 50; do
# # for var in mhttoy-2.5 mhttoy-3; do #mhtcorr-3


for var in mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do
# # #for var in mhtEvtCorr0 mhtEvtCorr0.5 mhtEvtCorr1 mhtEvtCorr1.5 mhtEvtCorr2 mhtEvtCorr2.5 #mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do
thisOutDir=$outDir/$sigName
python scripts/jetHtSuite.py $signal $minBias --jecs $jes --jer $jer \
    -w $cfg $thisOutDir/$var -v $var -e $eta -p $pt

# echo CHTEST writing to $thisOutDir/$var $eta

done
# done
# done

#return 0

rocFList=""
iso10FList=""
iso20FList=""
platroc50FList=""
platroc90FList=""
plateff50FList=""
plateff90FList=""


# for eta in 1.5 2.4 5; do
#     for sigName in ttbar vbfInv; do

if [[ $sigName == "ttbar" ]]; then
    signal=$ttbar
elif [[ $sigName == "vbfInv" ]]; then
    signal=$vbfInv
fi



#for var in met mht metCentral metBarrel; do
#for var in met mht metCentral metBarrel mhtCentral mhtBarrel; do
#for var in mht mhtCentral mhtBarrel; do
# for var in met metCentral metBarrel; do
# for pt in 10 15 20 25 30 40 50; do
    #var=mhtBarrel
    #var=mhtCentral
    # var=mht
# for var in mht-mhtsig0 mht-mhtsig2.5 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5; do
# for var in mht-mhtsig0 mht-mhtsig2 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5 mht-mhtsig6; do
#for var in mhtcorr-0 mhtcorr-0.5 mhtcorr-1 mhtcorr-1.5 mhtcorr-2 mhtcorr-2.5; do #mhtcorr-3
#for var in mht-mhtsig0 mht-mhtsig2 mht-mhtsig3 mht-mhtsig4 mht-mhtsig5; do
#for var in mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do #mhtcorr-3
for var in mhtcorr0 mhtcorr0.5 mhtcorr1 mhtcorr1.5 mhtcorr2 mhtcorr2.5; do
#for var in mhtEvtCorr0 mhtEvtCorr0.5 mhtEvtCorr1 mhtEvtCorr1.5 mhtEvtCorr2 mhtEvtCorr2.5 mhtEvtCorr3; do
#for var in mhtEvtCorr0 mhtEvtCorr0.5 mhtEvtCorr1 mhtEvtCorr1.5 mhtEvtCorr2; do
#for var in mhttoy0 mhttoy-1 mhttoy-2 mhttoy-2.5 mhttoy-3; do
    varL=$var
    # varL=pt$pt
    etaStr="eta"$eta
    ptStr="pt"$pt
    if [[ $var == *"Central"* ]]; then
        etaStr="eta2.4"
    elif [[ $var == *"Barrel"* ]]; then
        etaStr="eta1.5"
    fi
    # if [[ $var == *"mht"* ]]; then
    #     var="mht"
    # fi
    thisOutDir=$outDir/$sigName
    rocFList=$rocFList" ${varL}=${thisOutDir}/${var}/${var}roc-${cfg}_${etaStr}_${ptStr}_gen150.root"
    iso10FList=$iso10FList" ${varL}=${thisOutDir}/${var}/${var}isorate-${cfg}_${etaStr}_${ptStr}_10kHz.root"
    iso20FList=$iso20FList" ${varL}=${thisOutDir}/${var}/${var}isorate-${cfg}_${etaStr}_${ptStr}_20kHz.root"
    platroc50FList=$platroc50FList" ${varL}=${thisOutDir}/${var}/${var}platroc-${cfg}_${etaStr}_${ptStr}_eff0.500.root"
    platroc90FList=$platroc90FList" ${varL}=${thisOutDir}/${var}/${var}platroc-${cfg}_${etaStr}_${ptStr}_eff0.900.root"
    plateff50FList=$plateff50FList" ${varL}=${thisOutDir}/${var}/${var}plateff-${cfg}_${etaStr}_${ptStr}_eff0.50at150.root"
    plateff90FList=$plateff90FList" ${varL}=${thisOutDir}/${var}/${var}plateff-${cfg}_${etaStr}_${ptStr}_eff0.90at150.root"
    lt "${thisOutDir}/${var}/${var}isorate-${cfg}_${etaStr}_${ptStr}_20kHz.root"
done

if [[ $eta == "1.5" ]]; then
    etaLabel="Barrel"
elif [[ $eta == "2.4" ]]; then
    etaLabel="Central"
elif [[ $eta == "5" ]]; then
    etaLabel="Full"
elif [[ $eta == "5.0" ]]; then
    etaLabel="Full"
else
    echo "garbage eta value $eta"
    return 0
fi

suffix=""
if [[ $pt != "30" ]]; then
    suffix="_pt${pt}"
fi

# comboDir=$outDir/combo/${sigName}_metTypeComp/met_and_mht
# comboDir=$outDir/combo/${sigName}_metTypeComp/met_only
#comboDir=$outDir/combo/${sigName}_metTypeComp/mht_only
#comboDir=$outDir/combo/${sigName}_metTypeComp/mhtBarrelPt
#comboDir=$outDir/combo/${sigName}_metTypeComp/mhtCentralPt
#comboDir=$outDir/combo/${sigName}_metTypeComp/mhtFullPt
#comboDir=$outDir/combo/${sigName}_mht-mhtsig
# comboDir=$outDir/combo/${sigName}_mht-mhtsig/mhtFull
# comboDir=$outDir/combo/${sigName}_mht-mhtsig/mhtCentral
# comboDir=$outDir/combo/${sigName}_mht-mhtsig/mhtBarrel
# comboDir=$outDir/combo/${sigName}_mhtcorr/mhtFull
# combodir=$outDir/combo/${sigName}_mhtcorr/mhtCentral
# comboDir=$outDir/combo/${sigName}_mhtcorr/mhtBarrel
# comboDir=$outDir/combo/${sigName}_mhtoy/mhtFull
comboDir=$outDir/combo/${sigName}_mhtcorr${suffix}/mht${etaLabel}
#comboDir=$outDir/combo/${sigName}_mhtEvtCorr/mht$etaLabel
python scripts/stackPlotsFromFiles.py $comboDir/isorate10 L1Puppi_eff $iso10FList --suffix png,pdf
python scripts/stackPlotsFromFiles.py $comboDir/isorate20 L1Puppi_eff $iso20FList --suffix png,pdf
python scripts/stackPlotsFromFiles.py $comboDir/roc Graph --logy --legend TL $rocFList --suffix png,pdf
# python scripts/stackPlotsFromFiles.py $comboDir/platroc50 Graph --logy --legend TR $platroc50FList --suffix png,pdf
# python scripts/stackPlotsFromFiles.py $comboDir/platroc90 Graph --logy --legend TR $platroc90FList --suffix png,pdf
# python scripts/stackPlotsFromFiles.py $comboDir/plateff50 L1Puppi_eff $plateff50FList --suffix png,pdf
# python scripts/stackPlotsFromFiles.py $comboDir/plateff90 L1Puppi_eff $plateff90FList --suffix png,pdf

# done    
# done


    # met=$thisOutDir/met/metroc-ch_cfg_eta5.0_pt30_gen150.root \
    #     mht=$thisOutDir/mht/mhtroc-ch_cfg_eta5.0_pt30_gen150.root \
    #     mhtsig=$thisOutDir/mhtsig/mhtsigroc-ch_cfg_eta5.0_pt30_gen150.root \
    #     mhtsig-mht100=$thisOutDir/mhtsig-mht100/mhtsig-mht100roc-ch_cfg_eta5.0_pt30_gen150.root 



# for val in 30 50 80 100
# do
# echo $val
# python scripts/jetHtSuite.py $ttbar $minBias --jecs $jes --jer $jer \
#     -w ch_cfg webpage/test_210412/mhtsig-mht$val -v mhtsig-mht$val -e 5
# done








# plotting distributions
# python scripts/jetHtSuite.py data/perfNano_TTbar_PU200.210331.root data/perfNano_SingleNeutrino_PU200.210331.root --jecs jecs/jecs_TTbar_PU200.210331.root --jer resolutions/TTbar_PU200_210331/resolutions.root -w ch_cfg webpage/test_210412/dist -v metsig -e 5 -P dist


# for val in 0 0.5 1 2
# do
# echo $val
# python scripts/jetHtSuite.py data/perfNano_TTbar_PU200.210331.root data/perfNano_SingleNeutrino_PU200.210331.root \
#     --jecs jecs/jecs_TTbar_PU200.210331.root --jer resolutions/TTbar_PU200_210331/resolutions.root \
#     -w ch_cfg webpage/test_210412/mhtcorr$val -v mhtcorr$val -e 5
# done


# odir="combo_mhtcorr"
# for sfx in png pdf
# do
# python scripts/stackPlotsFromFiles.py webpage/test_210412/$odir/isorate20.$sfx L1Puppi_eff webpage/test_210412/mhtcorr*/mhtcorr*isorate*_20kHz.root
# #python scripts/stackPlotsFromFiles.py webpage/test_210412/combo1/roc.$sfx Graph webpage/test_210412/mhtcorr*/mhtcorr*roc*.root --logy
# python scripts/stackPlotsFromFiles.py webpage/test_210412/$odir/roc.$sfx Graph --logy --legend TL\
#     met=webpage/test_210412/mhtcorr0/mhtcorr0roc-ch_cfg_eta5.0_pt30_gen150.root \
#     0.5sigma=webpage/test_210412/mhtcorr0.5/mhtcorr0.5roc-ch_cfg_eta5.0_pt30_gen150.root \
#     1sigma=webpage/test_210412/mhtcorr1/mhtcorr1roc-ch_cfg_eta5.0_pt30_gen150.root \
#     2sigma=webpage/test_210412/mhtcorr2/mhtcorr2roc-ch_cfg_eta5.0_pt30_gen150.root 
# done





# ./scripts/prun.sh runPerformanceNTuple.py --110X_v2 TTbar_PU200 TTbar_PU200.210331

# ./scripts/prun.sh runPerformanceNTuple.py --110X_v2 SingleNeutrino_PU200 SingleNeutrino_PU200.210331

# python scripts/makeJecs.py data/perfNano_TTbar_PU200.210331.root -A -o jecs_TTbar_PU200.210331.root

# python scripts/jetHtSuite.py data/perfNano_TTbar_PU200.210331.root data/perfNano_SingleNeutrino_PU200.210331.root \
#     plots_dir -w l1pfpu -v met --jecs jecs/jecs_TTbar_PU200.210331.root

# # ntuple for resolution plots
# ./scripts/prun.sh runRespNTupler.py --110X_v2 TTbar_PU200 TTbar_PU200.110X_v2
# # make resolution plots
# python scripts/respPlots.py respTupleNew_TTbar_PU200.110X_v2.root resolutions/TTbar_PU200_210331 -w ch_cfg -p jet --no-eta --writeResolutions




# python scripts/jetHtSuite.py data/perfNano_TTbar_PU200.210331.root data/perfNano_SingleNeutrino_PU200.210331.root --jecs jecs/jecs_TTbar_PU200.210331.root -w ch_cfg webpage/test_210412/mhtsig -v mhtsig --jer resolutions/TTbar_PU200_210331/resolutions.root

# -v mhtsig-mht100

# #combine plots from multiple files
# python scripts/stackPlotsFromFiles.py webpage/test_210412/combo1/t.pdf L1Puppi_eff webpage/test_210412/mhtcorr*/mhtcorr*isorate-*_20kHz.root
