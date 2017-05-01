#!/bin/bash

totalNumberOfExecutions=1
#next variable can be used to continue previous series of simulations
#all new output files will have indexes beginning with this: _0001.out, _0002.out etc.
initialIndexForOutputs=1

#making links to folders
archived_data=./input_data/archived_data
code_folder=./code
input_data=./input_data
output_data=./output_data
temporary_files=./temporary_files
color_tables=$input_data/color_tables
cooling_tables=$input_data/cooling_tables
observational_data=$input_data/observational_data
parameter_grids=$input_data/parameter_grids
random_seeds=$input_data/random_seeds
DA_colors=$color_tables/DA/model_by_Renedo
DB_colors=$color_tables/DB/model_by_Bergeron
ONe_colors=$color_tables/ONe/model_D
DA_tracks=$cooling_tables/DA/model_A
DB_tracks=$cooling_tables/DB/model_B
ONe_tracks=$cooling_tables/ONe/model_C
observational_LF=$observational_data/luminosity_functions

# 'rm' $output_data/chisquare_test.out
# touch $output_data/chisquare_test.out
# 'rm' $output_data/chi_sq_test_maxregion.out
# touch $output_data/chi_sq_test_maxregion.out

currentNumberOfExecution=1
currentIndexForOutput=$initialIndexForOutputs

while test $currentNumberOfExecution -le $totalNumberOfExecutions
do
   #TODO: remove the files + temp folder in the end of all execs
   # 'rm' $temporary_files/seeds_line.in
   # 'rm' $temporary_files/grid_set_line.in
   # tail -$currentNumberOfExecution $random_seeds/seeds.dat  > $temporary_files/seeds_line.in
   # tail -$currentNumberOfExecution $parameter_grids/grid_set.dat > $temporary_files/grid_set_line.in

   echo "Current execution: " $currentNumberOfExecution "from total of" $totalNumberOfExecutions
   echo "Simulation Nº" $currentIndexForOutput

   gfortran main.f -o main.e

   
   #Parameters of the simulations 
   #TODO figure out what to do about $currentIndexForOutput

   #Various input parameters
   # ln -s $input_data/grid_set_line.in fort.10

   #TRACKS WD CO DA,  Z=0.001
   # ln -s $DA_tracks/Z0001/wd0505_z0001.trk fort.11
   # ln -s $DA_tracks/Z0001/wd0553_z0001.trk fort.12
   # ln -s $DA_tracks/Z0001/wd0593_z0001.trk fort.13
   # ln -s $DA_tracks/Z0001/wd0627_z0001.trk fort.14
   # ln -s $DA_tracks/Z0001/wd0660_z0001.trk fort.15
   # ln -s $DA_tracks/Z0001/wd0692_z0001.trk fort.16
   # ln -s $DA_tracks/Z0001/wd0863_z0001.trk fort.17
   # #TRACKS WD CO DA,  Z=0.01
   # ln -s $DA_tracks/Z001/wd0524_z001.trk fort.21
   # ln -s $DA_tracks/Z001/wd0570_z001.trk fort.22
   # ln -s $DA_tracks/Z001/wd0593_z001.trk fort.23
   # ln -s $DA_tracks/Z001/wd0609_z001.trk fort.24
   # ln -s $DA_tracks/Z001/wd0632_z001.trk fort.25
   # ln -s $DA_tracks/Z001/wd0659_z001.trk fort.26
   # ln -s $DA_tracks/Z001/wd0705_z001.trk fort.27
   # ln -s $DA_tracks/Z001/wd0767_z001.trk fort.28
   # ln -s $DA_tracks/Z001/wd0837_z001.trk fort.29
   # ln -s $DA_tracks/Z001/wd0877_z001.trk fort.30
   # #TRACKS WD CO DA,  Z=0.0382
   # ln -s $DA_tracks/Z003/0524_003_sflhdiff.trk fort.31
   # ln -s $DA_tracks/Z003/0570_003_sflhdiff.trk fort.32
   # ln -s $DA_tracks/Z003/0593_003_sflhdiff.trk fort.33
   # ln -s $DA_tracks/Z003/0610_003_sflhdiff.trk fort.34
   # ln -s $DA_tracks/Z003/0632_003_sflhdiff.trk fort.35
   # ln -s $DA_tracks/Z003/0659_003_sflhdiff.trk fort.36
   # ln -s $DA_tracks/Z003/0705_003_sflhdiff.trk fort.37
   # ln -s $DA_tracks/Z003/1000_003_sflhdiff.trk fort.38
   # #TRACKS WD CO DA,  Z=0.06
   # ln -s $DA_tracks/Z006/0524_006_sflhdiff.trk fort.41
   # ln -s $DA_tracks/Z006/0570_006_sflhdiff.trk fort.42
   # ln -s $DA_tracks/Z006/0593_006_sflhdiff.trk fort.43
   # ln -s $DA_tracks/Z006/0610_006_sflhdiff.trk fort.44
   # ln -s $DA_tracks/Z006/0632_006_sflhdiff.trk fort.45
   # ln -s $DA_tracks/Z006/0659_006_sflhdiff.trk fort.46
   # ln -s $DA_tracks/Z006/0705_006_sflhdiff.trk fort.47
   # ln -s $DA_tracks/Z006/1000_006_sflhdiff.trk fort.48

   # #TRACKS WD CO DB,  Z=0.001
   # ln -s $DB_tracks/Z0001/05047_db_Z=0.001.trk fort.91
   # ln -s $DB_tracks/Z0001/05527_db_Z=0.001.trk fort.92
   # ln -s $DB_tracks/Z0001/059328_db_Z=0.001.trk fort.93
   # ln -s $DB_tracks/Z0001/062738_db_Z=0.001.trk fort.94
   # ln -s $DB_tracks/Z0001/06602_db_Z=0.001.trk fort.95
   # ln -s $DB_tracks/Z0001/069289_db_Z=0.001.trk fort.96
   # ln -s $DB_tracks/Z0001/08637_db_Z=0.001.trk fort.97
   # #TRACKS WD CO DB,  Z=0.01
   # ln -s $DB_tracks/Z001/db_cool_0514.seq fort.101
   # ln -s $DB_tracks/Z001/db_cool_0530.seq fort.102
   # ln -s $DB_tracks/Z001/db_cool_0542.seq fort.103
   # ln -s $DB_tracks/Z001/db_cool_0565.seq fort.104
   # ln -s $DB_tracks/Z001/db_cool_0584.seq fort.105
   # ln -s $DB_tracks/Z001/db_cool_0609.seq fort.106
   # ln -s $DB_tracks/Z001/db_cool_0664.seq fort.107
   # ln -s $DB_tracks/Z001/db_cool_0741.seq fort.108
   # ln -s $DB_tracks/Z001/db_cool_0869.seq fort.109
   # #TRACKS WD CO DB,  Z=0.06
   # ln -s $DB_tracks/Z006/0524db_006_sflhdiff.trk fort.111
   # ln -s $DB_tracks/Z006/0570db_006_sflhdiff.trk fort.112
   # ln -s $DB_tracks/Z006/0593db_006_sflhdiff.trk fort.113
   # ln -s $DB_tracks/Z006/061db_006_sflhdiff.trk fort.114
   # ln -s $DB_tracks/Z006/0632db_006_sflhdiff.trk fort.115
   # ln -s $DB_tracks/Z006/0659db_006_sflhdiff.trk fort.116
   # ln -s $DB_tracks/Z006/070db_006_sflhdiff.trk fort.117
   # ln -s $DB_tracks/Z006/076db_006_sflhdiff.trk fort.118
   # ln -s $DB_tracks/Z006/087db_006_sflhdiff.trk fort.119

   # #TRACKS AND COLORS WD ONe
   # ln -s $ONe_colors/color_106.out fort.121
   # ln -s $ONe_colors/color_110.out fort.122
   # ln -s $ONe_colors/color_116.out fort.123
   # ln -s $ONe_colors/color_120.out fort.124
   # ln -s $ONe_colors/color_124.out fort.125
   # ln -s $ONe_colors/color_128.out fort.126
   # ln -s $ONe_tracks/t106_he.trk fort.127
   # ln -s $ONe_tracks/t110_he.trk fort.128
   # ln -s $ONe_tracks/t116_he.trk fort.129
   # ln -s $ONe_tracks/t120_he.trk fort.130
   # ln -s $ONe_tracks/t128_he.trk fort.131

   # #COLORS WD DB (Bergeron)
   # ln -s $DB_colors/Table_Mass_0.5_DB_sort   fort.132
   # ln -s $DB_colors/Table_Mass_0.6_DB_sort   fort.133
   # ln -s $DB_colors/Table_Mass_0.7_DB_sort   fort.134
   # ln -s $DB_colors/Table_Mass_0.8_DB_sort   fort.135
   # ln -s $DB_colors/Table_Mass_0.9_DB_sort   fort.136
   # ln -s $DB_colors/Table_Mass_1.0_DB_sort   fort.137
   # ln -s $DB_colors/Table_Mass_1.2_DB_sort   fort.138

   # #COLORS WD DA (Renedo)
   # ln -s $DA_colors/cox_0524.dat   fort.61
   # ln -s $DA_colors/cox_0570.dat   fort.62
   # ln -s $DA_colors/cox_0593.dat   fort.63
   # ln -s $DA_colors/cox_0609.dat   fort.64
   # ln -s $DA_colors/cox_0632.dat   fort.65
   # ln -s $DA_colors/cox_0659.dat   fort.66
   # ln -s $DA_colors/cox_0705.dat   fort.67
   # ln -s $DA_colors/cox_0767.dat   fort.68
   # ln -s $DA_colors/cox_0837.dat   fort.69
   # ln -s $DA_colors/cox_0877.dat   fort.70

   # #Observational Luminosity Function
   # ln -s $observational_LF/Althaus_40pc.out fort.71
   
   # #Random seeds
   # # ln -s $input_data/seeds_line.in fort.72
   # ln -s $random_seeds/seeds.dat fort.72



   #Outputs
   #TODO:check the name because it doesn't look appropriate
   # ln -s $output_data/tburstdata.out  fort.81
   # ln -s $output_data/extrapolcolores.out fort.1007
   # ln -s $output_data/average_velocities.out fort.815
   # ln -s $output_data/fl_prueba.out fort.155
   # ln -s $output_data/boot_rowell_thin_$currentNumberOfExecution.out fort.156
   # ln -s $output_data/parametersMC.out  fort.157
   # #Velocities
   # ln -s $output_data/fl_velocidades.out  fort.158
   # ln -s $output_data/fl_velocidades_u.out  fort.400
   # ln -s $output_data/fl_velocidades_v.out  fort.401
   # ln -s $output_data/fl_velocidades_w.out  fort.402
   # ln -s $output_data/flmass$currentNumberOfExecution.out fort.161
   # ln -s $output_data/mass_dis.out fort.162
   # ln -s $output_data/mimf.out fort.667
   # ln -s $output_data/coolda_test.out fort.700

   ln -s ./output.dat fort.420

   # db - DB fraction;
   # g - galactic disc age
   # mf - parameter of initial mass funct/ lifetime mass ratio
   # ifr - parameter of initial to final mass relation
   # bt - timeof burst in star generation process
   # mr - mass reduction factor
   # km - kinematic model 1 or 2
   time ./main.e -db 0.2 -g 9.2 -mf -2.35 -ifr 1 -bt 0.6 -mr 3.7 -km 1
   # time ./main.e 0.901
   # rm fort.*

   currentNumberOfExecution=`expr $currentNumberOfExecution + 1`
   currentIndexForOutput=`expr $currentIndexForOutput + 1`
done

# rm ./main.e
