#!/usr/bin/env python
# -*- coding: utf-8 -*-
###### Modules #################################################################
from unittest import TestProgram
from UDFManager import *

import argparse
import codecs
import os
import platform
import shutil
import sys

import dynamic_rheo_setup.values as val
################################################################################
## MAIN
################################################################################
def setup():
	# 各種条件を読み取り
	read_all()
	# セットアップ
	setup_all()

	# make_udfs()
	# summary_fname, job_dir = Setup(target_udf, base_udf, Setting)
	# #
	# Make_Series_Calc(Setting, Parameters, summary_fname, base_udf, job_dir)


################################################################################

################################################################################
# Functions
################################################################################
###################################
# 各種条件を読み取り
def read_all():
	read_arg()
	read_sim_cond()
	return


def read_arg():
	parser = argparse.ArgumentParser(description='Select udf file to read !')
	parser.add_argument('udf', help="udf file name to read previous simulation")
	args = parser.parse_args()
	if args.udf:
		if len(args.udf.split('.')) != 2 or args.udf.split('.')[1] != 'udf':
			print('\nthe file name you selected is not udf file !')
			sys.exit('select proper udf file to read.')
		elif not os.access(args.udf, os.R_OK):
			sys.exit('\nSelected udf of ', args.udf, ' seems not exist !\nbye now!!')
		else:
			val.read_udf = args.udf
			# print('Selected udf file is ' + val.read_udf)
	else:
		print('no udf file is selected')
		sys.exit('select proper udf file to read.')
	return

# シミュレーション条件を設定する。
def read_sim_cond():
	if not os.path.isfile('../dynamic_rheo.udf'):
		print('\nIn the parent directory, no "dynamic_rheo.udf" is found !')
		print('New one will be generated.')
		print('Please, modify and save it !\n')
		make_newudf()
		input('Press ENTER to continue...')
	else:
		read_and_set()
	return

# make new udf when not found.
def make_newudf():
	contents = '''
	\\begin{def}
		CalcConditions:{
			Cognac_ver:select{"cognac112"} "使用する Cognac のバージョン",
			Cores: int "計算に使用するコア数を指定"
			} "Cognac による計算の条件を設定"
		DeformationMode:{
			Deformation:select{"Stretch", "Shear"} "変形モードを選択"
			}
		Dynamics:{
			SweepMode:select{"none", "StrainSweep", "FrequencySweep"} "Sweep モードを選択",
				StrainSweep:{
					StrainSweepConditions[]:{
						Temperature:float "",
						MinStrain:float "最小ひずみ",
						MaxStrain:float "最大ひずみ",
						Frequency:float "測定周波数",
						Data_per_Digit:int "対数スケールで、一桁あたりの測定数"
						}
					}
				FrequencySweep:{
					FrequencySweepConditions[]:{
						Temperature:float "",
						MinFrequency:float "最小周波数",
						MaxFrequency:float "最大周波数",
						Strain:float "ひずみ",
						Data_per_Digit:int "対数スケールで、一桁あたりの測定数"
						}
					}
			} "計算ターゲットの条件を設定"		
		SimulationParameters:{
			Cycles:int "",
			Output_per_Cycle:int "",
			} "シミュレーション条件を設定"
		AnalysisParameters:{
			SkipCycles:int "",
			} "分析条件を設定"
	\end{def}	

	\\begin{data}
		CalcConditions:{"cognac112",2}
		DeformationMode:{"Shear"}
		Dynamics:{
		"FrequencySweep",
		{
		[{1.0000000,9.9999998e-03,1.0000000,1.0000000,5}]
		}
		{
		[{1.0000000,9.9999998e-03,1.0000000,9.9999998e-03,5}{1.2000000,9.9999998e-03,1.0000000,9.9999998e-03,5}]
		}
		}
		SimulationParameters:{25,20}
		AnalysisParameters:{5}
	\end{data}
	'''
	###
	with codecs.open('../dynamic_rheo.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

# Read udf and setup initial conditions
def read_and_set():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		read_condition()
		# select
		init_calc()
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		return
	else:
		sys.exit("##### \nQuit !!")

# Read condition udf
def read_condition():
	u = UDFManager('../dynamic_rheo.udf')
	u.jump(-1)
	# 使用するCognacのバージョン
	val.ver_Cognac = u.get('CalcConditions.Cognac_ver')
	# 計算に使用するコア数
	val.core = u.get('CalcConditions.Cores')
	# Deformation mode
	val.deform_mode  = u.get('DeformationMode.Deformation')
	# SweepConditions
	val.sweep_mode = u.get('Dynamics.SweepMode')
	if val.sweep_mode == 'StrainSweep':
		val.StrainSweepConditions = u.get('Dynamics.StrainSweep.StrainSweepConditions[]')
	elif val.sweep_mode == 'FrequencySweep':
		val.FrequencySweepConditions = u.get('Dynamics.FrequencySweep.FrequencySweepConditions[]')
	# SimulationParameters
	val.total_cycles = u.get('SymulationParameters.Cycles')
	val.output_cycle = u.get('SymulationParameters.Output_per_Cycle')
	print(val.total_cycles, val.output_cycle)
	# AnalysisParameters
	val.skipcycles = u.get('AnalysisParameters.SkipCycles')
	return
# 
def init_calc():
	text = "################################################\n"
	text += f"Cores used for simulation\t\t{val.core:.2g}\n"
	text += "################################################\n"
	text += f"Deform mode:\t\t\t\t{val.deform_mode:}\n"
	text += "################################################\n"
	if val.sweep_mode == 'StrainSweep':
		text += f"Sweep mode:\t\t\t{val.sweep_mode}\n"
		for i, data in enumerate(val.StrainSweepConditions):
			text += f'# {i}\n'
			text += f"temperature:\t\t\t\t{data[0]:.4g}\n"
			text += f"Minimum Strain:\t\t\t\t{data[1]:.4g}\n"
			text += f"Maximum Strain:\t\t\t\t{data[2]:.4g}\n"
			text += f"Frequency:\t\t\t\t{data[3]:.4g}\n"
			text += f"Data per Digit:\t\t\t\t{data[4]:.4g}\n"
	elif val.sweep_mode == 'FrequencySweep':
		text += "Sweep mode:\t\t\t" + str(val.sweep_mode) + "\n"
		for i, data in enumerate(val.FrequencySweepConditions):
			text += f'# {i}\n'
			text += f'  temperature:\t\t\t\t{data[0]:.2g}\n'
			text += f"  Minimum Frequency:\t\t\t{data[1]:.4g}\n"
			text += f"  Maximum Frequency:\t\t\t{data[2]:.4g}\n"
			text += f"  Strain:\t\t\t\t{data[3]:.4g}\n"
			text += f"  Data per Digit:\t\t\t{data[4]:.4g}\n"
	text += "################################################\n"
	text += 'Simulation Parameters:\n'
	text += f"  Cycles:\t\t\t\t{val.total_cycles:}\n"
	text += f"  Output_per_Cycle:\t\t\t{val.output_cycle:}\n"
	text += 'Analysis Parameters:\n'
	text += f"  Skip Cycles:\t\t\t\t{val.skipcycles:}\n"
	text += "################################################\n"
	print(text)
	return






def setup_all():
	val.base_dir = f"DynamicRheo_{val.deform_mode:}_Read_{val.read_udf.split('.')[0]:}"
	if val.sweep_mode == 'StrainSweep':
		for data in val.StrainSweepConditions:
			[val.temperature, val.min_strain, val.max_strain, val.frequency, val.data_per_digit] = data
			val.dir_name = f'{val.deform_mode:}_{val.sweep_mode:}_from_{val.min_strain:.3g}_to_{val.max_strain:.3g}_Freq_{val.frequency:.3g}_Temp_{val.temperature:.3g}'.replace('.', '_')
			val.subdir_list.append(val.dir_name)
			val.target_dir = os.path.join(val.base_dir, val.dir_name)

			val.strain_list = make_series(val.min_strain, val.max_strain, val.data_per_digit)
			val.freq_list = [val.frequency]

			setup_udfs()

	elif val.sweep_mode == 'FrequencySweep':
		for data in val.FrequencySweepConditions:
			[val.temperature, val.min_freq, val.max_freq, val.strain, val.data_per_digit] = data
			val.dir_name = f'{val.deform_mode:}_{val.sweep_mode:}_from_{val.min_freq:.3g}_to_{val.max_freq:.3g}_Strain_{val.strain:.3g}_Temp_{val.temperature:.3g}'.replace('.', '_')
			val.subdir_list.append(val.dir_name)
			val.target_dir = os.path.join(val.base_dir, val.dir_name)

			val.freq_list = make_series(val.min_freq, val.max_freq, val.data_per_digit)
			val.strain_list = [val.strain]

			setup_udfs()
	return



def make_series(min, max, per_digit):
	series = []
	data = min
	while data < max:
		for i in range(per_digit):
			series.append(f'{data*10.**float(i/per_digit):.3g}')
		data *= 10
	series.append(f'{max:.3g}')
	return series


def setup_udfs():
	make_dir()

	summary_fname = 'Result.txt'
	summary = "# analysis of dynamic viscoelasticity\n"
	summary += f"# read Structure from: {val.read_udf}"
	summary += "# file_name\tstrain_amplitude\tfrequency\tfitting_error\tdelta\tsigma_0\tomega\tG'\tG''\ttan_d\n\n"
	with open(os.path.join(val.target_dir, summary_fname), "w") as f:
		f.write(summary)

	udf = UDFManager(val.read_udf)
	udf.eraseRecord(record_pos=0, record_num=udf.totalRecord()-1)
	udf.put(['', -1], 'Initial_Structure.Read_Set_of_Molecules')
	loc = "Initial_Structure.Generate_Method"
	udf.put("Restart", loc + ".Method")
	udf.put(['', -1, 1, 0], loc + ".Restart")
	udf.put(0, "Initial_Structure.Relaxation.Relaxation")
	udf.write(os.path.join(val.target_dir, 'base_uin.udf'))
	
	make_udfs()
	return

def make_dir():
	if os.path.exists(val.target_dir):
		print("Use existing dir of ", val.target_dir)
	else:
		print("Make new dir of ", val.target_dir)
		os.makedirs(val.target_dir)
	return

def make_udfs():

	for i, freq in enumerate(val.freq_list):
		for j, strain in enumerate(val.strain_list):
			val.batch = "#!/bin/bash\n"

			job_id = f'Freq_{freq:}_Strain_{strain:}'.replace('.', '_')
			in_fname  = job_id + '_uin.udf'
			out_fname = job_id + '_out.udf'
			log_fname = job_id + '.log'

			time_setup(freq)

			# print(in_fname, out_fname, log_fname)



			# set_cyclic_dir()
			# make_cycle_batch(id)
			# バッチファイルを作成
			f_batch = os.path.join(val.target_dir, '_calc.bat')
			with open(f_batch, 'w') as f:
				f.write(val.batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
	batch_series()

	# 		bat.write(cognac_v + ' -I ' + in_fname + ' -O ' + out_fname + ' -n ' + str(core) +
	# 			' > ' + log_fname + '\n')
	# 		bat.write('python ' + py_fname + ' ' + out_fname +
	# 			' >> ' + summary_fname + '\n')



	# 		print("%10.5f  %10.5f  %8.5f  %10d  %8d  %s" %(freq, strain, dt, total_step, interval, in_fname))

	# 		create_in_udf(job_dir, base_udf, in_fname, strain, freq, dt, total_step, interval)
	# #
	# shutil.copy(py_fname, job_dir)

	# bat = open(os.path.join(job_dir, '_Run.bat'), "w")




def time_setup(freq):
	freq = round(float(freq), 3)
	# interval と total_step の算出
	# 1) 小数点計算により、dt, step_per_cycle の概算値を決める。
	dt = min(val.dt_max, 1.0/(freq*val.step_min))
	step_per_cycle = 1.0/(freq*dt)
	# 2) 整数値で、intervalを決める。
	interval = max(1, int(round(step_per_cycle/val.output_cycle, -1))) 
	# 3) 整数値 interval を用いて、step_per_cycle と dt を決める。
	dt = round(1.0/(freq*interval*val.output_cycle), 4)
	tmp = 1.0/(interval*dt*freq) # 0.025  400
	N = int(round(tmp))
	if abs(tmp - N) > 1e-5:
		interval = max(1, int(round((step_per_cycle/val.output_cycle)*10., -1)))
		dt = round(1.0/(freq*interval*val.output_cycle), 4)
	# print(1.0/(interval*dt*freq))
	step_per_cycle = interval*val.output_cycle
	total_step = val.total_cycles*step_per_cycle

	val.time = [dt, total_step, interval]
	return








# create input UDF file
def create_in_udf(job_dir, base_udf, in_fname, strain, freq, dt, total_step, interval):

	uobj = UDFManager(os.path.join(job_dir, base_udf))
	loc = "Simulation_Conditions.Dynamics_Conditions.Deformation"
	uobj.put("Lees_Edwards", loc + ".Method")
	uobj.put("Dynamic",    loc + ".Lees_Edwards.Method")
	uobj.put([strain, freq], loc + ".Lees_Edwards.Dynamic")

	loc = "Simulation_Conditions.Dynamics_Conditions.Time"
	uobj.put([dt, total_step, interval], loc)

	uobj.write(os.path.join(job_dir, in_fname))





#######################################################################
# ファイル名を設定し、バッチファイルを作成
def batch_series():
	batch_series = "#!/bin/bash\n"
	for subdir in val.subdir_list:
		if platform.system() == "Windows":
			batch_series += 'cd /d %~dp0\\' + subdir +'\n'
			batch_series += 'call _calc.bat\n'
		elif platform.system() == "Linux":
			batch_series += 'cd ./' + subdir +'\n'
			batch_series += './_calc.bat\n'
			batch_series += 'cd ../\n'
	if platform.system() == "Windows":
		batch_series += 'cd /d %~dp0\n'

	f_batch = os.path.join(val.base_dir, '_calc_all.bat')
	with open(f_batch, 'w') as f:
		f.write(batch_series)
		if platform.system() == "Linux":
			os.chmod(f_batch, 0o777)
	return