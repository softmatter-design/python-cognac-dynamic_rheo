#!/usr/bin/env python
# -*- coding: utf-8 -*-
###### Modules #################################################################
from unittest import TestProgram
from UDFManager import *

import argparse
import codecs
import os
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
	setup_base()
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
		CalcConditions:{"cognac112",1}
		DeformationMode:{"Shear"}
		Dynamics:{
		"FrequencySweep",
		{
		[{1.0000000,1.e-2,10.000000,1.0000000,5}{1.5000000,1.e-2,10.000000,1.0000000,5}]
		}
		{
		[{1.0000000,1.e-2,1.0000000,1.e-2,5}{1.5000000,1.e-2,1.0000000,1.e-2,5}]
		}
		}
		SymulationParameters:{20,20}
		AnalysisParameters:{10}
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
	val.cycles = u.get('SimulationParameters.Cycles')
	val.output_per_cycle = u.get('SimulationParameters.Output_per_Cycle')
	# AnalysisParameters
	val.skipcycles = u.get('AnalysisParameters.SkipCycles')
	return
# 
def init_calc():
	text = "################################################" + "\n"
	text += "Cores used for simulation\t\t" + str(val.core ) + "\n"
	text += "################################################" + "\n"
	text += "Deform mode:\t\t\t\t" + str(val.deform_mode) + "\n"
	text += "################################################" + "\n"
	if val.sweep_mode == 'StrainSweep':
		text += "Sweep mode:\t\t\t" + str(val.sweep_mode) + "\n"
		for data in val.StrainSweepConditions:
			text += 'temperature:\t\t\t\t' + str(data[0]) + '\n'
			text += "Minimum Strain:\t\t\t\t" + str(data[1]) + "\n"
			text += "Maximum Strain:\t\t\t\t" + str(data[2]) + "\n"
			text += "Frequency:\t\t\t\t" + str(data[3]) + "\n"
			text += "Data per Digit:\t\t\t\t" + str(data[4]) + "\n"
	elif val.sweep_mode == 'FrequencySweep':
		text += "Sweep mode:\t\t\t" + str(val.sweep_mode) + "\n"
		for data in val.FrequencySweepConditions:
			text += 'temperature:\t\t\t\t' + str(data[0]) + '\n'
			text += "Minimum Frequency:\t\t\t" + str(data[1]) + "\n"
			text += "Maximum Frequency:\t\t\t" + str(data[2]) + "\n"
			text += "Strain:\t\t\t\t\t" + str(data[3]) + "\n"
			text += "Data per Digit:\t\t\t\t" + str(data[4]) + "\n"
	text += "################################################" + "\n"
	text += 'SimulationParameters:\n'
	text += "\tCycles:\t\t\t\t" + str(val.cycles) + "\n"
	text += "\tOutput_per_Cycle:\t\t" + str(val.output_per_cycle) + "\n"
	text += 'AnalysisParameters:\n'
	text += "\tSkip Cycles:\t\t\t" + str(val.skipcycles) + "\n"
	text += "################################################" + "\n"
	print(text)
	return






def setup_base():
	if val.sweep_mode == 'StrainSweep':
		mode = val.sweep_mode
		for data in val.StrainSweepConditions:
			val.temperature = data[0]
			val.min_strain = data[1]
			val.max_strain = data[2]
			val.frequency = data[3]
			val.job_dir.append(mode + '_from_'+ '{:.3g}'.format(val.min_strain).replace('.', '_') + '_to_' +  '{:.3g}'.format(val.max_strain).replace('.', '_') + '_Freq_' + '{:.3g}'.format(val.frequency).replace('.', '_') + '_Temp_' + '{:.2g}'.format(val.temperature)).replace('.', '_')
	elif val.sweep_mode == 'FrequencySweep':
		mode = val.sweep_mode
		for data in val.FrequencySweepConditions:
			val.temperature = data[0]
			val.min_freq = data[1]
			val.max_freq = data[2]
			val.strain = data[3]
			val.job_dir.append(mode + '_from_' + '{:.2e}'.format(val.min_freq).replace('.', '_') + '_to_' + '{:.2e}'.format(val.max_freq).replace('.', '_')  + '_Strain_' + '{:.2e}'.format(val.strain).replace('.', '_')  + '_Temp_' + '{:.1e}'.format(val.temperature).replace('.', '_'))
	print(val.job_dir)

	# if Method == "Freq_Sweep":
	# 	job_dir = Method + '_from_' + str(int(100000*Freq_Base*Range[0]*Resol[0])/100000).replace('.', '_') + '_to_' + str(int(100000*Freq_Base*Range[-1]*Resol[-1])/100000).replace('.', '_') + "_Strain_" + str(Str_Base).replace('.', '_')
	# elif Method == "Str_Sweep":
	# 	job_dir = Method + '_from_' + str(int(100000*Str_Base*Range[0]*Resol[0])/100000).replace('.', '_') + '_to_' + str(int(100000*Str_Base*Range[-1]*Resol[-1])/100000).replace('.', '_') + "_Freq_" + str(Freq_Base).replace('.', '_')
	# summary_fname = 'Result.txt'
	# #
	# if os.path.exists(job_dir):
	# 	print(u"計算用のディレクトリーが既に存在します。")
	# 	sys.exit(1)
	# else:
	# 	os.mkdir(job_dir)
	# #
	# summary = open(os.path.join(job_dir, summary_fname),"w")
	# summary.write("# analysis of dynamic viscoelasticity\n")
	# summary.write("# file_name\tstrain_amplitude\tfrequency\tfitting_error\tdelta\tsigma_0\tomega\tG'\tG''\ttan_d\n\n")
	# #
	# shutil.copy(target_udf, job_dir)
	# udf = UDFManager(target_udf)
	# rec_len = udf.totalRecord() - 1
	# udf.eraseRecord(record_pos=0, record_num=udf.totalRecord()-1)
	# udf.put([target_udf, rec_len], 'Initial_Structure.Read_Set_of_Molecules')
	# loc = "Initial_Structure.Generate_Method"
	# udf.put("Restart", loc + ".Method")
	# udf.put([target_udf, rec_len, 1, 0], loc + ".Restart")
	# udf.put(0, "Initial_Structure.Relaxation.Relaxation")
	# udf.write(os.path.join(job_dir, base_udf))
	# # print info
	# print("###########################################################")
	# print("read Structure from: ", target_udf)
	# print("Calculation Directory Created:", job_dir)
	# print("summary file created:", summary_fname)
	# print("###########################################################")
	# # print("script file created: ", script_fname)
	# print("\n%10s  %10s  %8s  %10s  %8s  %s" % ("freq", "strain", "dt", "total_step", "interval", "fname"))
	#
	return
#
def Make_Series_Calc(Setting, Parameters, summary_fname, base_udf, job_dir):
	#
	core = Setting["core"]
	Method = Setting["Method"]
	Freq_Base = Setting["Freq_Base"]
	Str_Base = Setting["Str_Base"]
	Range = Setting["Range"]
	Resol = Setting["Resol"]
	dt_max = Setting["dt_max"]
	#
	Step_min = Parameters["Step_min"]
	output_per_cycle = Parameters["output_per_cycle"]
	num_cycles = Parameters["num_cycles"]
	py_fname = Parameters["py_fname"]
	cognac_v = Parameters["cognac_v"]
	#
	bat = open(os.path.join(job_dir, '_Run.bat'), "w")
	#
	Freq_list = []
	Strain_list = []
	if Method == "Freq_Sweep":
		for ratio in Range:
			for rate in Resol:
				Freq_list.append(Freq_Base*rate*ratio)
		Strain_list.append(Str_Base)
	elif Method == "Str_Sweep":
		for ratio in Range:
			for rate in Resol:
				Strain_list.append(Str_Base*rate*ratio)
		Freq_list.append(Freq_Base)
	#
	for i, freq in enumerate(Freq_list):
		for j, strain in enumerate(Strain_list):
			job_id = Method + '_' + ("{0:6.1e}".format(freq)).replace('.', '_') +'_Strain_'+ ("{0:6.1e}".format(strain)).replace('.', '_')
			in_fname  = job_id + '_uin.udf'
			out_fname = job_id + '_out.udf'
			log_fname = job_id + '.log'

			bat.write(cognac_v + ' -I ' + in_fname + ' -O ' + out_fname + ' -n ' + str(core) +
				' > ' + log_fname + '\n')
			bat.write('python ' + py_fname + ' ' + out_fname +
				' >> ' + summary_fname + '\n')

			# interval と total_step の算出
			# 1) 小数点計算により、dt, step_per_cycle の概算値を決める。
			dt = min(dt_max, 1.0/(freq*Step_min))
			step_per_cycle = 1.0/(freq*dt)	# float
			# 2) 整数値で、intervalを決める。
			interval = max(1, iround(step_per_cycle/output_per_cycle))
			# 3) 整数値 interval を用いて、step_per_cycle と dt を決める。
			dt = 1.0/(freq*interval*output_per_cycle)
			tmp = 1.0/(interval*dt*freq)
			N = int(round(tmp))
			if abs(tmp - N) > 1e-5:
				interval = max(1, iround((step_per_cycle/output_per_cycle)*10.))
				dt = 1.0/(freq*interval*output_per_cycle)
			# print(1.0/(interval*dt*freq))
			step_per_cycle = interval*output_per_cycle
			total_step = num_cycles*step_per_cycle

			print("%10.5f  %10.5f  %8.5f  %10d  %8d  %s" %(freq, strain, dt, total_step, interval, in_fname))

			create_in_udf(job_dir, base_udf, in_fname, strain, freq, dt, total_step, interval)
	#
	shutil.copy(py_fname, job_dir)

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

# integer closest to x
def iround(x):
	return int(round(x))
