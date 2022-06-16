#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#	動的粘弾性シミュレータ
#                     2007/2/20   H. Kodama, H. Sasaki
#                     2007/9/01   J. Takimoto
#                     2018/4/18   H. Sasaki
################################################################################
# #========== Change Here ===========
# # 計算に使用するコア数
# core = 7
# #
# # 計算条件
# Method = "Freq_Sweep"
# # Method = "Str_Sweep"
# #
# # 基準周波数
# Freq_Base = 0.01
# #
# # 基準ひずみ
# Str_Base = 0.1
# #
# # スイープレンジ
# Range = [1, 0.1]
# #
# # 解像度
# # Resol = [1.0, 0.79, 0.63, 0.5, 0.4, 0.32, 0.25, 0.2, 0.16, 0.13]
# # Resol = [1.0, 0.63, 0.4, 0.25, 0.16]
# Resol = [1., 0.56, 0.32, 0.18]
# # Resol = [1., 0.5, 0.2]
# #
# # シミュレーション時の最大 time step
# dt_max = 0.01
# #
# Setting = {"core":core, "Method":Method, "Freq_Base":Freq_Base,
# "Str_Base":Str_Base, "Range":Range, "Resol":Resol, "dt_max": dt_max}
# ################################################################################
# #========== parameters ==========
# # １サイクル当たりのステップ数の最小値
# # (maximum value of freq*dt) = 1.0/Step_min
# Step_min = 200
# # １サイクル当たりの出力回数（Step_minを割り切れる数字）
# output_per_cycle = 40
# # シミュレーションを行うサイクル数
# num_cycles = 30
# # name of analysis script
# py_fname = "Analysis_V2.py"
# # cognac version
# cognac_v = "cognac10"
# #
# Parameters = {"Step_min": Step_min, "output_per_cycle": output_per_cycle,
# "num_cycles": num_cycles, "py_fname": py_fname, "cognac_v": cognac_v}
# #========== end of parameters ==========

###### Modules #################################################################
from UDFManager import *
import sys
import os
import shutil

import dynamic_rheo_setup.values as val
################################################################################
## MAIN
################################################################################
def setup():
	# 各種条件を読み取り
	read_all()
	# # セットアップ
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
	read_nw_cond()
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

# 計算対象の条件を読み取る
def read_nw_cond():
	if not os.access('target_condition.udf', os.R_OK):
		sys.exit("\n'target_condition.udf' is not exists.")
	else:
		cond_u = UDFManager('target_condition.udf')
		val.func = cond_u.get('TargetCond.NetWork.N_Strands')
		val.nu = cond_u.get('TargetCond.System.Nu')
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
	Dynamics:{
		DeformationMode:select{"Stretch", "Shear"} "変形モードを選択",
		SweepMode:select{"StrainSweep", "FrequencySweep"} "Sweep モードを選択",
			StrainSweep[]:{
				Temperature:float "",
				BaseStrain:float "基準となる歪み",
				Frequency:float "最大ひずみ",
				DataperDigit:float "これは１ステップ計算での伸長度　Res = lambda/1_step"
				}
			FrequencySweep:{
				DeformRate[]:float "これらは変形レートのリスト",
				MaxDeformation:float "最大ひずみ",
				Resolution:float "これは１ステップ計算での伸長度　Res = lambda/1_step"
				}
		} "計算ターゲットの条件を設定"		
	CycleDeformation:{
		CyclicDeform:select{"none", "CyclicStretch", "CyclicShear"} "変形モードを選択",
		CyclicStretch:{
			StretchConditions[]:{
				MaxDeformation:float "最大ひずみ",
				Repeat:int "サイクルの繰り返し数",
				DeformRate[]:float "これらは変形レートのリスト",
				Resolution:float "これは１ステップ計算での伸長度　Res = lambda/1_step"
				}
			}
		CyclicShear:{
			ShearConditions[]:{
				MaxDeformation:float "最大ひずみ",
				Repeat:int "サイクルの繰り返し数",
				DeformRate[]:float "これらは変形レートのリスト",
				Resolution:float "これは１ステップ計算での伸長度　Res = lambda/1_step"
				}
			}
		} "計算ターゲットの条件を設定"
	\end{def}	

	\\begin{data}
	
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
	# Simple Deformation
	val.simple_def_mode  = u.get('SimpleDeformation.DeformMode')
	if val.simple_def_mode == 'Stretch':
		val.sim_rate_list = u.get('SimpleDeformation.Stretch.DeformRate[]')
		val.sim_deform_max = u.get('SimpleDeformation.Stretch.MaxDeformation')
		val.sim_resolution = u.get('SimpleDeformation.Stretch.Resolution')
		val.sim_deform = val.simple_def_mode
	elif val.simple_def_mode == 'Shear':
		val.sim_rate_list = u.get('SimpleDeformation.Shear.DeformRate[]')
		val.sim_deform_max = u.get('SimpleDeformation.Shear.MaxDeformation')
		val.sim_resolution = u.get('SimpleDeformation.Shear.Resolution')
		val.sim_deform = val.simple_def_mode
	elif val.simple_def_mode == 'both':
		val.sim_rate_list = u.get('SimpleDeformation.both.DeformRate[]')
		val.sim_deform_max = u.get('SimpleDeformation.both.MaxDeformation')
		val.sim_resolution = u.get('SimpleDeformation.both.Resolution')
	# Cyclic Deformation
	tmp = []
	val.cyc_deform_max = []
	val.cyc_repeat = []
	val.cyc_ratelist = []
	val.cyc_resolution = []
	val.cyclic_deform = u.get('CycleDeformation.CyclicDeform')
	if val.cyclic_deform == 'CyclicStretch':
		tmp = u.get('CycleDeformation.CyclicStretch.StretchConditions[]')
	elif val.cyclic_deform == 'CyclicShear':
		tmp = u.get('CycleDeformation.CyclicShear.ShearConditions[]')
	for data in tmp:
		val.cyc_deform_max.append(data[0])
		val.cyc_repeat.append(data[1])
		val.cyc_ratelist.append(data[2])
		val.cyc_resolution.append(data[3])
	if val.simple_def_mode == 'none' and val.cyclic_deform == 'none':
		sys.exit('No proper condition is selected.\nBye!')
	return
# 
def init_calc():
	text = "################################################" + "\n"
	text += "Cores used for simulation\t\t" + str(val.core ) + "\n"
	text += "################################################" + "\n"
	if val.simple_def_mode != 'none':
		text += "Deform mode:\t\t\t\t" + str(val.simple_def_mode) + "\n"
		text += "Deform Rate:\t\t" + ', '.join(["{0:4.0e}".format(x) for x in val.sim_rate_list]) + "\n"
		text += "Maximum Strain:\t\t\t\t" + str(val.sim_deform_max) + "\n"
		text += "Resolution:\t\t\t\t" + str(round(val.sim_resolution,4)) + "\n"
		text += "################################################" + "\n"
	if val.cyclic_deform != 'none':
		text += "Deform mode:\t\t\t" + str(val.cyclic_deform) + "\n"
		for i in range(len(val.cyc_deform_max)):
			text += 'Cyclic condition #' + str(i) + '\n'
			text += "\tMaximum Strain:\t\t\t" + str(val.cyc_deform_max[i]) + "\n"
			text += "\tRepeat:\t\t\t\t" + str(val.cyc_repeat[i]) + "\n"
			text += "\tCyclic Deform Rate:\t" + ', '.join(["{0:4.0e}".format(x) for x in val.cyc_ratelist[i]]) + "\n"
			text += "\tResolution:\t\t\t" + str(round(val.cyc_resolution[i], 4)) + "\n"
		text += "################################################" + "\n"
	print(text)
	return










































def Setup(target_udf, base_udf, Setting):
	#
	Method = Setting["Method"]
	Freq_Base = Setting["Freq_Base"]
	Str_Base = Setting["Str_Base"]
	Range = Setting["Range"]
	Resol = Setting["Resol"]
	#
	if Method == "Freq_Sweep":
		job_dir = Method + '_from_' + str(int(100000*Freq_Base*Range[0]*Resol[0])/100000).replace('.', '_') + '_to_' + str(int(100000*Freq_Base*Range[-1]*Resol[-1])/100000).replace('.', '_') + "_Strain_" + str(Str_Base).replace('.', '_')
	elif Method == "Str_Sweep":
		job_dir = Method + '_from_' + str(int(100000*Str_Base*Range[0]*Resol[0])/100000).replace('.', '_') + '_to_' + str(int(100000*Str_Base*Range[-1]*Resol[-1])/100000).replace('.', '_') + "_Freq_" + str(Freq_Base).replace('.', '_')
	summary_fname = 'Result.txt'
	#
	if os.path.exists(job_dir):
		print(u"計算用のディレクトリーが既に存在します。")
		sys.exit(1)
	else:
		os.mkdir(job_dir)
	#
	summary = open(os.path.join(job_dir, summary_fname),"w")
	summary.write("# analysis of dynamic viscoelasticity\n")
	summary.write("# file_name\tstrain_amplitude\tfrequency\tfitting_error\tdelta\tsigma_0\tomega\tG'\tG''\ttan_d\n\n")
	#
	shutil.copy(target_udf, job_dir)
	udf = UDFManager(target_udf)
	rec_len = udf.totalRecord() - 1
	udf.eraseRecord(record_pos=0, record_num=udf.totalRecord()-1)
	udf.put([target_udf, rec_len], 'Initial_Structure.Read_Set_of_Molecules')
	loc = "Initial_Structure.Generate_Method"
	udf.put("Restart", loc + ".Method")
	udf.put([target_udf, rec_len, 1, 0], loc + ".Restart")
	udf.put(0, "Initial_Structure.Relaxation.Relaxation")
	udf.write(os.path.join(job_dir, base_udf))
	# print info
	print("###########################################################")
	print("read Structure from: ", target_udf)
	print("Calculation Directory Created:", job_dir)
	print("summary file created:", summary_fname)
	print("###########################################################")
	# print("script file created: ", script_fname)
	print("\n%10s  %10s  %8s  %10s  %8s  %s" % ("freq", "strain", "dt", "total_step", "interval", "fname"))
	#
	return summary_fname, job_dir
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
