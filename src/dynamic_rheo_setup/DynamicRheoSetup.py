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
	print('OK')
	# 対象となるファイルの選択
	# target_udf = file_select()
	# base_udf = "Base.udf"
	# # セットアップ
	# summary_fname, job_dir = Setup(target_udf, base_udf, Setting)
	# #
	# Make_Series_Calc(Setting, Parameters, summary_fname, base_udf, job_dir)


################################################################################

################################################################################
# Functions
################################################################################
#----- ファイルを選択
def file_select():
	param = sys.argv
	if len(param) == 1:
		# print("to generate input files from XXX_eq_in.udf and XXX_eq_out.udf ...")
		print("usage: python",param[0],"Honya_out.udf")
		sys.exit(1)
	if not os.access(param[1],os.R_OK):
		print(param[1],"not exists.")
		sys.exit(1)
	return param[1]

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
