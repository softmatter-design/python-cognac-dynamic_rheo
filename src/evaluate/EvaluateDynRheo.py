#!/usr/bin/env python
# -*- coding: utf-8 -*-
###### Modules #################################################################
from UDFManager import *
from scipy.optimize import curve_fit

import numpy as np
import argparse
import os
import platform
import subprocess
import sys

import evaluate.values as val
#####################
# Main
def evaluate_all():
	read_arg()
	read_and_calc()
	fit_save()
	# save_ss()
	plot_ss()
	
	return


#---- initialize

# Read argument 
def read_arg():
	parser = argparse.ArgumentParser(description='Evaluate deformed simulations !')
	parser.add_argument('udf', help="UDF file name to evaluate")
	parser.add_argument('-s', '--skip', help="Skip time to ignore ss data")
	parser.add_argument('-f', '--freq', help="Frequency used for simulation")
	parser.add_argument('-a', '--amp', help="Amplitude of deform used for simulation")

	args = parser.parse_args()
	if args.udf:
		val.read_udf = args.udf
		val.def_mode = args.udf.split('_')[0].lower()
		val.base_name = val.read_udf.split('.')[0]
	else:
		print('\n#####\nUDF file is not specified')
		sys.exit('Set UDF file name!')
	if args.skip:
		val.skip = args.skip
	else:
		print('\n#####\nSkip time is not set!')
		print('all data will be evaluated !')
	if args.freq and args.amp:
		val.freq = args.freq
		val.amp = args.amp
	else:
		print('\n#####\nFrequency and/or Amplitude are not set!')
		sys.exit('Input proper data !')
	return





############################
# Calculate stress either for shear or stretch deformation
def read_and_calc():
	print("Readin file = ", val.read_udf)
	uobj = UDFManager(val.read_udf)
	#
	uobj.jump(0)
	cell = uobj.get("Structure.Unit_Cell.Cell_Size")
	area_init = cell[0]*cell[1]
	z_init = cell[2]
	for i in range(1, uobj.totalRecord()):
		uobj.jump(i)
		print("Reading Rec.=", i)
		time = uobj.get('Time')
		omega_t = round(2*np.pi*float(val.freq)*float(time), 5)
		if val.def_mode == 'shear':
			stress = uobj.get('Statistics_Data.Stress.Total.Batch_Average.xy')
			strain = uobj.get('Structure.Unit_Cell.Shear_Strain')
		elif val.def_mode == 'stretch':
			cell = uobj.get("Structure.Unit_Cell.Cell_Size")
			stress_list = uobj.get("Statistics_Data.Stress.Total.Batch_Average")
			stress = (cell[0]*cell[1])*(stress_list[2]-(stress_list[0] + stress_list[1])/2.)/area_init
			strain = uobj.get("Structure.Unit_Cell.Cell_Size.c")/z_init
		val.ss_data.append([time, omega_t, str(strain), stress])
	return


# ########################################
# # 計算結果をターゲットファイル名で保存
# def save_ss():
# 	val.f_name = val.base_name + '_ss.dat'
# 	with open(val.f_name, 'w') as f:
# 		f.write('# Time\tomega_t\tStrain\tStress\n\n')
# 		for line in val.ss_data:
# 			f.write(f'{line[0]:}\t{line[1]:}\t{line[2]}\t{line[3]}\n')
# 		cnt = 1
# 		while cnt*float(val.skip) < float(val.ss_data[-1][0]):
# 			f.write('\n\n#\n')
# 			f.write('# Time\tomega_t\tStrain\tStress\n\n')
# 			for line in val.ss_data:
# 				if float(line[0]) >= cnt*float(val.skip):
# 					f.write(f'{line[0]:}\t{line[1]:}\t{line[2]}\t{line[3]}\n')
# 			cnt += 1
# 	return


############################
# 結果をプロット
def plot_ss():
	script_content()
	with open(val.base_name + '_plot_ss.plt', 'w') as f:
		f.write(val.script)
	#
	if platform.system() == "Windows":
		subprocess.call([val.base_name + '_plot_ss.plt'], shell=True)
	elif platform.system() == "Linux":
		subprocess.call(['gnuplot ' + val.base_name + '_plot_ss.plt'], shell=True)
	return


# スクリプトの中身
def script_content():
	val.script = 'set term pngcairo font "Arial,14"\n'
	val.script += '#set mono\nset colorsequence classic\n\n'
	val.script += f'data = "{val.f_name}"\n'
	val.script += f'set output "{val.base_name}_lissajous.png"\n\n'
	val.script += 'set key left\nset size square\n'
	val.script += '#set xrange [1:3]\n#set yrange [0.:]\n#set xtics 0.5\n#set ytics 0.01\n'
	val.script += 'set xlabel "Strain"\nset ylabel "Stress"\n\n'
	val.script += 'plot '
	val.script += 'data ind 0 u 3:4 w l lw 2 lt 1 ti "all" ,\\\n'
	val.script += 'data ind 1 u 3:4 w l lw 2 lt 2 ti "skip=1" ,\\\n'
	val.script += 'data ind 2 u 3:4 w l lw 2 lt 3 ti "skip=2" ,\\\n'
	val.script += 'data ind 3 u 3:4 w l lw 2 lt 4 ti "skip=3" ,\\\n'
	val.script += 'data ind 0 u 3:5 w l lw 4 lt 5 ti "Fitted"'
	val.script += '\n\nreset\n\n'
	val.script += 'set term pngcairo font "Arial,14"\n'
	val.script += '#set mono\nset colorsequence classic\n\n'
	val.script += f'data = "{val.f_name}"\n'
	val.script += f'set output "{val.base_name}_time_stress.png"\n\n'
	val.script += 'set key left\nset size square\n'
	val.script += '#set xrange [1:3]\n#set yrange [0.:]\n#set xtics 0.5\n#set ytics 0.01\n'
	val.script += 'set xlabel "{/Symbol w}t"\nset ylabel "Stress"\n\n'
	val.script += 'sigma0=0.1\ndelta=1.\n'
	val.script += 'f(x) = sigma0*sin(x + delta)\n'
	val.script += 'fit f(x) data u 1:3 via sigma0, delta\n\n'
	val.script += 'plot '
	val.script += 'data ind 0 u 2:4 w l lw 2 lt 1 ti "all" ,\\\n'
	val.script += 'data ind 1 u 2:4 w l lw 2 lt 2 ti "skip=1" ,\\\n'
	val.script += 'data ind 2 u 2:4 w l lw 2 lt 3 ti "skip=2" ,\\\n'
	val.script += 'data ind 3 u 2:4 w l lw 2 lt 4 ti "skip=3",\\\n'
	val.script += 'data ind 0 u 2:5 w l lw 1 lt 5 ti "Fitted"'
	val.script += '\n\nreset'
	return




def fit_save():
	omega_t = []
	sigma_we = []
	for data in val.ss_data:
		if float(data[0]) >= float(val.skip):
			omega_t.append(float(data[1]))
			sigma_we.append(data[3])
	g_sigma0 = 0.5
	g_delta = 0.5
	guess = [g_sigma0, g_delta]
	popt, perr = fit_rheo(omega_t, sigma_we, guess)
	#
	f_sigma0, f_delta = popt
	error_sigma0, error_delta = perr
	omega = 2*np.pi*float(val.freq)
	g_prime = f_sigma0/float(val.amp)*np.cos(float(f_delta))
	g_dprime = f_sigma0/float(val.amp)*np.sin(float(f_delta))
	tan_d = np.tan(f_delta)
	with open('Result.txt', 'a') as f:
		f.write(f'{val.amp:}\t{val.freq:}\t{f_sigma0:.4g}\t{error_sigma0:.2e}\t{f_delta:.4g}\t{error_delta:.2e}\t{omega:.4g}\t{g_prime:.4g}\t{g_dprime:.4g}\t{tan_d:.4g}\n')
	#
	val.f_name = val.base_name + '_fitted_ss.dat'
	with open(val.f_name, 'w') as f:
		f.write('# Time\tomega_t\tStrain\tStress\tFitted Stress\n\n')
		for line in val.ss_data:
			f.write(f'{line[0]:}\t{line[1]:}\t{line[2]}\t{line[3]}\t{f_sigma0*np.sin(float(line[1]) + f_delta)}\n')
		cnt = 1
		while cnt*float(val.skip) < float(val.ss_data[-1][0]):
			print(cnt*float(val.skip), float(val.ss_data[-1][0]))
			f.write('\n\n#\n')
			f.write('# Time\tomega_t\tStrain\tStress\tFitted Stress\n\n')
			for line in val.ss_data:
				if float(line[0]) >= cnt*float(val.skip):
					f.write(f'{line[0]:}\t{line[1]:}\t{line[2]}\t{line[3]}\t{f_sigma0*np.sin(float(line[1]) + f_delta)}\n')
			cnt += 1
	return

def func(omega_t, sigma0, delta):
	sigma = sigma0*np.sin(omega_t + delta)
	return sigma

def fit_rheo(omega_t, sigma_we, guess):
	popt, pcov = curve_fit(func, np.array(omega_t), np.array(sigma_we), p0 = guess)
	perr = np.sqrt(np.diag(pcov))
	return popt, perr






















def _initialize(udf, verbose=0):
	global uobj, num_data, dt,interval, amp, freq, N, num_cycles, dphi
	uobj = udf
	# read simulation conditions
	loc = "Simulation_Conditions.Dynamics_Conditions."
	[dt, totalSteps, interval] = uobj.get(loc + "Time")
	[amp, freq] = uobj.get(loc + "Deformation.Lees_Edwards.Dynamic")

	# N = number of output per cycle
	num_data = uobj.totalRecord() - 1 # no stress in record 0
	tmp = 1.0/(interval*dt*freq)
	N = int(round(tmp))
	if abs(tmp - N) > 1e-3:
		print("1/(interval*dt*freq) =",tmp,"should be (almost) integer")
	if (num_data % N != 0):
		print("num_data =",num_data,"should be multipe of N =",N)
		return -1
	num_cycles = int(num_data/N)
	dphi = 2.0*pi/N

	if verbose:
		print("total_output:", uobj.totalRecord() - 1)
		print("output_per_cycle:",N)
		print("num_cycles:",num_cycles)
		print("skip_cycles:",skip_cycles)

	return 0

#---- read stress from UDF

def _read_stress():
	global stress, strain, time

	stress = [0.0]*num_data
	strain = [0.0]*num_data
	time = [0.0]*num_data

	# read stress and strain
	loc = "Statistics_Data.Stress."
	for rec in range(1,uobj.totalRecord()): # no Stress in rec=0
		uobj.jump(rec)
		stress[rec-1] = uobj.get(loc+"Total.Batch_Average.xy")
		strain_i = uobj.get("Structure.Unit_Cell.Shear_Strain")
		if abs(strain_i - amp*sin(dphi*rec)) > 1e-5:
			print("strain mismatch at rec =",rec)
			return -1

	for i in range(num_data):
		time[i] = dt*(i + 0.5 + 0.5/interval)
		strain[i] = amp*sin(dphi*(i + 0.5 + 0.5/interval))

	return 0

def _average_stress():
	global sxy, fit, shift, ss
	# allocate arrays
	sxy = numpy.zeros(N,float)
	fit = numpy.zeros(N,float)

	# calculate average stress sxy[]
	for i in range(skip_cycles*N,num_data):
		sxy[i % N] = sxy[i % N] + stress[i]
	sxy = sxy/(num_cycles - skip_cycles)

	# calculate average and amplitude of sxy
	shift = numpy.add.reduce(sxy)/N
	ss = numpy.dot(sxy-shift,sxy-shift)

#----- calculate sigma_0 and set fit[]

def _set_fit(delta):
	global fit, sigma_0
	for i in range(N):
		fit[i] = 0.
		for j in range(interval):
			fit[i] = fit[i] + sin(dphi*(i+(j+1.)/interval)+delta)
	fit /= interval
	sigma_0 = sqrt(ss/numpy.dot(fit,fit))
	fit = sigma_0*fit + shift

#----- find the best vallue of delta

def _find_delta(verbose=0):
	global delta, err
	# start from delta = 0
	delta = 0.0
	_set_fit(delta)
	err = numpy.dot(sxy-fit,sxy-fit)/numpy.dot(sxy,sxy)

	# increase delta while error decreaes
	prev_err = 0.
	frag_decrease = 0
	num_iter = 0
	while (err < prev_err or frag_decrease == 0) and num_iter < max_iter:
		num_iter += 1
		prev_err = err
		delta = delta + epsilon * err
		_set_fit(delta)

		err = numpy.dot(sxy-fit,sxy-fit)/numpy.dot(sxy,sxy)
		if err > 5.0:
			err = 5.0
		if (err < prev_err) and (err < 1.9):	# XXX why 1.9?
			frag_decrease = 1

	if num_iter >= max_iter:
		print("can't find delta after",num_iter,"iterations")
		print("err = ",err)
		return -1

	if verbose:
		print("\nerr = ",err,"(after",num_iter,"iterations)\n")

	return 0


#---------- public functions ----------

#----- the main function: do everything

def do_fit(udf_file, verbose=1):
	udf = UDFManager(udf_file)
	if _initialize(udf, verbose) != 0:
		return -2
	if _read_stress() != 0:
		return -3
	_average_stress()
	if _find_delta(verbose) != 0:
		return -4
	return 0

#----- output the results

# def output_table():
# 	print("Analysis of dynamic viscoelastisity")
# 	print("gamma_0\t",amp)
# 	print("frequency\t",freq)
# 	print("delta\t",delta)
# 	print("sigma_0\t",sigma_0)
# 	print("omega\t",2*pi*freq)
# 	print("G'\t",sigma_0/amp*cos(delta))
# 	print("G''\t",sigma_0/amp*sin(delta))
# 	print("strain\ts_batch_xy\tfitting")

# 	for i in range(N):
# 		print(strain[i],"\t",sxy[i],"\t",fit[i])


def output_summary():
	print(os.path.basename(uobj.udfFile()), "\t",amp,"\t",freq,"\t",\
		err,"\t",delta,"\t",sigma_0,"\t",2*pi*freq,"\t",\
		sigma_0/amp*cos(delta), "\t", sigma_0/amp*sin(delta), "\t", tan(delta))

def save_data(udf_file):
	#
	if udf_file.split("_")[0] == "Freq":
		title = "Freq_" + ("{0:6.1e}".format(freq)).replace('.', '_')
	elif udf_file.split("_")[0] == "Str":
		title = "Str_" + ("{0:6.1e}".format(amp)).replace('.', '_')
	else:
		print("wrong")
		exit(1)
	sim_dat = title + "_sim.dat"
	rheo_dat = "Rheo.dat"
	#
	if os.path.isfile("summary.dat"):
		with open("summary.dat", 'a') as f:
			f.write(str(os.path.basename(uobj.udfFile())) + "\t" + str(amp) + "\t" + str(freq) 
			+ "\t" + str(err) + "\t" + str(delta) + "\t" + str(sigma_0) + "\t" + str(2*pi*freq)
			+ "\t" + str(sigma_0/amp*cos(delta)) + "\t" + str(sigma_0/amp*sin(delta)) + "\t" 
			+ str(tan(delta)) + "\n")
	else:
		with open("summary.dat", 'w') as f:
			f.write("# file_name\tstrain_amplitude\tfrequency\tfitting_error\tdelta\tsigma_0\tomega\tG'\tG''\ttan_d\n\n")
			f.write(str(os.path.basename(uobj.udfFile())) + "\t" + str(amp) + "\t" + str(freq) 
			+ "\t" + str(err) + "\t" + str(delta) + "\t" + str(sigma_0) + "\t" + str(2*pi*freq)
			+ "\t" + str(sigma_0/amp*cos(delta)) + "\t" + str(sigma_0/amp*sin(delta)) + "\t" 
			+ str(tan(delta)) + "\n")
	#
	if os.path.isfile(rheo_dat):
		with open(rheo_dat, 'a') as f:
			f.write(str(round(amp, 4)) + "\t" + str(round(freq, 6)) + "\t" + str(round(sigma_0/amp*cos(delta), 6)) 
			+ "\t" + str(round(sigma_0/amp*sin(delta), 6)) + "\t" + str(round(tan(delta), 6)) + "\n")
	else:
		with open(rheo_dat, 'w') as f:
			f.write("# strain\tfrequency\tG'\tG''\ttan_d\n\n")
			f.write(str(round(amp, 4)) + "\t" + str(round(freq, 6)) + "\t" + str(round(sigma_0/amp*cos(delta), 6)) 
			+ "\t" + str(round(sigma_0/amp*sin(delta), 6)) + "\t" + str(round(tan(delta), 6)) + "\n")
	#
	with open(sim_dat, 'w') as f:
		f.write("# Strain\tStress\tfit\n\n")
		#
		fit2 = fit.tolist()*num_cycles
		for i, tt in enumerate(time):
			f.write(str(tt) + "\t" + str(strain[i]) + "\t" + str(stress[i]) + "\t" 
			+ str(fit2[i]) + "\n")
		#
		fit3 = fit.tolist()
		fit3.append(fit[0])
		f.write("\n\n\n# Data for lissajous\n# Strain\tfit\n\n")
		for i in range(N+1):
			f.write(str(strain[i]) + "\t" + str(fit3[i]) + "\n")

#----- plot

def plot_stress():
	pre = 'set ytics nomirror; set y2tics; set zeroaxis'
	pre += '\nset ylabel "strain"; set y2label "stress"; set xlabel "t"'
	pre += '\nset title "' + os.path.basename(uobj.udfFile()) + '"'
	data1 = [time,strain,'axis x1y1 w l ti "strain"']
	data2 = [time,stress,'axis x1y2 w lp ti "stress"']
	data3 = [time,fit.tolist()*num_cycles,'axis x1y2 w l ti "fit"']
	gnuplot2.plot(pre,[data1, data2, data3])
	sleep(1)

def plot_lissajous(skip=1):
	pre = 'set xlabel "strain"; set ylabel "stress"; set zeroaxis'
	pre += '\nset title "' + os.path.basename(uobj.udfFile()) + '"'
	if skip:
		cmd = [strain[skip_cycles*N:],
			stress[skip_cycles*N:], 'title "stress" w p']
	else:
		cmd = [strain, stress, 'title "stress" w p']
	fit2 = fit.tolist()
	fit2.append(fit[0])
	gnuplot2.plot(pre, [cmd, [strain[0:N+1], fit2, 'title "fit" w l'] ])

# integer closest to x
def iround(x):
	return int(round(x))

#========== main script =====================================================

if __name__ == '__main__':
	from UDFManager import *
	from sys import *

	if len(argv) == 1:
		print("usage: python", argv[0],"XXX_osc_nn_out.udf >> XXX.txt")
		exit(1)
	elif not os.access(argv[1], os.R_OK):
		print(argv[1],"not exists.")
		exit(1)
	#
	udf_file = argv[1]
	if do_fit(udf_file, 0) == 0:
		output_summary()
		plot_stress()
		plot_lissajous()
		save_data(udf_file)
		
