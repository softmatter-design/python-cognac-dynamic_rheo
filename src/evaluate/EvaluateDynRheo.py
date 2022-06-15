#!/usr/bin/env python
######################################################
#    Python script for analysis of dynamic rheology
#    from Cognac Calculations
#                     2007/2/20   H. Kodama, H. Sasaki
#                     2007/9/01   J. Takimoto
#                     2016/4/18   H.Sasaki
######################################################
import numpy
from math import sin, cos, tan, pi, sqrt
import os
import gnuplot2
from time import sleep

import evaluate.values as val
#========== parameters ======================================================
# number of cycles to be omitted from the analysis
skip_cycles = 10
# make this smaller for higher accuracy
epsilon = 0.001
# maximum number of iteration
max_iter = 1000000

#####################
#
def evaluate_all():
	print('kk')
	return


#---- initialize

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
		
