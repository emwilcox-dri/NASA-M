import numpy as np

def get_roca_ram_sizes():

	area = np.array([
			500., 1000., 2000., 3300., 6000.,
						10000., 20000., 33000., 60000.,
						100000., 200000., 330000., 600000.,
						1000000.])

	class1_num_jan = np.array([
			1.8, 1.1, 1.8, 2.3, 3.8,
					9., 9.2, 7.1, 4.1,
					3.3, 2.1, 1., 0.4,
					.048])

	class2_num_jan = np.array([
			25., 17., 20., 25., 30.,
			31., 13., 5., 2.,
			0.7, 0.055, 0.02, 0.,
			0.])
	
	class3_num_jan = np.array([
			250., 41., 28., 20., 13.,
			12., 5., 2.3, 1.,
			0.65, 0.2, 0.08, 0.014,
			0.])

	class1_num_jul = np.array([
			1., 0.7, 1., 1.3, 2.3,
			5., 7., 6., 4.1,
			3.2, 1.3, 0.6, 0.23,
			0.058])

	class2_num_jul = np.array([
			20., 11., 15., 20., 21.,
			28., 11., 5., 2.1,
			1., 0.12, 0.035, 0.007,
			0.01])

	class3_num_jul = np.array([
			210., 37., 22., 18., 11.,
			10., 5.1, 3.1, 1.4,
			0.8, 0.23, 0.048, 0.01,
			0.])
	
	num_jan = class1_num_jan + class2_num_jan + class3_num_jan
	num_jul = class1_num_jul + class2_num_jul + class3_num_jul
			
	return (area, num_jan, num_jul)
