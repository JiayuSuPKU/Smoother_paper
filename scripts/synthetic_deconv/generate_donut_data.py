import os
import itertools

def main():
	# The following code generates synthetic spatial data for the donut scenario

	# output directory
	OUT_DIR = '/Users/jysumac/Projects/Smoother_paper/data/synthetic_deconv/donut/'
	# simulation scripts directory (from the smoother package, see the github repo)
	SIM_DIR = '/Users/jysumac/Projects/Smoother/simulation/'

	os.chdir(OUT_DIR)
	n_experiments = 10
	n_spots = [50, 50]
	n_regional_zones = 15
	n_ubiquitous_zones = 0

	# sample spatial abundance patterns
	print("===== generate_st_data.sh: generate spatial patterns (donut) =====")
	exec_cmd = f"python {SIM_DIR}/001b_donut2zone.py " + \
	 		   f"-ne {n_experiments} -ns {n_spots[0]} {n_spots[1]} " + \
			   f"-rz {n_regional_zones} -uz {n_ubiquitous_zones} -o ."
	os.system(exec_cmd)

	# make directory for each scenario and simulate spatial counts
	print("===== generate_st_data.sh: generate synthetic counts =====")
	for nm, sm in itertools.product(['0', '20', '50'], ['0', '0.1', '0.5']):
		dir_name = f"ne10_rz15_nm{nm}_sm{sm}"
		os.makedirs(OUT_DIR + dir_name, exist_ok=True)
		os.chdir(OUT_DIR + dir_name)

		print(f"generate spatial data for {dir_name} ...")
		exec_cmd = f"python {SIM_DIR}/002_zone2counts.py " + \
				   f"-d {SIM_DIR}/example_adata/E-MTAB-11115_subset.h5ad -l annotation_1 " + \
				   f"-i '../zone_abundances.csv' -o {OUT_DIR + dir_name} " + \
				   f"-ne {n_experiments} -ns {n_spots[0]} {n_spots[1]} " + \
				   f"-nm {nm} -sm {sm} "
		if nm == '0': # return full spatial data for benchmarking against other methods
			exec_cmd += '--output-anndata '

		exec_cmd += '>/dev/null 2>&1'
		os.system(exec_cmd)
		os.chdir('..')

	print("===== generate_st_data.sh: finished =====")

if __name__ == '__main__':
	main()
