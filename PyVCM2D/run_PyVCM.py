import time
from PyVCM2D import PyVCM

case = 'problem_setup.yml'  # Write a for to run several cases in the future

initial_time = time.time()

VCM = PyVCM.PyVCM(case)

final_time = time.time()
print("\nThe simulation lasted {0}s".format(final_time - initial_time))
