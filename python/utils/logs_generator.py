import sys, logging

# auxiliary function to initiatite a logger
def init_commandline_logger(logger):
	logger.setLevel(logging.DEBUG)
	# create console handler and set level to debug
	ch = logging.StreamHandler(sys.stdout)
	ch.setLevel(logging.DEBUG)
	# create formatter
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	# add formatter to ch
	ch.setFormatter(formatter)
	# add ch to logger
	logger.addHandler(ch)
